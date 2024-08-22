#include "InscribedManifold.h"

#include "geometry/GridUtil.h"
#include "geometry/GeometryConversionUtils.h"
#include "geometry/GeometryUtil.h"

// ---------------------------------------------------------------------------

std::vector<Circle2D> NaiveInscribedCircleCalculator::Calculate(const InscribedCircleInputData& data)
{
    pmp::BoundingBox2 bbox(data.Points);
    pmp::Point2 bboxCenter = (bbox.min() + bbox.max()) * 0.5f;

    // Find the closest point to the bounding box center
    pmp::Scalar minDist = std::numeric_limits<pmp::Scalar>::max();
    for (const auto& point : data.Points)
    {
        pmp::Scalar dist = norm(bboxCenter - point);
        if (dist < minDist)
        {
            minDist = dist;
        }
    }

    Circle2D circle;
    circle.Center = bboxCenter;
    circle.Radius = minDist;

    return { circle };
}

// ---------------------------------------------------------------------------

namespace
{
    std::vector<Circle2D> FilterOverlappingCircles(const std::vector<Circle2D>& circles)
    {
        std::vector<Circle2D> filteredCircles;
        filteredCircles.reserve(circles.size());

        for (const auto& circle : circles)
        {
            bool intersects = false;

            for (auto& filteredCircle : filteredCircles)
            {
                if (Geometry::CircleIntersectsCircle2D(circle.Center, circle.Radius, filteredCircle.Center, filteredCircle.Radius))
                {
                    intersects = true;

                    // If the current circle has a larger radius, replace the existing one
                    if (circle.Radius > filteredCircle.Radius)
                    {
                        filteredCircle = circle;
                    }
                    break;
                }
            }

            // If it doesn't intersect with any filtered circle, add it to the filtered list
            if (!intersects)
            {
                filteredCircles.push_back(circle);
            }
        }

        filteredCircles.shrink_to_fit();
        return filteredCircles;
    }

    std::vector<Sphere3D> FilterOverlappingSpheres(const std::vector<Sphere3D>& spheres)
    {
        std::vector<Sphere3D> filteredSpheres;
        filteredSpheres.reserve(spheres.size());

        for (const auto& sphere : spheres)
        {
            bool intersects = false;

            for (auto& filteredSphere : filteredSpheres)
            {
                if (Geometry::SphereIntersectsSphere3D(sphere.Center, sphere.Radius, filteredSphere.Center, filteredSphere.Radius))
                {
                    intersects = true;

                    // If the current sphere has a larger radius, replace the existing one
                    if (sphere.Radius > filteredSphere.Radius)
                    {
                        filteredSphere = sphere;
                    }
                    break;
                }
            }

            // If it doesn't intersect with any filtered sphere, add it to the filtered list
            if (!intersects)
            {
                filteredSpheres.push_back(sphere);
            }
        }

        filteredSpheres.shrink_to_fit();
        return filteredSpheres;
    }

} // anonymous namespace

// ...........................................................................

std::vector<Circle2D> DistanceFieldInscribedCircleCalculator::Calculate(const InscribedCircleInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("DistanceFieldInscribedCircleCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Circle2D> circles;
    auto& grid = *data.DistanceField;
    ApplyNarrowGaussianBlur2D(grid);

    const auto& dim = grid.Dimensions();
    const auto& orig = grid.Box().min();
    const float cellSize = grid.CellSize();

    const auto Nx = static_cast<unsigned int>(dim.Nx);
    const auto Ny = static_cast<unsigned int>(dim.Ny);

    // Calculate grid bounds within the bounding box
    const auto bbox = pmp::BoundingBox2(data.Points);
    auto minX = static_cast<unsigned int>((bbox.min()[0] - orig[0]) / cellSize);
    auto maxX = static_cast<unsigned int>((bbox.max()[0] - orig[0]) / cellSize);
    auto minY = static_cast<unsigned int>((bbox.min()[1] - orig[1]) / cellSize);
    auto maxY = static_cast<unsigned int>((bbox.max()[1] - orig[1]) / cellSize);

    minX = std::max(1u, minX);
    maxX = std::min(Nx - 2, maxX);
    minY = std::max(1u, minY);
    maxY = std::min(Ny - 2, maxY);

    // Build the KD-tree
    Geometry::PointCloud2D pointCloud;
    pointCloud.points = data.Points;
    Geometry::PointCloud2DTree kdTree(2, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    // Find local maxima of the distance field
    for (unsigned int iy = minY; iy <= maxY; iy++)
    {
        for (unsigned int ix = minX; ix <= maxX; ix++)
        {
            auto localMaxPt = FindLocalMaximumNearScalarGridCell(grid, ix, iy);
            if (!localMaxPt.has_value())
                continue;

            auto closestPtDistSq = Geometry::GetDistanceToClosestPoint2DSquared(kdTree, *localMaxPt);
            if (!closestPtDistSq.has_value())
                continue;

            circles.push_back({ *localMaxPt, std::sqrt(*closestPtDistSq) });
        }
    }

    return FilterOverlappingCircles(circles);
}

// ------------------------------------------------------------------------------------

namespace
{
    void RecursiveSearchForMaxima(
        std::vector<Circle2D>& circles,
        const Geometry::ScalarGrid2D& grid,
        const Geometry::VectorGrid2D& gridGradient,
        unsigned int minX,
        unsigned int minY,
        unsigned int maxX,
        unsigned int maxY,
        unsigned int radius,
        const std::function<void(const pmp::Point2& localMaxPt)>& onTerminate)
    {
        // Base case: Stop recursion when the radius is 2 or the region is too small
        if (radius < 1)
        {
            return;
        }

        const unsigned int centerX = (minX + maxX) / 2;
        const unsigned int centerY = (minY + maxY) / 2;

        if (radius == 1)
        {
            const auto localMaxPt = FindLocalMaximumNearScalarGridCell(grid, centerX, centerY, radius);
            if (!localMaxPt.has_value())
                return;

            onTerminate(*localMaxPt);
            return;
        }

        // Check if there's a local extreme near the center of this region
        if (IsConvergentOrDivergentNearCell(gridGradient, centerX, centerY, radius))
        {
            // Subdivide the region into four quadrants and search each
            const unsigned int midX = (minX + maxX) / 2;
            const unsigned int midY = (minY + maxY) / 2;
            const auto newRadius = static_cast<unsigned int>(std::floor(static_cast<float>(radius) / 2));

            RecursiveSearchForMaxima(circles, grid, gridGradient, minX, minY, midX, midY, newRadius, onTerminate);
            RecursiveSearchForMaxima(circles, grid, gridGradient, midX, minY, maxX, midY, newRadius, onTerminate);
            RecursiveSearchForMaxima(circles, grid, gridGradient, minX, midY, midX, maxY, newRadius, onTerminate);
            RecursiveSearchForMaxima(circles, grid, gridGradient, midX, midY, maxX, maxY, newRadius, onTerminate);
        }
    }

    void RecursiveSearchForMaxima3D(
        std::vector<Sphere3D>& spheres,
        const Geometry::ScalarGrid& grid,
        const Geometry::VectorGrid& gridGradient,
        unsigned int minX,
        unsigned int minY,
        unsigned int minZ,
        unsigned int maxX,
        unsigned int maxY,
        unsigned int maxZ,
        unsigned int radius,
        const std::function<void(const pmp::Point& localMaxPt)>& onTerminate)
    {
        // Base case: Stop recursion when the radius is less than 1
        if (radius < 1)
        {
            return;
        }

        const unsigned int centerX = (minX + maxX) / 2;
        const unsigned int centerY = (minY + maxY) / 2;
        const unsigned int centerZ = (minZ + maxZ) / 2;

        if (radius == 1)
        {
            const auto localMaxPt = FindLocalMaximumNearScalarGridCell(grid, centerX, centerY, centerZ, radius);
            if (!localMaxPt.has_value())
                return;

            onTerminate(*localMaxPt);
            return;
        }

        // Check if there's a local extreme near the center of this region
        if (IsConvergentOrDivergentNearCell(gridGradient, centerX, centerY, centerZ, radius))
        {
            // Subdivide the region into eight octants and search each
            const unsigned int midX = (minX + maxX) / 2;
            const unsigned int midY = (minY + maxY) / 2;
            const unsigned int midZ = (minZ + maxZ) / 2;
            const auto newRadius = static_cast<unsigned int>(std::floor(static_cast<float>(radius) / 2));

            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, minX, minY, minZ, midX, midY, midZ, newRadius, onTerminate);
            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, midX, minY, minZ, maxX, midY, midZ, newRadius, onTerminate);
            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, minX, midY, minZ, midX, maxY, midZ, newRadius, onTerminate);
            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, midX, midY, minZ, maxX, maxY, midZ, newRadius, onTerminate);

            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, minX, minY, midZ, midX, midY, maxZ, newRadius, onTerminate);
            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, midX, minY, midZ, maxX, midY, maxZ, newRadius, onTerminate);
            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, minX, midY, midZ, midX, maxY, maxZ, newRadius, onTerminate);
            RecursiveSearchForMaxima3D(spheres, grid, gridGradient, midX, midY, midZ, maxX, maxY, maxZ, newRadius, onTerminate);
        }
    }

} // anonymous namespace

// ...........................................................................

std::vector<Circle2D> HierarchicalDistanceFieldInscribedCircleCalculator::Calculate(const InscribedCircleInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("HierarchicalDistanceFieldInscribedCircleCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Circle2D> circles;
    auto& grid = *data.DistanceField;
    ApplyNarrowGaussianBlur2D(grid);
    const auto gridGradient = ComputeGradient(grid);

    const auto& dim = grid.Dimensions();
    const auto& orig = grid.Box().min();
    const float cellSize = grid.CellSize();
    const auto Nx = static_cast<unsigned int>(dim.Nx);
    const auto Ny = static_cast<unsigned int>(dim.Ny);

    // Calculate grid bounds within the bounding box
    const auto bbox = pmp::BoundingBox2(data.Points);
    auto minX = static_cast<unsigned int>((bbox.min()[0] - orig[0]) / cellSize);
    auto maxX = static_cast<unsigned int>((bbox.max()[0] - orig[0]) / cellSize);
    auto minY = static_cast<unsigned int>((bbox.min()[1] - orig[1]) / cellSize);
    auto maxY = static_cast<unsigned int>((bbox.max()[1] - orig[1]) / cellSize);

    minX = std::max(1u, minX);
    maxX = std::min(Nx - 2, maxX);
    minY = std::max(1u, minY);
    maxY = std::min(Ny - 2, maxY);

    // Build the KD-tree
    Geometry::PointCloud2D pointCloud;
    pointCloud.points = data.Points;
    Geometry::PointCloud2DTree kdTree(2, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    // Initial search radius - start with half of the largest dimension
    const unsigned int initialRadius = std::max(Nx, Ny) / 2 - 1;

    // Recursively find local maxima in the grid
    RecursiveSearchForMaxima(circles, grid, gridGradient, minX, minY, maxX, maxY, initialRadius,
        // Termination function:
    [&](const pmp::Point2& localMaxPt) {
        const auto nearestPointDistanceSq = GetDistanceToClosestPoint2DSquared(kdTree, localMaxPt);
        if (!nearestPointDistanceSq.has_value())
            return;
        circles.push_back({ localMaxPt, sqrt(*nearestPointDistanceSq) });
    });

    return FilterOverlappingCircles(circles);
}

// ------------------------------------------------------------------------------------

namespace
{
    /// \brief A very simple index wrapper for a 2D grid particle acting in the "particle swarm optimization" (PSO).
    struct Grid2DParticle
    {
        int ix{ -1 }, iy{ -1 }; // -1 means invalid
    };

    /// \brief A very simple index wrapper for a 3D grid particle acting in the "particle swarm optimization" (PSO).
    struct Grid3DParticle
    {
        int ix{ -1 }, iy{ -1 }, iz{ -1 }; // -1 means invalid
    };

} // anonymous namespace

// ....................................................................................

std::vector<Circle2D> ParticleSwarmDistanceFieldInscribedCircleCalculator::Calculate(const InscribedCircleInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("ParticleSwarmDistanceFieldInscribedCircleCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Circle2D> circles;
    auto& grid = *data.DistanceField;
    ApplyNarrowGaussianBlur2D(grid);
    const auto gridGradient = ComputeGradient(grid);
    constexpr float GRAD_EPSILON = 1e-5f;

    const auto& dim = grid.Dimensions();
    const auto Nx = static_cast<unsigned int>(dim.Nx);
    const auto Ny = static_cast<unsigned int>(dim.Ny);

    // Derive the sample interval from grid dimensions
    constexpr float factor = 0.2f; // TODO: make this value a parameter
    const auto sampleInterval = static_cast<unsigned int>(std::max(1u, static_cast<unsigned int>(Nx * factor)));

    // Build the KD-tree
    Geometry::PointCloud2D pointCloud;
    pointCloud.points = data.Points;
    Geometry::PointCloud2DTree kdTree(2, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    std::vector<Grid2DParticle> particles;

    // A container for grid pts where particles have terminated
    std::unordered_set<unsigned int> processedGridPts{};

    // Initialize the particles
    for (unsigned int iy = sampleInterval; iy < Ny - 1; iy += sampleInterval)
    {
        for (unsigned int ix = sampleInterval; ix < Nx - 1; ix += sampleInterval)
        {
            particles.push_back({ static_cast<int>(ix), static_cast<int>(iy) });
        }
    }

    // Iterate over active particles
    while (!particles.empty())
    {
        std::vector<Grid2DParticle> activeParticles;
        for (const auto& particle : particles)
        {
            int ix = particle.ix;
            int iy = particle.iy;

            // Check if the position has already been covered
            if (processedGridPts.contains(iy * Nx + ix))
                continue;

            processedGridPts.insert(iy * Nx + ix);  // Mark this position as covered

            // Get the gradient at the particle's position
            const pmp::Scalar gradientX = gridGradient.ValuesX()[Nx * iy + ix];
            const pmp::Scalar gradientY = gridGradient.ValuesY()[Nx * iy + ix];

            if (gradientX * gradientX + gradientY * gradientY < GRAD_EPSILON)
                continue; // the particle can't move with a zero gradient

            // Determine the direction of movement based on the gradient
            const int moveX = (gradientX > GRAD_EPSILON) ? 1 : (std::abs(gradientX) < GRAD_EPSILON ? 0 : -1);
            const int moveY = (gradientY > GRAD_EPSILON) ? 1 : (std::abs(gradientY) < GRAD_EPSILON ? 0 : -1);

            // ignore stuck particles
            if (moveX == 0 && moveY == 0)
                continue;

            // Update particle position
            ix = std::clamp(ix + moveX, 1, static_cast<int>(Nx) - 2);
            iy = std::clamp(iy + moveY, 1, static_cast<int>(Ny) - 2);

            // Check if the particle has reached the boundary
            if (ix <= 1 || ix >= static_cast<int>(Nx) - 2 || iy <= 1 || iy >= static_cast<int>(Ny) - 2)
            {
                continue; // Particle has reached the boundary and is discarded
            }

            // Check if the particle's position is a local maximum
            const auto localMaxPt = FindLocalMaximumNearScalarGridCell(grid, ix, iy, 1);
            if (localMaxPt.has_value())
            {
                const auto nearestPointDistanceSq = GetDistanceToClosestPoint2DSquared(kdTree, *localMaxPt);
                if (nearestPointDistanceSq.has_value())
                {
                    circles.push_back({ *localMaxPt, sqrt(*nearestPointDistanceSq) });
                }
                continue; // Particle stops if it finds a local maximum
            }

            // Otherwise, keep the particle active for the next iteration
            activeParticles.push_back({ ix, iy });
        }

        particles.swap(activeParticles); // Update the active particles for the next iteration
    }

    return FilterOverlappingCircles(circles);
}

std::vector<Sphere3D> NaiveInscribedSphereCalculator::Calculate(const InscribedSphereInputData& data)
{
    pmp::BoundingBox bbox(data.Points);
    pmp::Point bboxCenter = (bbox.min() + bbox.max()) * 0.5f;

    // Find the closest point to the bounding box center
    pmp::Scalar minDist = std::numeric_limits<pmp::Scalar>::max();
    for (const auto& point : data.Points)
    {
        pmp::Scalar dist = norm(bboxCenter - point);
        if (dist < minDist)
        {
            minDist = dist;
        }
    }

    Sphere3D sphere;
    sphere.Center = bboxCenter;
    sphere.Radius = minDist;

    return { sphere };
}

std::vector<Sphere3D> DistanceFieldInscribedSphereCalculator::Calculate(const InscribedSphereInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("DistanceFieldInscribedSphereCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Sphere3D> spheres;
    auto& grid = *data.DistanceField;
    ApplyNarrowGaussianBlur(grid);

    const auto& dim = grid.Dimensions();
    const auto& orig = grid.Box().min();
    const float cellSize = grid.CellSize();

    const auto Nx = static_cast<unsigned int>(dim.Nx);
    const auto Ny = static_cast<unsigned int>(dim.Ny);
    const auto Nz = static_cast<unsigned int>(dim.Nz);

    // Calculate grid bounds within the bounding box
    const auto bbox = pmp::BoundingBox(data.Points);
    auto minX = static_cast<unsigned int>((bbox.min()[0] - orig[0]) / cellSize);
    auto maxX = static_cast<unsigned int>((bbox.max()[0] - orig[0]) / cellSize);
    auto minY = static_cast<unsigned int>((bbox.min()[1] - orig[1]) / cellSize);
    auto maxY = static_cast<unsigned int>((bbox.max()[1] - orig[1]) / cellSize);
    auto minZ = static_cast<unsigned int>((bbox.min()[2] - orig[2]) / cellSize);
    auto maxZ = static_cast<unsigned int>((bbox.max()[2] - orig[2]) / cellSize);

    minX = std::max(1u, minX);
    maxX = std::min(Nx - 2, maxX);
    minY = std::max(1u, minY);
    maxY = std::min(Ny - 2, maxY);
    minZ = std::max(1u, minZ);
    maxZ = std::min(Nz - 2, maxZ);

    // Build the KD-tree
    Geometry::PointCloud3D pointCloud;
    pointCloud.points = data.Points;
    Geometry::PointCloud3DTree kdTree(3, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    // Find local maxima of the distance field
    for (unsigned int iz = minZ; iz <= maxZ; iz++)
    {
        for (unsigned int iy = minY; iy <= maxY; iy++)
        {
            for (unsigned int ix = minX; ix <= maxX; ix++)
            {
                auto localMaxPt = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz);
                if (!localMaxPt.has_value())
                    continue;

                auto closestPtDistSq = Geometry::GetDistanceToClosestPoint3DSquared(kdTree, *localMaxPt);
                if (!closestPtDistSq.has_value())
                    continue;

                spheres.push_back({ *localMaxPt, std::sqrt(*closestPtDistSq) });
            }
        }
    }

    return FilterOverlappingSpheres(spheres);
}

std::vector<Sphere3D> HierarchicalDistanceFieldInscribedSphereCalculator::Calculate(const InscribedSphereInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("HierarchicalDistanceFieldInscribedSphereCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Sphere3D> spheres;
    auto& grid = *data.DistanceField;
    ApplyNarrowGaussianBlur(grid);
    const auto gridGradient = ComputeGradient(grid);

    const auto& dim = grid.Dimensions();
    const auto& orig = grid.Box().min();
    const float cellSize = grid.CellSize();
    const auto Nx = static_cast<unsigned int>(dim.Nx);
    const auto Ny = static_cast<unsigned int>(dim.Ny);
    const auto Nz = static_cast<unsigned int>(dim.Nz);

    // Calculate grid bounds within the bounding box
    const auto bbox = pmp::BoundingBox(data.Points);
    auto minX = static_cast<unsigned int>((bbox.min()[0] - orig[0]) / cellSize);
    auto maxX = static_cast<unsigned int>((bbox.max()[0] - orig[0]) / cellSize);
    auto minY = static_cast<unsigned int>((bbox.min()[1] - orig[1]) / cellSize);
    auto maxY = static_cast<unsigned int>((bbox.max()[1] - orig[1]) / cellSize);
    auto minZ = static_cast<unsigned int>((bbox.min()[2] - orig[2]) / cellSize);
    auto maxZ = static_cast<unsigned int>((bbox.max()[2] - orig[2]) / cellSize);

    minX = std::max(1u, minX);
    maxX = std::min(Nx - 2, maxX);
    minY = std::max(1u, minY);
    maxY = std::min(Ny - 2, maxY);
    minZ = std::max(1u, minZ);
    maxZ = std::min(Nz - 2, maxZ);

    // Build the KD-tree
    Geometry::PointCloud3D pointCloud;
    pointCloud.points = data.Points;
    Geometry::PointCloud3DTree kdTree(3, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    // Initial search radius - start with half of the largest dimension
    const unsigned int initialRadius = std::max({ Nx, Ny, Nz }) / 2 - 1;

    // Recursively find local maxima in the grid
    RecursiveSearchForMaxima3D(spheres, grid, gridGradient, minX, minY, minZ, maxX, maxY, maxZ, initialRadius,
        // Termination function:
        [&](const pmp::Point& localMaxPt) {
            const auto nearestPointDistanceSq = GetDistanceToClosestPoint3DSquared(kdTree, localMaxPt);
            if (!nearestPointDistanceSq.has_value())
                return;
            spheres.push_back({ localMaxPt, sqrt(*nearestPointDistanceSq) });
        });

    return FilterOverlappingSpheres(spheres);
}

std::vector<Sphere3D> ParticleSwarmDistanceFieldInscribedSphereCalculator::Calculate(const InscribedSphereInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("ParticleSwarmDistanceFieldInscribedSphereCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Sphere3D> spheres;
    auto& grid = *data.DistanceField;
    ApplyNarrowGaussianBlur(grid);
    const auto gridGradient = ComputeGradient(grid);
    constexpr float GRAD_EPSILON = 1e-5f;

    const auto& dim = grid.Dimensions();
    const auto Nx = static_cast<unsigned int>(dim.Nx);
    const auto Ny = static_cast<unsigned int>(dim.Ny);
    const auto Nz = static_cast<unsigned int>(dim.Nz);

    // Derive the sample interval from grid dimensions
    constexpr float factor = 0.2f; // TODO: make this value a parameter
    const auto sampleInterval = static_cast<unsigned int>(std::max(1u, static_cast<unsigned int>(Nx * factor)));

    // Build the KD-tree
    Geometry::PointCloud3D pointCloud;
    pointCloud.points = data.Points;
    Geometry::PointCloud3DTree kdTree(3, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    std::vector<Grid3DParticle> particles;

    // A container for grid points where particles have terminated
    std::unordered_set<unsigned int> processedGridPts{};

    // Initialize the particles
    for (unsigned int iz = sampleInterval; iz < Nz - 1; iz += sampleInterval)
    {
        for (unsigned int iy = sampleInterval; iy < Ny - 1; iy += sampleInterval)
        {
            for (unsigned int ix = sampleInterval; ix < Nx - 1; ix += sampleInterval)
            {
                particles.push_back({ static_cast<int>(ix), static_cast<int>(iy), static_cast<int>(iz) });
            }
        }
    }

    // Iterate over active particles
    while (!particles.empty())
    {
        std::vector<Grid3DParticle> activeParticles;
        for (const auto& particle : particles)
        {
            int ix = particle.ix;
            int iy = particle.iy;
            int iz = particle.iz;

            // Check if the position has already been covered
            if (processedGridPts.contains(iz * Nx * Ny + iy * Nx + ix))
                continue;

            processedGridPts.insert(iz * Nx * Ny + iy * Nx + ix);  // Mark this position as covered

            // Get the gradient at the particle's position
            const pmp::Scalar gradientX = gridGradient.ValuesX()[Nx * Ny * iz + Nx * iy + ix];
            const pmp::Scalar gradientY = gridGradient.ValuesY()[Nx * Ny * iz + Nx * iy + ix];
            const pmp::Scalar gradientZ = gridGradient.ValuesZ()[Nx * Ny * iz + Nx * iy + ix];

            if (gradientX * gradientX + gradientY * gradientY + gradientZ * gradientZ < GRAD_EPSILON)
                continue; // the particle can't move with a zero gradient

            // Determine the direction of movement based on the gradient
            const int moveX = (gradientX > GRAD_EPSILON) ? 1 : (std::abs(gradientX) < GRAD_EPSILON ? 0 : -1);
            const int moveY = (gradientY > GRAD_EPSILON) ? 1 : (std::abs(gradientY) < GRAD_EPSILON ? 0 : -1);
            const int moveZ = (gradientZ > GRAD_EPSILON) ? 1 : (std::abs(gradientZ) < GRAD_EPSILON ? 0 : -1);

            // Ignore stuck particles
            if (moveX == 0 && moveY == 0 && moveZ == 0)
                continue;

            // Update particle position
            ix = std::clamp(ix + moveX, 1, static_cast<int>(Nx) - 2);
            iy = std::clamp(iy + moveY, 1, static_cast<int>(Ny) - 2);
            iz = std::clamp(iz + moveZ, 1, static_cast<int>(Nz) - 2);

            // Check if the particle has reached the boundary
            if (ix <= 1 || ix >= static_cast<int>(Nx) - 2 || iy <= 1 || iy >= static_cast<int>(Ny) - 2 || iz <= 1 || iz >= static_cast<int>(Nz) - 2)
            {
                continue; // Particle has reached the boundary and is discarded
            }

            // Check if the particle's position is a local maximum
            const auto localMaxPt = FindLocalMaximumNearScalarGridCell(grid, ix, iy, iz, 1);
            if (localMaxPt.has_value())
            {
                const auto nearestPointDistanceSq = GetDistanceToClosestPoint3DSquared(kdTree, *localMaxPt);
                if (nearestPointDistanceSq.has_value())
                {
                    spheres.push_back({ *localMaxPt, sqrt(*nearestPointDistanceSq) });
                }
                continue; // Particle stops if it finds a local maximum
            }

            // Otherwise, keep the particle active for the next iteration
            activeParticles.push_back({ ix, iy, iz });
        }

        particles.swap(activeParticles); // Update the active particles for the next iteration
    }

    return FilterOverlappingSpheres(spheres);
}
