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

        for (const auto& circle : circles)
        {
            bool intersects = false;

            // Check if the current circle intersects with any circle in the filtered list
            for (const auto& filteredCircle : filteredCircles)
            {
                if (Geometry::CircleIntersectsCircle2D(circle.Center, circle.Radius, filteredCircle.Center, filteredCircle.Radius))
                {
                    intersects = true;
                    break;
                }
            }

            // If it doesn't intersect with any filtered circle, add it to the filtered list
            if (!intersects)
            {
                filteredCircles.push_back(circle);
            }
        }

        return filteredCircles;
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
    /// \brief A very simple index warpper for a 2D grid particle acting in the "particle swarm optimization" (PSO).
    struct Grid2DParticle
    {
        int ix, iy; // -1 means invalid
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

    // Initialize the covered grid with false values
    std::vector coveredGrid(Nx, std::vector(Ny, false));

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
            if (coveredGrid[ix][iy])
                continue;

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
                    coveredGrid[ix][iy] = true;  // Mark this position as covered
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
