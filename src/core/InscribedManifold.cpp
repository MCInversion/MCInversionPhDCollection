#include "InscribedManifold.h"

#include "geometry/GridUtil.h"
#include "geometry/GeometryConversionUtils.h"

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

std::vector<Circle2D> DistanceFieldInscribedCircleCalculator::Calculate(const InscribedCircleInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("DistanceFieldInscribedCircleCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Circle2D> circles;
    auto grid = *data.DistanceField;
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

    return circles;
}

// ------------------------------------------------------------------------------------

//namespace
//{
    void RecursiveSearchForMaxima(
        std::vector<Circle2D>& circles,
        const Geometry::ScalarGrid2D& grid,
        unsigned int minX,
        unsigned int minY,
        unsigned int maxX,
        unsigned int maxY,
        unsigned int radius,
        const std::function<void(const pmp::Point2& localMaxPt)>& onTerminate)
    {
        // Base case: Stop recursion when the radius is 1 or the region is too small
        if (radius < 1 || maxX - minX < 2 || maxY - minY < 2)
        {
            return;
        }

        const unsigned int centerX = (minX + maxX) / 2;
        const unsigned int centerY = (minY + maxY) / 2;

        // Check if there's a local extreme near the center of this region
        if (ContainsLocalExtremesNearScalarGridCell(grid, centerX, centerY, radius))
        {
            if (radius == 1)
            {
                const auto localMaxPt = FindLocalMaximumNearScalarGridCell(grid, centerX, centerY, radius);
                if (!localMaxPt.has_value())
                    return;

                onTerminate(*localMaxPt);
            }

            // Subdivide the region into four quadrants and search each
            const unsigned int midX = (minX + maxX) / 2;
            const unsigned int midY = (minY + maxY) / 2;
            const unsigned int newRadius = radius / 2 + /* slightly inflate */ 1;

            RecursiveSearchForMaxima(circles, grid, minX, minY, midX, midY, newRadius, onTerminate);
            RecursiveSearchForMaxima(circles, grid, midX, minY, maxX, midY, newRadius, onTerminate);
            RecursiveSearchForMaxima(circles, grid, minX, midY, midX, maxY, newRadius, onTerminate);
            RecursiveSearchForMaxima(circles, grid, midX, midY, maxX, maxY, newRadius, onTerminate);
        }
    }

//} // anonymous namespace

std::vector<Circle2D> HierarchicalDistanceFieldInscribedCircleCalculator::Calculate(const InscribedCircleInputData& data)
{
    if (!data.DistanceField)
    {
        throw std::invalid_argument("HierarchicalDistanceFieldInscribedCircleCalculator::Calculate: Distance field is not provided!");
    }

    std::vector<Circle2D> circles;
    const auto& grid = *data.DistanceField;
    const auto Nx = static_cast<unsigned int>(grid.Dimensions().Nx);
    const auto Ny = static_cast<unsigned int>(grid.Dimensions().Ny);

    // Build the KD-tree
    Geometry::PointCloud2D pointCloud;
    pointCloud.points = data.Points;
    Geometry::PointCloud2DTree kdTree(2, pointCloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    kdTree.buildIndex();

    // Initial search radius - start with half of the largest dimension
    const unsigned int initialRadius = std::max(Nx, Ny) / 2;

    // Recursively find local maxima in the grid
    RecursiveSearchForMaxima(circles, grid, 0, 0, Nx, Ny, initialRadius, [&](const pmp::Point2& localMaxPt) {
        const auto nearestPointDistanceSq = Geometry::GetDistanceToClosestPoint2DSquared(kdTree, localMaxPt);
        if (!nearestPointDistanceSq.has_value())
            return;
        circles.push_back({ localMaxPt, sqrt(*nearestPointDistanceSq) });
     });

    return circles;
}
