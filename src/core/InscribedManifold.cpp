#include "InscribedManifold.h"

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
    auto& grid = *data.DistanceField;
    auto& values = grid.Values();
    const auto& dim = grid.Dimensions();
    const auto& orig = grid.Box().min();
    const float cellSize = grid.CellSize();

    const auto Nx = static_cast<unsigned int>(dim.Nx);
    const auto Ny = static_cast<unsigned int>(dim.Ny);

    // Calculate grid bounds within the bounding box
    const auto bbox = pmp::BoundingBox2(data.Points);
    unsigned int minX = static_cast<unsigned int>((bbox.min()[0] - orig[0]) / cellSize);
    unsigned int maxX = static_cast<unsigned int>((bbox.max()[0] - orig[0]) / cellSize);
    unsigned int minY = static_cast<unsigned int>((bbox.min()[1] - orig[1]) / cellSize);
    unsigned int maxY = static_cast<unsigned int>((bbox.max()[1] - orig[1]) / cellSize);

    minX = std::max(1u, minX);
    maxX = std::min(Nx - 2, maxX);
    minY = std::max(1u, minY);
    maxY = std::min(Ny - 2, maxY);

    // Find local maxima of the distance field
    for (unsigned int iy = minY; iy <= maxY; iy++)
    {
        for (unsigned int ix = minX; ix <= maxX; ix++)
        {
            const unsigned int gridPos = Nx * iy + ix;
            pmp::Scalar value = values[gridPos];

            bool isLocalMaximum = true;
            for (int di = -1; di <= 1; ++di)
            {
                for (int dj = -1; dj <= 1; ++dj)
                {
                    if (di == 0 && dj == 0) continue;

                    const unsigned int gridNeighborPos = Nx * (iy + dj) + ix + di;
                    if (values[gridNeighborPos] > value)
                    {
                        isLocalMaximum = false;
                        break;
                    }
                }
                if (!isLocalMaximum) break;
            }

            if (isLocalMaximum)
            {
                Circle2D circle;
                circle.Center = pmp::Point2{
                    orig[0] + static_cast<float>(ix) * cellSize,
                    orig[1] + static_cast<float>(iy) * cellSize
                };
                circle.Radius = value;
                circles.push_back(circle);
            }
        }
    }

    return circles;
}
