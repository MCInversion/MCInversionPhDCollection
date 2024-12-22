#include "TerrainBuilder.h"

#define STB_PERLIN_IMPLEMENTATION
#include "utils/STBPerlin.h"

#include <vector>
#include <cmath>

namespace
{
	/// \brief Utility function to generate a random float between min and max
    pmp::Scalar RandomFloat(const pmp::Scalar& min, const pmp::Scalar& max)
	{
		static std::random_device rd;
		static std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(min, max);
		return dis(gen);
	}

	/// \brief Utility function to compute Perlin noise value
    pmp::Scalar PerlinNoise(
		const pmp::Scalar& x, const pmp::Scalar& y, const pmp::Scalar& scale,
		const size_t& octaves, const pmp::Scalar& persistence, const pmp::Scalar& lacunarity)
	{
        pmp::Scalar amplitude = 1.0;
        pmp::Scalar frequency = 1.0;
        pmp::Scalar noiseHeight = 0.0;

		for (size_t i = 0; i < octaves; ++i)
		{
			const pmp::Scalar sampleX = x * frequency / scale;
			const pmp::Scalar sampleY = y * frequency / scale;

			const pmp::Scalar perlinValue = stb_perlin_noise3(sampleX, sampleY, 0.0, 0, 0, 0);
			noiseHeight += perlinValue * amplitude;

			amplitude *= persistence;
			frequency *= lacunarity;
		}

		return noiseHeight;
	}

	/// \brief Utility function to determine if a point is within the boundary polygon
	bool IsPointInPolygon(const pmp::Point& point, const std::vector<pmp::Point>& polygon)
	{
		// Simple polygon point-inclusion test
		bool inside = false;
		for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++)
		{
			if (((polygon[i][1] > point[1]) != (polygon[j][1] > point[1])) &&
				(point[0] < (polygon[j][0] - polygon[i][0]) * (point[1] - polygon[i][1]) / (polygon[j][1] - polygon[i][1]) + polygon[i][0]))
			{
				inside = !inside;
			}
		}
		return inside;
	}

} // anonymous namespace

namespace Geometry
{

	void TerrainBuilder::PoissonSamplePointsInXYPlane()
	{
        std::vector<pmp::Point> samples;
        std::vector<pmp::Point> activeList;
        const pmp::Scalar radius = m_Settings.SamplingRadius;
        const pmp::Scalar radiusSquared = radius * radius;
        const pmp::Scalar cellSize = radius / std::sqrt(2.0);

        // Calculate the bounding box of the boundary polygon
        pmp::BoundingBox bbox(m_Settings.BoundaryLoopPolyline);
        const pmp::Point minPoint = bbox.min();
        const pmp::Point maxPoint = bbox.max();
        const pmp::Scalar minX = minPoint[0];
        const pmp::Scalar maxX = maxPoint[0];
        const pmp::Scalar minY = minPoint[1];
        const pmp::Scalar maxY = maxPoint[1];

        const int gridWidth = static_cast<int>((maxX - minX) / cellSize) + 1;
        const int gridHeight = static_cast<int>((maxY - minY) / cellSize) + 1;

        std::vector grid(gridWidth, std::vector(gridHeight, pmp::Point(-1, -1, 0)));
        std::vector occupied(gridWidth, std::vector(gridHeight, false));

        auto AddPoint = [&](const pmp::Point& point)
        {
            samples.push_back(point);
            activeList.push_back(point);
            const int gridX = static_cast<int>((point[0] - minX) / cellSize);
            const int gridY = static_cast<int>((point[1] - minY) / cellSize);
            if (gridX >= 0 && gridX < gridWidth && gridY >= 0 && gridY < gridHeight)
            {
                grid[gridX][gridY] = point;
                occupied[gridX][gridY] = true;
            }
        };

        // Start with a random point
        auto initialPoint = pmp::Point(RandomFloat(minX, maxX), RandomFloat(minY, maxY), 0.0);
        while (!IsPointInPolygon(initialPoint, m_Settings.BoundaryLoopPolyline))
        {
            initialPoint = pmp::Point(RandomFloat(minX, maxX), RandomFloat(minY, maxY), 0.0);
        }
        AddPoint(initialPoint);

        while (!activeList.empty())
        {
            const auto randomIndex = static_cast<size_t>(RandomFloat(0.0, static_cast<pmp::Scalar>(activeList.size())));
            pmp::Point point = activeList[randomIndex];
            bool found = false;

            for (size_t k = 0; k < m_Settings.SamplingAttempts; ++k)
            {
                const pmp::Scalar angle = RandomFloat(0.0, 2.0 * static_cast<pmp::Scalar>(M_PI));
                const pmp::Scalar distance = RandomFloat(radius, 2.0 * radius);
                auto newPoint = pmp::Point(point[0] + distance * std::cos(angle), point[1] + distance * std::sin(angle), 0.0);

                if (!IsPointInPolygon(newPoint, m_Settings.BoundaryLoopPolyline))
                {
                    continue;
                }

                const int gridX = static_cast<int>((newPoint[0] - minX) / cellSize);
                const int gridY = static_cast<int>((newPoint[1] - minY) / cellSize);

                bool valid = true;
                for (int i = std::max(0, gridX - 2); i <= std::min(gridWidth - 1, gridX + 2) && valid; ++i)
                {
                    for (int j = std::max(0, gridY - 2); j <= std::min(gridHeight - 1, gridY + 2) && valid; ++j)
                    {
                        if (occupied[i][j] && sqrnorm(newPoint - grid[i][j]) < radiusSquared)
                        {
                            valid = false;
                        }
                    }
                }

                if (valid)
                {
                    AddPoint(newPoint);
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                activeList.erase(activeList.begin() + randomIndex);
            }
        }

        m_Result.Vertices = std::move(samples);
	}

	void TerrainBuilder::GeneratePerlinZElevation()
	{
		for (auto& point : m_Result.Vertices)
		{
            pmp::Scalar elevation = PerlinNoise(point[0], point[1], m_Settings.NoiseScale, m_Settings.Octaves, m_Settings.Persistence, m_Settings.Lacunarity);
			elevation = m_Settings.MinElevation + (m_Settings.MaxElevation - m_Settings.MinElevation) * (elevation * 0.5 + 0.5); // Normalize to [MinElevation, MaxElevation]
			point[2] = elevation;
		}
	}
}
