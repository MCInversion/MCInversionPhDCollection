#include "TerrainBuilder.h"

#define STB_PERLIN_IMPLEMENTATION
#include "utils/STBPerlin.h"

#include <vector>
#include <cmath>

namespace
{
	/// \brief Utility function to generate a random float between min and max
	float RandomFloat(const float& min, const float& max)
	{
		static std::random_device rd;
		static std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(min, max);
		return dis(gen);
	}

	/// \brief Utility function to compute Perlin noise value
	float PerlinNoise(
		const float& x, const float& y, const float& scale, 
		const size_t& octaves, const float& persistence, const float& lacunarity)
	{
		float amplitude = 1.0f;
		float frequency = 1.0f;
		float noiseHeight = 0.0f;

		for (size_t i = 0; i < octaves; ++i)
		{
			const float sampleX = x * frequency / scale;
			const float sampleY = y * frequency / scale;

			const float perlinValue = stb_perlin_noise3(sampleX, sampleY, 0.0f, 0, 0, 0);
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
        const float radius = m_Settings.SamplingRadius;
        const float radiusSquared = radius * radius;
        const float cellSize = radius / std::sqrt(2.0f);

        // Calculate the bounding box of the boundary polygon
        pmp::BoundingBox bbox(m_Settings.BoundaryLoopPolyline);
        const pmp::Point minPoint = bbox.min();
        const pmp::Point maxPoint = bbox.max();
        const float minX = minPoint[0];
        const float maxX = maxPoint[0];
        const float minY = minPoint[1];
        const float maxY = maxPoint[1];

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
        auto initialPoint = pmp::Point(RandomFloat(minX, maxX), RandomFloat(minY, maxY), 0.0f);
        while (!IsPointInPolygon(initialPoint, m_Settings.BoundaryLoopPolyline))
        {
            initialPoint = pmp::Point(RandomFloat(minX, maxX), RandomFloat(minY, maxY), 0.0f);
        }
        AddPoint(initialPoint);

        while (!activeList.empty())
        {
            const auto randomIndex = static_cast<size_t>(RandomFloat(0.0f, static_cast<float>(activeList.size())));
            pmp::Point point = activeList[randomIndex];
            bool found = false;

            for (size_t k = 0; k < m_Settings.SamplingAttempts; ++k)
            {
                const float angle = RandomFloat(0.0f, 2.0f * 3.14159265359f);
                const float distance = RandomFloat(radius, 2.0f * radius);
                auto newPoint = pmp::Point(point[0] + distance * std::cos(angle), point[1] + distance * std::sin(angle), 0.0f);

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
			float elevation = PerlinNoise(point[0], point[1], m_Settings.NoiseScale, m_Settings.Octaves, m_Settings.Persistence, m_Settings.Lacunarity);
			elevation = m_Settings.MinElevation + (m_Settings.MaxElevation - m_Settings.MinElevation) * (elevation * 0.5f + 0.5f); // Normalize to [MinElevation, MaxElevation]
			point[2] = elevation;
		}
	}
}
