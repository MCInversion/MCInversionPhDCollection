#include "TerrainBuilder.h"

#include "Utils/STBPerlin.h"

#include <random>
#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

namespace
{
	/// \brief Utility function to generate a random float between min and max
	float RandomFloat(float min, float max)
	{
		static std::random_device rd;
		static std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(min, max);
		return dis(gen);
	}

	/// \brief Utility function to compute Perlin noise value
	float PerlinNoise(float x, float y, float scale, int octaves, float persistence, float lacunarity)
	{
		float amplitude = 1.0f;
		float frequency = 1.0f;
		float noiseHeight = 0.0f;

		for (int i = 0; i < octaves; ++i)
		{
			float sampleX = x * frequency / scale;
			float sampleY = y * frequency / scale;

			float perlinValue = static_cast<float>(stb_perlin_noise3(sampleX, sampleY, 0.0f, 0, 0, 0));
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
		float radius = m_Settings.SamplingRadius;
		float radiusSquared = radius * radius;
		float cellSize = radius / std::sqrt(2.0f);
		int gridWidth = static_cast<int>((1.0f / cellSize)) + 1;
		int gridHeight = static_cast<int>((1.0f / cellSize)) + 1;

		std::vector grid(gridWidth, std::vector(gridHeight, pmp::Point(-1, -1, 0)));

		auto AddPoint = [&](const pmp::Point& point)
		{
			samples.push_back(point);
			activeList.push_back(point);
			int gridX = static_cast<int>(point[0] / cellSize);
			int gridY = static_cast<int>(point[1] / cellSize);
			grid[gridX][gridY] = point;
		};

		// Start with a random point
		pmp::Point initialPoint = pmp::Point(RandomFloat(0.0f, 1.0f), RandomFloat(0.0f, 1.0f), 0.0f);
		while (!IsPointInPolygon(initialPoint, m_Settings.BoundaryLoopPolyline))
		{
			initialPoint = pmp::Point(RandomFloat(0.0f, 1.0f), RandomFloat(0.0f, 1.0f), 0.0f);
		}
		AddPoint(initialPoint);

		while (!activeList.empty())
		{
			size_t randomIndex = static_cast<size_t>(RandomFloat(0.0f, static_cast<float>(activeList.size())));
			pmp::Point point = activeList[randomIndex];
			bool found = false;

			for (size_t k = 0; k < m_Settings.SamplingAttempts; ++k)
			{
				float angle = RandomFloat(0.0f, 2.0f * 3.14159265359f);
				float distance = RandomFloat(radius, 2.0f * radius);
				pmp::Point newPoint = pmp::Point(point[0] + distance * std::cos(angle), point[1] + distance * std::sin(angle), 0.0f);

				if (!IsPointInPolygon(newPoint, m_Settings.BoundaryLoopPolyline))
				{
					continue;
				}

				int gridX = static_cast<int>(newPoint[0] / cellSize);
				int gridY = static_cast<int>(newPoint[1] / cellSize);

				bool valid = true;
				for (int i = std::max(0, gridX - 2); i <= std::min(gridWidth - 1, gridX + 2) && valid; ++i)
				{
					for (int j = std::max(0, gridY - 2); j <= std::min(gridHeight - 1, gridY + 2) && valid; ++j)
					{
						pmp::Point neighbor = grid[i][j];
						if (neighbor[0] != -1 && sqrnorm(newPoint - neighbor) < radiusSquared)
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
