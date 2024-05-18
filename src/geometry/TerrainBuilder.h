#pragma once

#include "GeometryConversionUtils.h"

namespace Geometry
{
	/**
	 * \brief A wrapper for the settings of the terrain builder.
	 * \struct TerrainSettings
	 */
	struct TerrainSettings
	{
		std::vector<pmp::Point> BoundaryLoopPolyline{}; //>! points determining the boundary of the generated terrain
		pmp::Scalar MinElevation{ 0.0f }; //>! minimum elevation.
		pmp::Scalar MaxElevation{ 100.0f }; //>! maximum elevation.
		
		// Poisson Disk Sampling settings
		pmp::Scalar SamplingRadius{ 1.0f }; //>! Example value, adjust as needed
		size_t SamplingAttempts{ 30 }; //>! Typical value, can be adjusted

		// Perlin Noise settings
		pmp::Scalar NoiseScale{ 1.0f }; //>! Scale of the noise, affects frequency
		size_t Octaves{ 4 }; //>! Number of noise layers
		pmp::Scalar Persistence{ 0.5f }; //>! Amplitude decrease
		pmp::Scalar Lacunarity{ 2.0f }; //>! Frequency increase
	};

	/// \brief Function to verify the validity of TerrainSettings
	inline [[nodiscard]] bool TerrainBuilderSettingsValid(const TerrainSettings& settings)
	{
		// Ensure MinElevation is less than or equal to MaxElevation
		if (settings.MinElevation > settings.MaxElevation)
		{
			std::cerr << "MinElevation cannot be greater than MaxElevation.\n";
			return false;
		}

		// Ensure SamplingRadius is positive
		if (settings.SamplingRadius <= 0)
		{
			std::cerr << "SamplingRadius must be positive.\n";
			return false;
		}

		// Ensure SamplingAttempts is positive
		if (settings.SamplingAttempts == 0)
		{
			std::cerr << "SamplingAttempts must be greater than zero.\n";
			return false;
		}

		// Ensure NoiseScale is positive
		if (settings.NoiseScale <= 0)
		{
			std::cerr << "NoiseScale must be positive.\n";
			return false;
		}

		// Ensure Octaves is a reasonable number (typically between 1 and 10)
		if (settings.Octaves < 1 || settings.Octaves > 10)
		{
			std::cerr << "Octaves must be between 1 and 10.\n";
			return false;
		}

		// Ensure Persistence is between 0 and 1
		if (settings.Persistence < 0 || settings.Persistence > 1)
		{
			std::cerr << "Persistence must be between 0 and 1.\n";
			return false;
		}

		// Ensure Lacunarity is positive
		if (settings.Lacunarity <= 0)
		{
			std::cerr << "Lacunarity must be positive.\n";
			return false;
		}

		return true;
	}

	/// \brief A changeable utility for triangulating a set of points with unique z-coordinates
	using TriangulationFunction = std::function<void(BaseMeshGeometryData&)>;

	/**
	 * \brief A builder for generating terrain mesh as BaseMeshGeometryData using Perlin noise and Poisson disk sampling.
	 * \class TerrainBuilder
	 */
	class TerrainBuilder
	{
	public:
		/// \brief Constructor
		explicit TerrainBuilder(TerrainSettings settings)
			: m_Settings(std::move(settings))
		{
			if (!TerrainBuilderSettingsValid(m_Settings))
			{
				throw std::invalid_argument("Geometry::TerrainBuilder::TerrainBuilder: TerrainBuilderSettingsValid(m_Settings) == false!\n");
			}
		}

		/// \brief Generate terrain mesh vertices using Perlin noise and Poisson disk sampling.
		void GeneratePoints()
		{
			PoissonSamplePointsInXYPlane();
			GeneratePerlinZElevation();
		}

		/// ==========================================================================================
		///	\brief Triangulate terrain points (generates PolyIndices as triples).
		/// \param[in] triangulate        a custom triangulation function.
		///	==========================================================================================
		void Triangulate(const TriangulationFunction& triangulate)
		{
			if (m_Result.Vertices.empty())
			{
				std::cerr << "Geometry::TerrainBuilder::Triangulate: nothing to triangulate! m_Result.Vertices.empty()!\n";
				return;
			}
			triangulate(m_Result);
		}

	private:

		/// \brief Generates the proper amount of terrain points, and computes their x,y coordinates within the boundary loop.
		void PoissonSamplePointsInXYPlane();

		/// \brief Computes the z-coordinate of each generated point using Perlin noise within range [m_Settings.MinElevation, m_Settings.MaxElevation].
		void GeneratePerlinZElevation();

		TerrainSettings m_Settings; //>! settings for this builder.

		BaseMeshGeometryData m_Result{}; //>! result
	};
	
} // namespace Geometry