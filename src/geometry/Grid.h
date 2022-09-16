#pragma once

#include "pmp/BoundingBox.h"

#include <vector>


namespace Geometry
{
	/// \brief dimensions wrapper item for any grid object.
	struct GridDimensions
	{
		size_t Nx;
		size_t Ny;
		size_t Nz;
	};

	/// \brief A 3D grid object containing scalar values and additional flags for individual voxels.
	class ScalarGrid
	{
	public:
		ScalarGrid(const float& cellSize, const pmp::BoundingBox& box);

		ScalarGrid(const float& cellSize, const pmp::BoundingBox& box, const double& initVal);

		// ====== Getters ======================

		pmp::BoundingBox& Box()
		{
			return m_Box;
		}

		[[nodiscard]] const pmp::BoundingBox& Box() const
		{
			return m_Box;
		}

		GridDimensions& Dimensions()
		{
			return m_Dimensions;
		}

		[[nodiscard]] const GridDimensions& Dimensions() const
		{
			return m_Dimensions;
		}

		float& CellSize()
		{
			return m_CellSize;
		}

		[[nodiscard]] const float& CellSize() const
		{
			return m_CellSize;
		}

		std::vector<double>& Values()
		{
			return m_Values;
		}

		[[nodiscard]] const std::vector<double>& Values() const
		{
			return m_Values;
		}

		std::vector<bool>& FrozenValues()
		{
			return m_FrozenValues;
		}

		[[nodiscard]] const std::vector<bool>& FrozenValues() const
		{
			return m_FrozenValues;
		}

	private:
		pmp::BoundingBox m_Box{};
		GridDimensions m_Dimensions{};
		float m_CellSize{};
		std::vector<double> m_Values{};
		std::vector<bool> m_FrozenValues{};
	};

	/// \brief A 3D grid object containing vector values and additional flags for individual voxels.
	class VectorGrid
	{
	public:
		VectorGrid(const float& cellSize, const pmp::BoundingBox& box);

		VectorGrid(const float& cellSize, const pmp::BoundingBox& box, const pmp::vec3& initVal);

		// ====== Getters ======================

		pmp::BoundingBox& Box()
		{
			return m_Box;
		}

		[[nodiscard]] const pmp::BoundingBox& Box() const
		{
			return m_Box;
		}

		GridDimensions& Dimensions()
		{
			return m_Dimensions;
		}

		[[nodiscard]] const GridDimensions& Dimensions() const
		{
			return m_Dimensions;
		}

		float& CellSize()
		{
			return m_CellSize;
		}

		[[nodiscard]] const float& CellSize() const
		{
			return m_CellSize;
		}

		std::vector<double>& ValuesX()
		{
			return m_ValuesX;
		}

		[[nodiscard]] const std::vector<double>& ValuesX() const
		{
			return m_ValuesX;
		}

		std::vector<double>& ValuesY()
		{
			return m_ValuesY;
		}

		[[nodiscard]] const std::vector<double>& ValuesY() const
		{
			return m_ValuesY;
		}

		std::vector<double>& ValuesZ()
		{
			return m_ValuesZ;
		}

		[[nodiscard]] const std::vector<double>& ValuesZ() const
		{
			return m_ValuesZ;
		}

		std::vector<bool>& FrozenValues()
		{
			return m_FrozenValues;
		}

		[[nodiscard]] const std::vector<bool>& FrozenValues() const
		{
			return m_FrozenValues;
		}

	private:
		pmp::BoundingBox m_Box{};
		GridDimensions m_Dimensions{};
		float m_CellSize{};

		std::vector<double> m_ValuesX{};
		std::vector<double> m_ValuesY{};
		std::vector<double> m_ValuesZ{};

		std::vector<bool> m_FrozenValues{};
	};
	
} // namespace Geometry