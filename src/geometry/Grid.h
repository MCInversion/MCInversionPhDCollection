#pragma once

#include "pmp/BoundingBox.h"

#include <vector>


namespace Geometry
{
	struct GridDimensions
	{
		size_t Nx;
		size_t Ny;
		size_t Nz;
	};


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
	
} // namespace Geometry