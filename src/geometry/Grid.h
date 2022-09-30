#pragma once

#include "pmp/BoundingBox.h"

#include <vector>


namespace Geometry
{
	/// \brief default value to initialize a scalar grid with.
	constexpr double DEFAULT_SCALAR_GRID_INIT_VAL = 1e+9;
	/// \brief default value to initialize a vector grid with.
	constexpr double DEFAULT_VECTOR_GRID_INIT_VAL = 0.0;

	/// \brief dimensions wrapper item for any grid object.
	struct GridDimensions
	{
		size_t Nx{ 0 };
		size_t Ny{ 0 };
		size_t Nz{ 0 };

		[[nodiscard]] bool Valid() const
		{
			return Nx > 0 && Ny > 0 && Nz > 0;
		}
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

		// ====== Operators ======================

		bool operator==(const ScalarGrid& other) const;
		bool operator!=(const ScalarGrid& other) const;
		ScalarGrid& operator*= (const double& scalar);
		ScalarGrid& operator*= (const ScalarGrid& other);
		ScalarGrid& operator/= (const double& scalar);
		ScalarGrid& operator/= (const ScalarGrid& other);
		ScalarGrid& operator+= (const double& scalar);
		ScalarGrid& operator+= (const ScalarGrid& other);
		ScalarGrid& operator-= (const double& scalar);
		ScalarGrid& operator-= (const ScalarGrid& other);
		ScalarGrid operator* (const double& scalar) const;
		ScalarGrid operator* (const ScalarGrid& other) const;
		ScalarGrid operator/ (const double& scalar) const;
		ScalarGrid operator/ (const ScalarGrid& other) const;
		ScalarGrid operator+ (const double& scalar) const;
		ScalarGrid operator+ (const ScalarGrid& other) const;
		ScalarGrid operator- (const double& scalar) const;
		ScalarGrid operator- (const ScalarGrid& other) const;

		ScalarGrid& operator*= (const pmp::mat4& mat);

		// ====== Validity ======================

		[[nodiscard]] bool IsValid() const;

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
		/**
		 * \brief Constructor. Initializes from a scalar grid with default initialization vector value.
		 * \param scalarGrid      scalar grid to initialize from.
		 */
		explicit VectorGrid(const ScalarGrid& scalarGrid);

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

		// ====== Operators ======================

		bool operator==(const VectorGrid& other) const;
		bool operator!=(const VectorGrid& other) const;
		VectorGrid& operator*= (const double& scalar);
		VectorGrid& operator*= (const VectorGrid& other);
		VectorGrid& operator/= (const double& scalar);
		VectorGrid& operator/= (const VectorGrid& other);
		VectorGrid& operator+= (const double& scalar);
		VectorGrid& operator+= (const VectorGrid& other);
		VectorGrid& operator-= (const double& scalar);
		VectorGrid& operator-= (const VectorGrid& other);
		VectorGrid operator* (const double& scalar) const;
		VectorGrid operator* (const VectorGrid& other) const;
		VectorGrid operator/ (const double& scalar) const;
		VectorGrid operator/ (const VectorGrid& other) const;
		VectorGrid operator+ (const double& scalar) const;
		VectorGrid operator+ (const VectorGrid& other) const;
		VectorGrid operator- (const double& scalar) const;
		VectorGrid operator- (const VectorGrid& other) const;

		// ====== Validity ======================

		[[nodiscard]] bool IsValid() const;

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