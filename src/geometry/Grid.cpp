#include "Grid.h"

namespace Geometry
{
	ScalarGrid::ScalarGrid(const float& cellSize, const pmp::BoundingBox& box)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		const int nZMinus = std::floor(box.min()[2] / m_CellSize);
		const int nZPlus = std::ceil(box.max()[2] / m_CellSize);

		// adjust box
		const auto minVec = pmp::vec3(nXMinus, nYMinus, nZMinus) * m_CellSize;
		const auto maxVec = pmp::vec3(nXPlus, nYPlus, nZPlus) * m_CellSize;
		m_Box = pmp::BoundingBox(minVec, maxVec);

		m_Dimensions = GridDimensions{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus),
			static_cast<size_t>(nZPlus - nZMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny * m_Dimensions.Nz;
		m_Values = std::vector(nValues, DEFAULT_SCALAR_GRID_INIT_VAL);
		m_FrozenValues = std::vector(nValues, false);
	}

	ScalarGrid::ScalarGrid(const float& cellSize, const pmp::BoundingBox& box, const double& initVal)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		const int nZMinus = std::floor(box.min()[2] / m_CellSize);
		const int nZPlus = std::ceil(box.max()[2] / m_CellSize);

		// adjust box
		const auto maxVec = pmp::vec3(nXPlus, nYPlus, nZPlus) * m_CellSize;
		const auto minVec = pmp::vec3(nXMinus, nYMinus, nZMinus) * m_CellSize;
		m_Box = pmp::BoundingBox(minVec, maxVec);

		m_Dimensions = GridDimensions{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus),
			static_cast<size_t>(nZPlus - nZMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny * m_Dimensions.Nz;
		m_Values = std::vector(nValues, initVal);
		m_FrozenValues = std::vector(nValues, false);
	}

	bool ScalarGrid::operator==(const ScalarGrid& other) const
	{
		if (m_Values.size() != other.m_Values.size())
			return false;

		return memcmp(&m_Values.front(), &other.m_Values.front(), sizeof(m_Values[0]) * m_Values.size()) == 0;
	}

	bool ScalarGrid::operator!=(const ScalarGrid& other) const
	{
		if (m_Values.size() != other.m_Values.size())
			return true;

		return memcmp(&m_Values.front(), &other.m_Values.front(), sizeof(m_Values[0]) * m_Values.size()) != 0;
	}

	ScalarGrid& ScalarGrid::operator*=(const double& scalar)
	{
		for (auto& val : m_Values)
			val = scalar * val;
		return *this;
	}

	ScalarGrid& ScalarGrid::operator*=(const ScalarGrid& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
			m_Values[i] *= other.m_Values[i];
		return *this;
	}

	ScalarGrid& ScalarGrid::operator/=(const double& scalar)
	{
		if (scalar == 0.0)
		{
			assert(false);
			return *this;
		}

		for (auto& val : m_Values)
			val /= scalar;
		return *this;
	}

	ScalarGrid& ScalarGrid::operator/=(const ScalarGrid& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
		{
			assert(other.m_Values[i] != 0.0);
			m_Values[i] /= other.m_Values[i];
		}
		return *this;
	}

	ScalarGrid& ScalarGrid::operator+=(const double& scalar)
	{
		for (auto& val : m_Values)
			val += scalar;
		return *this;
	}

	ScalarGrid& ScalarGrid::operator+=(const ScalarGrid& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
			m_Values[i] += other.m_Values[i];
		return *this;
	}

	ScalarGrid& ScalarGrid::operator-=(const double& scalar)
	{
		for (auto& val : m_Values)
			val -= scalar;
		return *this;
	}

	ScalarGrid& ScalarGrid::operator-=(const ScalarGrid& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
			m_Values[i] -= other.m_Values[i];
		return *this;
	}

	ScalarGrid ScalarGrid::operator*(const double& scalar) const
	{
		ScalarGrid result(*this);
		result *= scalar;
		return result;
	}

	ScalarGrid ScalarGrid::operator*(const ScalarGrid& other) const
	{
		ScalarGrid result(*this);
		result *= other;
		return result;
	}

	ScalarGrid ScalarGrid::operator/(const double& scalar) const
	{
		ScalarGrid result(*this);
		result /= scalar;
		return result;
	}

	ScalarGrid ScalarGrid::operator/(const ScalarGrid& other) const
	{
		ScalarGrid result(*this);
		result /= other;
		return result;
	}

	ScalarGrid ScalarGrid::operator+(const double& scalar) const
	{
		ScalarGrid result(*this);
		result += scalar;
		return result;
	}

	ScalarGrid ScalarGrid::operator+(const ScalarGrid& other) const
	{
		ScalarGrid result(*this);
		result += other;
		return result;
	}

	ScalarGrid ScalarGrid::operator-(const double& scalar) const
	{
		ScalarGrid result(*this);
		result -= scalar;
		return result;
	}

	ScalarGrid ScalarGrid::operator-(const ScalarGrid& other) const
	{
		ScalarGrid result(*this);
		result -= other;
		return result;
	}

	ScalarGrid& ScalarGrid::operator*=(const pmp::mat4& mat)
	{
		// assuming mat is a uniform scaling + translation matrix only.
		m_CellSize *= mat(0, 0);
		m_Box *= mat;
		return *this;
	}

	bool ScalarGrid::IsValid() const
	{
		if (m_CellSize <= 0.0f)
			return false;

		if (m_Box.is_empty())
			return false;

		if (m_Values.empty())
			return false;

		return m_Dimensions.Valid();
	}

	VectorGrid::VectorGrid(const ScalarGrid& scalarGrid)
		: m_Box(scalarGrid.Box()), m_Dimensions(scalarGrid.Dimensions()), m_CellSize(scalarGrid.CellSize())
	{
		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny * m_Dimensions.Nz;
		m_ValuesX = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_ValuesY = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_ValuesZ = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_FrozenValues = std::vector(nValues, false);
	}

	VectorGrid::VectorGrid(const float& cellSize, const pmp::BoundingBox& box)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		const int nZMinus = std::floor(box.min()[2] / m_CellSize);
		const int nZPlus = std::ceil(box.max()[2] / m_CellSize);

		// adjust box
		const auto minVec = pmp::vec3(nXMinus, nYMinus, nZMinus) * m_CellSize;
		const auto maxVec = pmp::vec3(nXPlus, nYPlus, nZPlus) * m_CellSize;
		m_Box = pmp::BoundingBox(minVec, maxVec);

		m_Dimensions = GridDimensions{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus),
			static_cast<size_t>(nZPlus - nZMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny * m_Dimensions.Nz;
		m_ValuesX = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_ValuesY = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_ValuesZ = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_FrozenValues = std::vector(nValues, false);
	}

	VectorGrid::VectorGrid(const float& cellSize, const pmp::BoundingBox& box, const pmp::vec3& initVal)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		const int nZMinus = std::floor(box.min()[2] / m_CellSize);
		const int nZPlus = std::ceil(box.max()[2] / m_CellSize);

		// adjust box
		const auto maxVec = pmp::vec3(nXPlus, nYPlus, nZPlus) * m_CellSize;
		const auto minVec = pmp::vec3(nXMinus, nYMinus, nZMinus) * m_CellSize;
		m_Box = pmp::BoundingBox(minVec, maxVec);

		m_Dimensions = GridDimensions{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus),
			static_cast<size_t>(nZPlus - nZMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny * m_Dimensions.Nz;
		m_ValuesX = std::vector(nValues, static_cast<double>(initVal[0]));
		m_ValuesY = std::vector(nValues, static_cast<double>(initVal[1]));
		m_ValuesZ = std::vector(nValues, static_cast<double>(initVal[2]));
		m_FrozenValues = std::vector(nValues, false);
	}

	bool VectorGrid::operator==(const VectorGrid& other) const
	{
		if (m_ValuesX.size() != other.m_ValuesX.size())
			return false;

		return 
			(memcmp(&m_ValuesX.front(), &other.m_ValuesX.front(), sizeof(m_ValuesX[0]) * m_ValuesX.size()) == 0) &&
			(memcmp(&m_ValuesX.front(), &other.m_ValuesY.front(), sizeof(m_ValuesY[0]) * m_ValuesY.size()) == 0) &&
			(memcmp(&m_ValuesX.front(), &other.m_ValuesZ.front(), sizeof(m_ValuesZ[0]) * m_ValuesZ.size()) == 0);
	}

	bool VectorGrid::operator!=(const VectorGrid& other) const
	{
		if (m_ValuesX.size() != other.m_ValuesX.size())
			return true;

		return
			(memcmp(&m_ValuesX.front(), &other.m_ValuesX.front(), sizeof(m_ValuesX[0]) * m_ValuesX.size()) != 0) ||
			(memcmp(&m_ValuesX.front(), &other.m_ValuesY.front(), sizeof(m_ValuesY[0]) * m_ValuesY.size()) != 0) ||
			(memcmp(&m_ValuesX.front(), &other.m_ValuesZ.front(), sizeof(m_ValuesZ[0]) * m_ValuesZ.size()) != 0);
	}

	VectorGrid& VectorGrid::operator*=(const double& scalar)
	{
		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] *= scalar;
			m_ValuesY[i] *= scalar;
			m_ValuesZ[i] *= scalar;
		}
		return *this;
	}

	VectorGrid& VectorGrid::operator*=(const VectorGrid& other)
	{
		if (m_ValuesX.size() != other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] *= other.m_ValuesX[i];
			m_ValuesY[i] *= other.m_ValuesY[i];
			m_ValuesZ[i] *= other.m_ValuesZ[i];
		}
		return *this;
	}

	VectorGrid& VectorGrid::operator/=(const double& scalar)
	{
		if (scalar == 0.0)
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] /= scalar;
			m_ValuesY[i] /= scalar;
			m_ValuesZ[i] /= scalar;
		}
		return *this;
	}

	VectorGrid& VectorGrid::operator/=(const VectorGrid& other)
	{
		if (m_ValuesX.size() == other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			assert(other.m_ValuesX[i] != 0.0);
			m_ValuesX[i] /= other.m_ValuesX[i];
			assert(other.m_ValuesY[i] != 0.0);
			m_ValuesY[i] /= other.m_ValuesY[i];
			assert(other.m_ValuesZ[i] != 0.0);
			m_ValuesZ[i] /= other.m_ValuesZ[i];
		}
		return *this;
	}

	VectorGrid& VectorGrid::operator+=(const double& scalar)
	{
		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] *= scalar;
			m_ValuesY[i] *= scalar;
			m_ValuesZ[i] *= scalar;
		}
		return *this;
	}

	VectorGrid& VectorGrid::operator+=(const VectorGrid& other)
	{
		if (m_ValuesX.size() == other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] += other.m_ValuesX[i];
			m_ValuesY[i] += other.m_ValuesY[i];
			m_ValuesZ[i] += other.m_ValuesZ[i];
		}
		return *this;
	}

	VectorGrid& VectorGrid::operator-=(const double& scalar)
	{
		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] -= scalar;
			m_ValuesY[i] -= scalar;
			m_ValuesZ[i] -= scalar;
		}
		return *this;
	}

	VectorGrid& VectorGrid::operator-=(const VectorGrid& other)
	{
		if (m_ValuesX.size() == other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] -= other.m_ValuesX[i];
			m_ValuesY[i] -= other.m_ValuesY[i];
			m_ValuesZ[i] -= other.m_ValuesZ[i];
		}
		return *this;
	}

	VectorGrid VectorGrid::operator*(const double& scalar) const
	{
		VectorGrid result(*this);
		result *= scalar;
		return result;
	}

	VectorGrid VectorGrid::operator*(const VectorGrid& other) const
	{
		VectorGrid result(*this);
		result *= other;
		return result;
	}

	VectorGrid VectorGrid::operator/(const double& scalar) const
	{
		VectorGrid result(*this);
		result /= scalar;
		return result;
	}

	VectorGrid VectorGrid::operator/(const VectorGrid& other) const
	{
		VectorGrid result(*this);
		result /= other;
		return result;
	}

	VectorGrid VectorGrid::operator+(const double& scalar) const
	{
		VectorGrid result(*this);
		result += scalar;
		return result;
	}

	VectorGrid VectorGrid::operator+(const VectorGrid& other) const
	{
		VectorGrid result(*this);
		result += other;
		return result;
	}

	VectorGrid VectorGrid::operator-(const double& scalar) const
	{
		VectorGrid result(*this);
		result -= scalar;
		return result;
	}

	VectorGrid VectorGrid::operator-(const VectorGrid& other) const
	{
		VectorGrid result(*this);
		result -= other;
		return result;
	}

	VectorGrid& VectorGrid::operator*=(const pmp::mat4& mat)
	{
		// assuming mat is a uniform scaling + translation matrix only.
		m_CellSize *= mat(0, 0);
		m_Box *= mat;
		return *this;
	}

	bool VectorGrid::IsValid() const
	{
		if (m_CellSize <= 0.0f)
			return false;

		if (m_Box.is_empty())
			return false;

		if (m_ValuesX.empty())
			return false;
		if (m_ValuesY.empty())
			return false;
		if (m_ValuesZ.empty())
			return false;

		if (m_ValuesX.size() != m_ValuesY.size())
			return false;
		if (m_ValuesY.size() != m_ValuesZ.size())
			return false;

		return m_Dimensions.Valid();
	}

	ScalarGrid2D::ScalarGrid2D(const float& cellSize, const pmp::BoundingBox2& box)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		// adjust box
		const auto minVec = pmp::vec2(nXMinus, nYMinus) * m_CellSize;
		const auto maxVec = pmp::vec2(nXPlus, nYPlus) * m_CellSize;
		m_Box = pmp::BoundingBox2(minVec, maxVec);

		m_Dimensions = GridDimensions2D{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny;
		m_Values = std::vector(nValues, DEFAULT_SCALAR_GRID_INIT_VAL);
		m_FrozenValues = std::vector(nValues, false);
	}

	ScalarGrid2D::ScalarGrid2D(const float& cellSize, const pmp::BoundingBox2& box, const double& initVal)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		// adjust box
		const auto maxVec = pmp::vec2(nXPlus, nYPlus) * m_CellSize;
		const auto minVec = pmp::vec2(nXMinus, nYMinus) * m_CellSize;
		m_Box = pmp::BoundingBox2(minVec, maxVec);

		m_Dimensions = GridDimensions2D{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny;
		m_Values = std::vector(nValues, initVal);
		m_FrozenValues = std::vector(nValues, false);
	}

	bool ScalarGrid2D::operator==(const ScalarGrid2D& other) const
	{
		if (m_Values.size() != other.m_Values.size())
			return false;

		return memcmp(&m_Values.front(), &other.m_Values.front(), sizeof(m_Values[0]) * m_Values.size()) == 0;
	}

	bool ScalarGrid2D::operator!=(const ScalarGrid2D& other) const
	{
		if (m_Values.size() != other.m_Values.size())
			return true;

		return memcmp(&m_Values.front(), &other.m_Values.front(), sizeof(m_Values[0]) * m_Values.size()) != 0;
	}

	ScalarGrid2D& ScalarGrid2D::operator*=(const double& scalar)
	{
		for (auto& val : m_Values)
			val = scalar * val;
		return *this;
	}

	ScalarGrid2D& ScalarGrid2D::operator*=(const ScalarGrid2D& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
			m_Values[i] *= other.m_Values[i];
		return *this;
	}

	ScalarGrid2D& ScalarGrid2D::operator/=(const double& scalar)
	{
		if (scalar == 0.0)
		{
			assert(false);
			return *this;
		}

		for (auto& val : m_Values)
			val /= scalar;
		return *this;
	}

	ScalarGrid2D& ScalarGrid2D::operator/=(const ScalarGrid2D& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
		{
			assert(other.m_Values[i] != 0.0);
			m_Values[i] /= other.m_Values[i];
		}
		return *this;
	}

	ScalarGrid2D& ScalarGrid2D::operator+=(const double& scalar)
	{
		for (auto& val : m_Values)
			val += scalar;
		return *this;
	}

	ScalarGrid2D& ScalarGrid2D::operator+=(const ScalarGrid2D& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
			m_Values[i] += other.m_Values[i];
		return *this;
	}

	ScalarGrid2D& ScalarGrid2D::operator-=(const double& scalar)
	{
		for (auto& val : m_Values)
			val -= scalar;
		return *this;
	}

	ScalarGrid2D& ScalarGrid2D::operator-=(const ScalarGrid2D& other)
	{
		if (m_Values.size() == other.m_Values.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_Values.size(); i++)
			m_Values[i] -= other.m_Values[i];
		return *this;
	}

	ScalarGrid2D ScalarGrid2D::operator*(const double& scalar) const
	{
		ScalarGrid2D result(*this);
		result *= scalar;
		return result;
	}

	ScalarGrid2D ScalarGrid2D::operator*(const ScalarGrid2D& other) const
	{
		ScalarGrid2D result(*this);
		result *= other;
		return result;
	}

	ScalarGrid2D ScalarGrid2D::operator/(const double& scalar) const
	{
		ScalarGrid2D result(*this);
		result /= scalar;
		return result;
	}

	ScalarGrid2D ScalarGrid2D::operator/(const ScalarGrid2D& other) const
	{
		ScalarGrid2D result(*this);
		result /= other;
		return result;
	}

	ScalarGrid2D ScalarGrid2D::operator+(const double& scalar) const
	{
		ScalarGrid2D result(*this);
		result += scalar;
		return result;
	}

	ScalarGrid2D ScalarGrid2D::operator+(const ScalarGrid2D& other) const
	{
		ScalarGrid2D result(*this);
		result += other;
		return result;
	}

	ScalarGrid2D ScalarGrid2D::operator-(const double& scalar) const
	{
		ScalarGrid2D result(*this);
		result -= scalar;
		return result;
	}

	ScalarGrid2D ScalarGrid2D::operator-(const ScalarGrid2D& other) const
	{
		ScalarGrid2D result(*this);
		result -= other;
		return result;
	}

	ScalarGrid2D& ScalarGrid2D::operator*=(const pmp::mat3& mat)
	{
		// assuming mat is a uniform scaling + translation matrix only.
		m_CellSize *= mat(0, 0);
		m_Box *= mat;
		return *this;
	}

	bool ScalarGrid2D::IsValid() const
	{
		if (m_CellSize <= 0.0f)
			return false;

		if (m_Box.is_empty())
			return false;

		if (m_Values.empty())
			return false;

		return m_Dimensions.Valid();
	}

	VectorGrid2D::VectorGrid2D(const ScalarGrid2D& scalarGrid)
		: m_Box(scalarGrid.Box()), m_Dimensions(scalarGrid.Dimensions()), m_CellSize(scalarGrid.CellSize())
	{
		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny;
		m_ValuesX = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_ValuesY = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_FrozenValues = std::vector(nValues, false);
	}

	VectorGrid2D::VectorGrid2D(const float& cellSize, const pmp::BoundingBox2& box)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		// adjust box
		const auto minVec = pmp::vec2(nXMinus, nYMinus) * m_CellSize;
		const auto maxVec = pmp::vec2(nXPlus, nYPlus) * m_CellSize;
		m_Box = pmp::BoundingBox2(minVec, maxVec);

		m_Dimensions = GridDimensions2D{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny;
		m_ValuesX = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_ValuesY = std::vector(nValues, DEFAULT_VECTOR_GRID_INIT_VAL);
		m_FrozenValues = std::vector(nValues, false);
	}

	VectorGrid2D::VectorGrid2D(const float& cellSize, const pmp::BoundingBox2& box, const pmp::vec2& initVal)
		: m_CellSize(cellSize)
	{
		// compute global grid bounds
		const int nXMinus = std::floor(box.min()[0] / m_CellSize);
		const int nXPlus = std::ceil(box.max()[0] / m_CellSize);

		const int nYMinus = std::floor(box.min()[1] / m_CellSize);
		const int nYPlus = std::ceil(box.max()[1] / m_CellSize);

		// adjust box
		const auto minVec = pmp::vec2(nXMinus, nYMinus) * m_CellSize;
		const auto maxVec = pmp::vec2(nXPlus, nYPlus) * m_CellSize;
		m_Box = pmp::BoundingBox2(minVec, maxVec);

		m_Dimensions = GridDimensions2D{
			static_cast<size_t>(nXPlus - nXMinus),
			static_cast<size_t>(nYPlus - nYMinus)
		};

		const size_t nValues = m_Dimensions.Nx * m_Dimensions.Ny;
		m_ValuesX = std::vector(nValues, static_cast<double>(initVal[0]));
		m_ValuesY = std::vector(nValues, static_cast<double>(initVal[1]));
		m_FrozenValues = std::vector(nValues, false);
	}

	bool VectorGrid2D::operator==(const VectorGrid2D& other) const
	{
		if (m_ValuesX.size() != other.m_ValuesX.size())
			return false;

		return
			(memcmp(&m_ValuesX.front(), &other.m_ValuesX.front(), sizeof(m_ValuesX[0]) * m_ValuesX.size()) == 0) &&
			(memcmp(&m_ValuesX.front(), &other.m_ValuesY.front(), sizeof(m_ValuesY[0]) * m_ValuesY.size()) == 0);
	}

	bool VectorGrid2D::operator!=(const VectorGrid2D& other) const
	{
		if (m_ValuesX.size() != other.m_ValuesX.size())
			return true;

		return
			(memcmp(&m_ValuesX.front(), &other.m_ValuesX.front(), sizeof(m_ValuesX[0]) * m_ValuesX.size()) != 0) ||
			(memcmp(&m_ValuesX.front(), &other.m_ValuesY.front(), sizeof(m_ValuesY[0]) * m_ValuesY.size()) != 0);
	}

	VectorGrid2D& VectorGrid2D::operator*=(const double& scalar)
	{
		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] *= scalar;
			m_ValuesY[i] *= scalar;
		}
		return *this;
	}

	VectorGrid2D& VectorGrid2D::operator*=(const VectorGrid2D& other)
	{
		if (m_ValuesX.size() != other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] *= other.m_ValuesX[i];
			m_ValuesY[i] *= other.m_ValuesY[i];
		}
		return *this;
	}


	VectorGrid2D& VectorGrid2D::operator/=(const double& scalar)
	{
		if (scalar == 0.0)
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] /= scalar;
			m_ValuesY[i] /= scalar;
		}
		return *this;
	}

	VectorGrid2D& VectorGrid2D::operator/=(const VectorGrid2D& other)
	{
		if (m_ValuesX.size() == other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			assert(other.m_ValuesX[i] != 0.0);
			m_ValuesX[i] /= other.m_ValuesX[i];
			assert(other.m_ValuesY[i] != 0.0);
			m_ValuesY[i] /= other.m_ValuesY[i];
		}
		return *this;
	}

	VectorGrid2D& VectorGrid2D::operator+=(const double& scalar)
	{
		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] *= scalar;
			m_ValuesY[i] *= scalar;
		}
		return *this;
	}

	VectorGrid2D& VectorGrid2D::operator+=(const VectorGrid2D& other)
	{
		if (m_ValuesX.size() == other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] += other.m_ValuesX[i];
			m_ValuesY[i] += other.m_ValuesY[i];
		}
		return *this;
	}

	VectorGrid2D& VectorGrid2D::operator-=(const double& scalar)
	{
		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] -= scalar;
			m_ValuesY[i] -= scalar;
		}
		return *this;
	}

	VectorGrid2D& VectorGrid2D::operator-=(const VectorGrid2D& other)
	{
		if (m_ValuesX.size() == other.m_ValuesX.size())
		{
			assert(false);
			return *this;
		}

		for (unsigned int i = 0; i < m_ValuesX.size(); i++)
		{
			m_ValuesX[i] -= other.m_ValuesX[i];
			m_ValuesY[i] -= other.m_ValuesY[i];
		}
		return *this;
	}

	VectorGrid2D VectorGrid2D::operator*(const double& scalar) const
	{
		VectorGrid2D result(*this);
		result *= scalar;
		return result;
	}

	VectorGrid2D VectorGrid2D::operator*(const VectorGrid2D& other) const
	{
		VectorGrid2D result(*this);
		result *= other;
		return result;
	}

	VectorGrid2D VectorGrid2D::operator/(const double& scalar) const
	{
		VectorGrid2D result(*this);
		result /= scalar;
		return result;
	}

	VectorGrid2D VectorGrid2D::operator/(const VectorGrid2D& other) const
	{
		VectorGrid2D result(*this);
		result /= other;
		return result;
	}

	VectorGrid2D VectorGrid2D::operator+(const double& scalar) const
	{
		VectorGrid2D result(*this);
		result += scalar;
		return result;
	}

	VectorGrid2D VectorGrid2D::operator+(const VectorGrid2D& other) const
	{
		VectorGrid2D result(*this);
		result += other;
		return result;
	}

	VectorGrid2D VectorGrid2D::operator-(const double& scalar) const
	{
		VectorGrid2D result(*this);
		result -= scalar;
		return result;
	}

	VectorGrid2D VectorGrid2D::operator-(const VectorGrid2D& other) const
	{
		VectorGrid2D result(*this);
		result -= other;
		return result;
	}

	VectorGrid2D& VectorGrid2D::operator*=(const pmp::mat3& mat)
	{
		// assuming mat is a uniform scaling + translation matrix only.
		m_CellSize *= mat(0, 0);
		m_Box *= mat;
		return *this;
	}

	bool VectorGrid2D::IsValid() const
	{
		if (m_CellSize <= 0.0f)
			return false;

		if (m_Box.is_empty())
			return false;

		if (m_ValuesX.empty())
			return false;
		if (m_ValuesY.empty())
			return false;

		if (m_ValuesX.size() != m_ValuesY.size())
			return false;

		return m_Dimensions.Valid();
	}

} // namespace Geometry
