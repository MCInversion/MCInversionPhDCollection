#pragma once

#include "pmp/SurfaceMesh.h"
#include "pmp/ManifoldCurve2D.h"
#include "GeometryConversionUtils.h"

namespace Geometry
{
	/**
	 * \brief Abstract class that defines a common interface for mesh adapters.
	 * \class MeshAdapter
	 */
	class MeshAdapter
	{
	public:
		/**
		 * \brief Creates a unique pointer to a new MeshAdapter that is a copy of the current object.
		 * \return A unique pointer to a new MeshAdapter.
		 */
		virtual [[nodiscard]] std::unique_ptr<MeshAdapter> Clone() const = 0;

		/// \brief Verifies whether this adapter applies to an empty mesh.
		virtual [[nodiscard]] bool IsEmpty() const = 0;

		/**
		 * \brief Retrieves the vertices of the mesh.
		 * \return A vector containing the vertices of the mesh.
		 */
		virtual [[nodiscard]] std::vector<pmp::Point> GetVertices() = 0;

		/**
		 * \brief Retrieves the polygon indices of the mesh.
		 * \return A vector of vectors, where each inner vector contains the indices of a polygon.
		 */
		virtual [[nodiscard]] std::vector<std::vector<unsigned int>> GetPolyIndices() = 0;

		/**
		 * \brief Checks if the mesh consists solely of triangles.
		 * \return True if the mesh is composed only of triangles, false otherwise.
		 */
		virtual [[nodiscard]] bool IsTriangle() const = 0;

		/**
		 * \brief Retrieves the bounding box of the mesh.
		 * \return A BoundingBox object representing the bounding box of the mesh.
		 */
		virtual [[nodiscard]] pmp::BoundingBox GetBounds() const = 0;

		/**
		 * \brief Virtual destructor for MeshAdapter.
		 */
		virtual ~MeshAdapter() = default;
	};


	/**
	 * \brief Concrete class extending MeshAdapter to provide specific functionality for BaseMeshGeometryData.
	 * \class BaseMeshAdapter
	 */
	class BaseMeshAdapter : public MeshAdapter
	{
	public:
		/**
		 * \brief Constructor initializing the BaseMeshAdapter with given mesh data.
		 * \param mesh Shared pointer to BaseMeshGeometryData.
		 */
		explicit BaseMeshAdapter(std::shared_ptr<BaseMeshGeometryData> mesh) : m_BaseMesh(std::move(mesh)) {}

		/// \brief Verifies whether this adapter applies to an empty mesh.
		[[nodiscard]] bool IsEmpty() const override
		{
			if (!m_BaseMesh)
				return true;

			return m_BaseMesh->Vertices.empty() || m_BaseMesh->PolyIndices.empty();
		}

		/// \copydoc MeshAdapter::Clone
		[[nodiscard]] std::unique_ptr<MeshAdapter> Clone() const override
		{
			return std::make_unique<BaseMeshAdapter>(*this);
		}

		/**
		 * \brief Retrieves the base mesh geometry data.
		 * \return Reference to BaseMeshGeometryData.
		 */
		BaseMeshGeometryData& GetBaseMesh()
		{
			return *m_BaseMesh;
		}

		/**
		 * \brief Retrieves the base mesh geometry data (const version).
		 * \return Constant reference to BaseMeshGeometryData.
		 */
		[[nodiscard]] const BaseMeshGeometryData& GetBaseMesh() const
		{
			return *m_BaseMesh;
		}

		/// \copydoc MeshAdapter::GetVertices
		[[nodiscard]] std::vector<pmp::Point> GetVertices() override
		{
			return m_BaseMesh->Vertices;
		}

		/// \copydoc MeshAdapter::GetPolyIndices
		[[nodiscard]] std::vector<std::vector<unsigned int>> GetPolyIndices() override
		{
			return m_BaseMesh->PolyIndices;
		}

		/// \copydoc MeshAdapter::IsTriangle
		[[nodiscard]] bool IsTriangle() const override
		{
			return std::ranges::all_of(
				m_BaseMesh->PolyIndices, [](const auto& poly) { return poly.size() == 3; });
		}

		/// \copydoc MeshAdapter::GetBounds
		[[nodiscard]] pmp::BoundingBox GetBounds() const override
		{
			return pmp::BoundingBox(m_BaseMesh->Vertices);
		}

	private:
		std::shared_ptr<BaseMeshGeometryData> m_BaseMesh; //>! actual geometry instance
	};

	/**
	 * \brief Concrete class extending MeshAdapter to provide specific functionality for pmp::SurfaceMesh.
	 * \class PMPSurfaceMeshAdapter
	 */
	class PMPSurfaceMeshAdapter : public MeshAdapter
	{
	public:
		/**
		 * \brief Constructor initializing the PMPSurfaceMeshAdapter with given mesh data.
		 * \param mesh Shared pointer to pmp::SurfaceMesh.
		 */
		explicit PMPSurfaceMeshAdapter(std::shared_ptr<pmp::SurfaceMesh> mesh) : m_SurfaceMesh(std::move(mesh)) {}

		/// \brief Verifies whether this adapter applies to an empty mesh.
		[[nodiscard]] bool IsEmpty() const override
		{
			if (!m_SurfaceMesh)
				return true;

			return m_SurfaceMesh->is_empty();
		}

		/// \copydoc MeshAdapter::Clone
		[[nodiscard]] std::unique_ptr<MeshAdapter> Clone() const override
		{
			return std::make_unique<PMPSurfaceMeshAdapter>(*this);
		}

		/// \copydoc MeshAdapter::GetVertices
		[[nodiscard]] std::vector<pmp::Point> GetVertices() override
		{
			std::vector<pmp::Point> vertices;
			vertices.reserve(m_SurfaceMesh->n_vertices());
			for (const auto v : m_SurfaceMesh->vertices())
			{
				vertices.push_back(m_SurfaceMesh->position(v));
			}
			return vertices;
		}

		/// \copydoc MeshAdapter::GetPolyIndices
		[[nodiscard]] std::vector<std::vector<unsigned int>> GetPolyIndices() override
		{
			std::vector<std::vector<unsigned int>> polyIndices;
			for (const auto f : m_SurfaceMesh->faces())
			{
				std::vector<unsigned int> faceIndices;
				for (const auto v : m_SurfaceMesh->vertices(f))
				{
					faceIndices.push_back(v.idx());
				}
				polyIndices.push_back(faceIndices);
			}
			return polyIndices;
		}

		/// \copydoc MeshAdapter::IsTriangle
		[[nodiscard]] bool IsTriangle() const override
		{
			return m_SurfaceMesh->is_triangle_mesh();
		}

		/**
		 * \brief Retrieves the PMP surface mesh data.
		 * \return Reference to pmp::SurfaceMesh.
		 */
		pmp::SurfaceMesh& GetMesh()
		{
			return *m_SurfaceMesh;
		}

		/**
		 * \brief Retrieves the PMP surface mesh data (const version).
		 * \return Constant reference to pmp::SurfaceMesh.
		 */
		[[nodiscard]] const pmp::SurfaceMesh& GetMesh() const
		{
			return *m_SurfaceMesh;
		}

		/// \copydoc MeshAdapter::GetBounds
		[[nodiscard]] pmp::BoundingBox GetBounds() const override
		{
			return m_SurfaceMesh->bounds();
		}

	private:
		std::shared_ptr<pmp::SurfaceMesh> m_SurfaceMesh; //>! actual geometry instance
	};

	//
	// ==============================================================================
	//        2D Adapters

	/**
	 * \brief Abstract class that defines a common interface for curve adapters.
	 * \class CurveAdapter
	 */
	class CurveAdapter
	{
	public:
		/**
		 * \brief Creates a unique pointer to a new CurveAdapter that is a copy of the current object.
		 * \return A unique pointer to a new CurveAdapter.
		 */
		virtual [[nodiscard]] std::unique_ptr<CurveAdapter> Clone() const = 0;

		/// \brief Verifies whether this adapter applies to an empty curve.
		virtual [[nodiscard]] bool IsEmpty() const = 0;

		/**
		 * \brief Retrieves the vertices of the curve.
		 * \return A vector containing the vertices of the curve.
		 */
		virtual [[nodiscard]] std::vector<pmp::Point2> GetVertices() = 0;

		/**
		 * \brief Retrieves the polygon indices of the curve.
		 * \return A vector of vectors, where each inner vector contains the indices of a polygon.
		 */
		virtual [[nodiscard]] std::vector<std::pair<unsigned int, unsigned int>> GetEdgeIndices() = 0;

		/**
		 * \brief Retrieves the bounding box of the curve.
		 * \return A BoundingBox2 object representing the bounding box of the curve.
		 */
		virtual [[nodiscard]] pmp::BoundingBox2 GetBounds() const = 0;

		/**
		 * \brief Virtual destructor for CurveAdapter.
		 */
		virtual ~CurveAdapter() = default;
	};

	/**
	 * \brief Adapter class for BaseCurveGeometryData implementing CurveAdapter interface for 2D data.
	 * \class BaseCurveAdapter
	 */
	class BaseCurveAdapter : public CurveAdapter
	{
	public:
		/**
		 * \brief Constructor initializing the BaseCurveAdapter with given curve data.
		 * \param curve Shared pointer to BaseCurveGeometryData.
		 */
		explicit BaseCurveAdapter(std::shared_ptr<BaseCurveGeometryData> curve) : m_BaseCurve(std::move(curve)) {}

		/// \brief Verifies whether this adapter applies to an empty curve.
		[[nodiscard]] bool IsEmpty() const override
		{
			if (!m_BaseCurve)
				return true;

			return m_BaseCurve->Vertices.empty() || m_BaseCurve->EdgeIndices.empty();
		}

		/// \copydoc CurveAdapter::Clone
		[[nodiscard]] std::unique_ptr<CurveAdapter> Clone() const override
		{
			return std::make_unique<BaseCurveAdapter>(*this);
		}

		/// \copydoc CurveAdapter::GetVertices
		[[nodiscard]] std::vector<pmp::Point2> GetVertices() override
		{
			return m_BaseCurve->Vertices;
		}

		/// \copydoc CurveAdapter::GetEdgeIndices
		[[nodiscard]] std::vector<std::pair<unsigned int, unsigned int>> GetEdgeIndices() override
		{
			return m_BaseCurve->EdgeIndices;
		}

		/// \copydoc CurveAdapter::GetBounds
		[[nodiscard]] pmp::BoundingBox2 GetBounds() const override
		{
			pmp::BoundingBox2 bb;
			bb += m_BaseCurve->Vertices;
			return bb;
		}

	private:
		std::shared_ptr<BaseCurveGeometryData> m_BaseCurve;
	};

	/**
	 * \brief Adapter class for ManifoldCurve2D implementing CurveAdapter interface for 2D data.
	 * \class ManifoldCurve2DAdapter
	 */
	class ManifoldCurve2DAdapter : public CurveAdapter
	{
	public:
		/**
		 * \brief Constructor initializing the ManifoldCurve2DAdapter with given curve.
		 * \param curve Shared pointer to ManifoldCurve2D.
		 */
		explicit ManifoldCurve2DAdapter(std::shared_ptr<pmp::ManifoldCurve2D> curve) : m_ManifoldCurve(std::move(curve)) {}

		/// \brief Verifies whether this adapter applies to an empty curve.
		[[nodiscard]] bool IsEmpty() const override
		{
			if (!m_ManifoldCurve)
				return true;

			return m_ManifoldCurve->is_empty();
		}

		/// \copydoc CurveAdapter::Clone
		[[nodiscard]] std::unique_ptr<CurveAdapter> Clone() const override
		{
			return std::make_unique<ManifoldCurve2DAdapter>(*this);
		}

		/// \copydoc CurveAdapter::GetVertices
		[[nodiscard]] std::vector<pmp::Point2> GetVertices() override
		{
			std::vector<pmp::Point2> vertices;
			vertices.reserve(m_ManifoldCurve->n_vertices());
			for (const auto v : m_ManifoldCurve->vertices())
			{
				vertices.push_back(m_ManifoldCurve->position(v));
			}
			return vertices;
		}

		/// \copydoc CurveAdapter::GetEdgeIndices
		[[nodiscard]] std::vector<std::pair<unsigned int, unsigned int>> GetEdgeIndices() override
		{
			std::vector<std::pair<unsigned int, unsigned int>> edgeIndices;
			for (const auto e : m_ManifoldCurve->edges())
			{
				edgeIndices.emplace_back(m_ManifoldCurve->from_vertex(e).idx(), m_ManifoldCurve->to_vertex(e).idx());
			}
			return edgeIndices;
		}

		/// \copydoc CurveAdapter::GetBounds
		[[nodiscard]] pmp::BoundingBox2 GetBounds() const override
		{
			return m_ManifoldCurve->bounds();
		}

	private:
		std::shared_ptr<pmp::ManifoldCurve2D> m_ManifoldCurve;
	};

	
} // namespace Geometry
