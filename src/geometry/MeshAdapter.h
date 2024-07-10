#pragma once

#include "pmp/SurfaceMesh.h"
#include "GeometryConversionUtils.h"

namespace Geometry
{
	class MeshAdapter
	{
	public:
		virtual [[nodiscard]] std::unique_ptr<MeshAdapter> Clone() const = 0;
	    virtual [[nodiscard]] std::vector<pmp::vec3> GetVertices() = 0;
	    virtual [[nodiscard]] std::vector<std::vector<unsigned int>> GetPolyIndices() = 0;
		virtual [[nodiscard]] bool IsTriangle() const = 0;
		virtual [[nodiscard]] pmp::BoundingBox GetBounds() const = 0;
	    virtual ~MeshAdapter() = default;
	};

	class BaseMeshAdapter : public MeshAdapter
	{
	public:
	    explicit BaseMeshAdapter(std::shared_ptr<BaseMeshGeometryData> mesh) : m_BaseMesh(std::move(mesh)) {}

		[[nodiscard]] std::unique_ptr<MeshAdapter> Clone() const override
		{
			return std::make_unique<BaseMeshAdapter>(*this);
		}

		BaseMeshGeometryData& GetBaseMesh()
	    {
		    return *m_BaseMesh;
	    }

		[[nodiscard]] const BaseMeshGeometryData& GetBaseMesh() const
	    {
		    return *m_BaseMesh;
	    }

		[[nodiscard]] std::vector<pmp::vec3> GetVertices() override
		{
	        return m_BaseMesh->Vertices;
	    }

		[[nodiscard]] std::vector<std::vector<unsigned int>> GetPolyIndices() override
		{
			return m_BaseMesh->PolyIndices;
	    }

		[[nodiscard]] bool IsTriangle() const override
	    {
			return std::ranges::all_of(
				m_BaseMesh->PolyIndices, [](const auto& poly) { return poly.size() == 3; });
	    }

		[[nodiscard]] pmp::BoundingBox GetBounds() const override
	    {
			return pmp::BoundingBox(m_BaseMesh->Vertices);
	    }
	private:
	    std::shared_ptr<BaseMeshGeometryData> m_BaseMesh;
	};

	class PMPSurfaceMeshAdapter : public MeshAdapter
	{
	public:
		explicit PMPSurfaceMeshAdapter(std::shared_ptr<pmp::SurfaceMesh> mesh) : m_SurfaceMesh(std::move(mesh)) {}

		[[nodiscard]] std::unique_ptr<MeshAdapter> Clone() const override
		{
			return std::make_unique<PMPSurfaceMeshAdapter>(*this);
		}

		[[nodiscard]] std::vector<pmp::vec3> GetVertices() override
		{
			std::vector<pmp::vec3> vertices;
			vertices.reserve(m_SurfaceMesh->n_vertices());
			for (const auto v : m_SurfaceMesh->vertices())
			{
				vertices.push_back(m_SurfaceMesh->position(v));
			}
			return vertices;
		}

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

		[[nodiscard]] bool IsTriangle() const override
		{
			return m_SurfaceMesh->is_triangle_mesh();
		}

		pmp::SurfaceMesh& GetMesh()
		{
			return *m_SurfaceMesh;
		}

		[[nodiscard]] const pmp::SurfaceMesh& GetMesh() const
		{
			return *m_SurfaceMesh;
		}

		[[nodiscard]] pmp::BoundingBox GetBounds() const override
		{
			return m_SurfaceMesh->bounds();
		}
	private:

		std::shared_ptr<pmp::SurfaceMesh> m_SurfaceMesh;
	};
	
} // namespace Geometry
