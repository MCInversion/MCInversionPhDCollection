#include "pmp/ManifoldCurve2D.h"

namespace pmp
{

    void ManifoldCurve2D::clear()
    {
    }

    void ManifoldCurve2D::free_memory()
    {
    }

    void ManifoldCurve2D::reserve(size_t nvertices, size_t nedges)
    {
    }

    void ManifoldCurve2D::garbage_collection()
    {
    }

    ManifoldCurve2D::ManifoldCurve2D()
	{
        vpoint_ = add_vertex_property<Point2>("v:point");
        vconn_ = add_vertex_property<VertexConnectivity>("v:connectivity");
        econn_ = add_edge_property<EdgeConnectivity>("e:connectivity");
        
        vdeleted_ = add_vertex_property<bool>("v:deleted", false);
        edeleted_ = add_edge_property<bool>("e:deleted", false);
    }

    ManifoldCurve2D& ManifoldCurve2D::operator=(const ManifoldCurve2D& rhs)
    {
        if (this != &rhs)
        {
            // deep copy of property containers
            vprops_ = rhs.vprops_;
            eprops_ = rhs.eprops_;

            // property handles contain pointers, have to be reassigned
            vpoint_ = vertex_property<Point2>("v:point");
            vconn_ = vertex_property<VertexConnectivity>("v:connectivity");

            vdeleted_ = vertex_property<bool>("v:deleted");
            edeleted_ = edge_property<bool>("e:deleted");

            // how many elements are deleted?
            deleted_vertices_ = rhs.deleted_vertices_;
            deleted_edges_ = rhs.deleted_edges_;

            has_garbage_ = rhs.has_garbage_;
        }

        return *this;
    }

    ManifoldCurve2D& ManifoldCurve2D::operator*=(const mat3& mat)
    {
        for (const auto v : vertices())
            position(v) = affine_transform(mat, position(v));

        return *this;
    }

    Vertex ManifoldCurve2D::add_vertex(const Point2& p)
    {
        return Vertex();
    }

    Vertex ManifoldCurve2D::add_vertex_to_edge(Edge e, const Point2& p)
    {
        return Vertex();
    }

    Edge ManifoldCurve2D::add_edge(Vertex v0, Vertex v1)
    {
        return Edge();
    }

    void ManifoldCurve2D::delete_vertex(Vertex v)
    {
    }

    void ManifoldCurve2D::delete_edge(Edge e)
    {
    }

    void ManifoldCurve2D::collapse_edge(Vertex v0, Vertex v1)
    {
    }

    void ManifoldCurve2D::split_edge(Vertex v0, Vertex v1)
    {
    }

    void ManifoldCurve2D::remove_vertex(Vertex v)
    {
    }

} // namespace pmp
