#include "pmp/ManifoldCurve2D.h"

namespace pmp
{
    void ManifoldCurve2D::clear()
    {
        // remove all properties
        vprops_.clear();
        eprops_.clear();

        // really free their memory
        free_memory();

        // add the standard properties back
        vpoint_ = add_vertex_property<Point2>("v:point");
        vconn_ = add_vertex_property<VertexConnectivity>("v:connectivity");
        econn_ = add_edge_property<EdgeConnectivity>("e:connectivity");
        vdeleted_ = add_vertex_property<bool>("v:deleted", false);
        edeleted_ = add_edge_property<bool>("e:deleted", false);

        // set initial status (as in constructor)
        deleted_vertices_ = 0;
        deleted_edges_ = 0;
        has_garbage_ = false;
    }

    void ManifoldCurve2D::free_memory() const
    {
        vprops_.free_memory();
        eprops_.free_memory();
    }

    void ManifoldCurve2D::reserve(size_t nVertices, size_t nEdges) const
    {
        vprops_.reserve(nVertices);
        eprops_.reserve(nEdges);
    }

    void ManifoldCurve2D::garbage_collection()
    {
        if (!has_garbage_)
            return;

        auto nV = vertices_size();
        auto nE = edges_size();

        // setup handle mapping
        auto vmap = add_vertex_property<Vertex>("v:garbage-collection");
        auto emap = add_edge_property<Edge>("e:garbage-collection");
        for (size_t i = 0; i < nV; ++i)
            vmap[Vertex(i)] = Vertex(i);
        for (size_t i = 0; i < nE; ++i)
            emap[Edge(i)] = Edge(i);

        // remove deleted vertices
        if (nV > 0)
        {
            int i0 = 0;
            int i1 = nV - 1;

            while (true)
            {
                // find first deleted and last un-deleted
                while (!vdeleted_[Vertex(i0)] && i0 < i1)
                    ++i0;
                while (vdeleted_[Vertex(i1)] && i0 < i1)
                    --i1;
                if (i0 >= i1)
                    break;

                // swap
                vprops_.swap(i0, i1);
            }

            // remember new size
            nV = vdeleted_[Vertex(i0)] ? i0 : i0 + 1;
        }

        // remove deleted edges
        if (nE > 0)
        {
            int i0 = 0;
            int i1 = nE - 1;

            while (true)
            {
                // find first deleted and last un-deleted
                while (!edeleted_[Edge(i0)] && i0 < i1)
                    ++i0;
                while (edeleted_[Edge(i1)] && i0 < i1)
                    --i1;
                if (i0 >= i1)
                    break;

                // swap
                eprops_.swap(i0, i1);
            }

            // remember new size
            nE = edeleted_[Edge(i0)] ? i0 : i0 + 1;
        }

        // update vertex connectivity
        for (size_t i = 0; i < nV; ++i)
        {
            auto v = Vertex(i);
            if (is_isolated(v))
                continue;

            auto eTo = edge_to(v);
            auto eFrom = edge_from(v);
            if (eTo.is_valid())
            {
	            set_edge_to(v, emap[eTo]);
            }
            if (eFrom.is_valid())
            {
	            set_edge_from(v, emap[eFrom]);
            }
        }

        // update edge connectivity
        for (size_t i = 0; i < nE; ++i)
        {
            auto e = Edge(i);
            set_start_vertex(e, vmap[from_vertex(e)]);
            set_end_vertex(e, vmap[to_vertex(e)]);
        }

        // remove handle maps
        remove_vertex_property(vmap);
        remove_edge_property(emap);

        // finally resize arrays
        vprops_.resize(nV);
        vprops_.free_memory();
        eprops_.resize(nE);
        eprops_.free_memory();

        deleted_vertices_ = deleted_edges_ = 0;
        has_garbage_ = false;
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
            econn_ = edge_property<EdgeConnectivity>("e:connectivity");

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

    ManifoldCurve2D& ManifoldCurve2D::assign(const ManifoldCurve2D& rhs)
    {
        if (this != &rhs)
        {
            // clear properties
            vprops_.clear();
            eprops_.clear();

            // allocate standard properties
            vpoint_ = add_vertex_property<Point2>("v:point");
            vconn_ = add_vertex_property<VertexConnectivity>("v:connectivity");
            econn_ = add_edge_property<EdgeConnectivity>("e:connectivity");

            vdeleted_ = add_vertex_property<bool>("v:deleted", false);
            edeleted_ = add_edge_property<bool>("e:deleted", false);

            // copy properties from other curve
            vpoint_.array() = rhs.vpoint_.array();
            vconn_.array() = rhs.vconn_.array();
            econn_.array() = rhs.econn_.array();

            vdeleted_.array() = rhs.vdeleted_.array();
            edeleted_.array() = rhs.edeleted_.array();

            // resize (needed by property containers)
            vprops_.resize(rhs.vertices_size());
            eprops_.resize(rhs.edges_size());

            // how many elements are deleted?
            deleted_vertices_ = rhs.deleted_vertices_;
            deleted_edges_ = rhs.deleted_edges_;
            has_garbage_ = rhs.has_garbage_;
        }

        return *this;
    }

    Vertex ManifoldCurve2D::add_vertex(const Point2& p)
    {
        Vertex v = new_vertex();
        if (v.is_valid())
            vpoint_[v] = p;
        return v;
    }

    Edge ManifoldCurve2D::add_edge(Vertex v0, Vertex v1)
    {
        return new_edge(v0, v1);
    }

    void ManifoldCurve2D::delete_vertex(Vertex v, bool reconnect)
    {
        if (is_deleted(v))
            return;

        if (!is_isolated(v))
        {
	        const Edge eTo = vconn_[v].to_;
        	const Edge eFrom = vconn_[v].from_;

            if (!is_boundary(v) && reconnect)
            {
                // the first edge eT before v will remain connected to vertices vPrev and vNext
                const Vertex vNext = econn_[eFrom].end_;
                const Vertex vPrev = econn_[eTo].start_;
            	if (vNext.is_valid() && vPrev.is_valid())
	            {
	                vconn_[vNext].to_ = eTo;
                    vconn_[vPrev].from_ = eTo;
	            }
                if (eTo.is_valid())
                {
                    econn_[eTo].end_ = vNext;
                }
            }

            if (!reconnect)
            {
                // eTo is also cleared if reconnection is not desired
                edeleted_[eTo] = true;
                deleted_edges_++;
            }

            edeleted_[eFrom] = true;
            deleted_edges_++;
            has_garbage_ = true;
        }

        vdeleted_[v] = true;
        deleted_vertices_++;
        has_garbage_ = true;
    }

    void ManifoldCurve2D::delete_edge(Edge e)
    {
        if (is_deleted(e))
            return;

        const Vertex vStart = econn_[e].start_;
        const Vertex vEnd = econn_[e].end_;

        vconn_[vStart].from_.reset();
        vconn_[vEnd].to_.reset();

        edeleted_[e] = true;
        deleted_edges_++;
        has_garbage_ = true;
    }

    void ManifoldCurve2D::remove_edge(Edge e)
    {
        delete_edge(e);
    }

    void ManifoldCurve2D::remove_vertex(Vertex v, bool reconnect)
    {
        delete_vertex(v, reconnect);
    }

    void ManifoldCurve2D::collapse_edge(Edge e, bool keepStartVertex)
    {
        // before:
		//
		// ePrev    vFrom     e      vTo    eNext
		// ----------o--------------->o----------
		//
		// after (if keepStartVertex == true):
		//
        // ePrev    vFrom    eNext
        // ----------o-------------
        //
        // after (if keepStartVertex == false):
        //
        // ePrev    vTo    eNext
        // ---------->o----------
		//
    	const Vertex vFrom = econn_[e].start_;
        const Vertex vTo = econn_[e].end_;

        if (keepStartVertex)
        {
            const Edge eNext = vconn_[vTo].from_;
            vconn_[vFrom].from_ = eNext;
            if (is_valid(eNext))
            {
				econn_[eNext].start_ = vFrom;
            }

            vconn_[vTo].from_.reset();
            vconn_[vTo].to_.reset();
            vdeleted_[vTo] = true;
            edeleted_[e] = true;
            has_garbage_ = true;

            return;
        }

        const Edge ePrev = vconn_[vFrom].to_;
        vconn_[vTo].to_ = ePrev;
        if (is_valid(ePrev))
        {
            econn_[ePrev].end_ = vTo;
        }

        vconn_[vFrom].from_.reset();
        vconn_[vFrom].to_.reset();
        vdeleted_[vFrom] = true;
        edeleted_[e] = true;
        has_garbage_ = true;
    }

    void ManifoldCurve2D::collapse_edge(Vertex v0, Vertex v1, bool keepStartVertex)
    {
        // before:
        //
        // ePrev    v0     e         v1    eNext
        // ----------o--------------->o----------
        //
        // after (if keepStartVertex == true):
        //
        // ePrev    v0    eNext
        // ----------o-------------
        //
        // after (if keepStartVertex == false):
        //
        // ePrev    v1     eNext
        // ---------->o----------
        //
        const Edge e = vconn_[v0].from_;

        if (keepStartVertex)
        {
            const Edge eNext = vconn_[v0].from_;
            vconn_[v1].from_ = eNext;
            if (is_valid(eNext))
            {
                econn_[eNext].start_ = v1;
            }

            vconn_[v0].from_.reset();
            vconn_[v0].to_.reset();
            vdeleted_[v0] = true;
            edeleted_[e] = true;
            has_garbage_ = true;

            return;
        }

        const Edge ePrev = vconn_[v1].to_;
        vconn_[v0].to_ = ePrev;
        if (is_valid(ePrev))
        {
            econn_[ePrev].end_ = v1;
        }

        vconn_[v1].from_.reset();
        vconn_[v1].to_.reset();
        vdeleted_[v1] = true;
        edeleted_[e] = true;
        has_garbage_ = true;
    }

    Vertex ManifoldCurve2D::split_edge(Edge e, float param)
    {
        // before:
		//
		// v0      e        v1
		//  o--------------->o
		//
		// after:
		//
		// v0  e   v   eNew  v1
		//  o------>o------->o
        //

        const Vertex v0 = econn_[e].start_;
        const Vertex v1 = econn_[e].end_;

        const Point2 newPos = (1.0f - param) * position(v0) + param * position(v1);
        const Vertex v = add_vertex(newPos);

        econn_[e].end_ = v;
        vconn_[v1].to_.reset();
        const Edge eNew = new_edge(v, v1);
        if (!is_valid(eNew))
            return {};

        return v;
    }

    Vertex ManifoldCurve2D::split_edge(Vertex v0, Vertex v1, float param)
    {
        // before:
        //
        // v0      e        v1
        //  o--------------->o
        //
        // after:
        //
        // v0  e   v   eNew  v1
        //  o------>o------->o
        //

        const Edge e = vconn_[v0].from_;
        const Point2 newPos = (1.0f - param) * position(v0) + param * position(v1);
        const Vertex v = add_vertex(newPos);

        econn_[e].end_ = v;
        vconn_[v1].to_.reset();
        const Edge eNew = new_edge(v, v1);
        if (!is_valid(eNew))
            return {};

        return v;
    }

    Vertex ManifoldCurve2D::split_edge_with_vertex(Edge e, const Point2& pNew)
    {
        // before:
        //
        // v0      e        v1
        //  o--------------->o
        //
        // after:
        //
        // v0  e   v   eNew  v1
        //  o------>o------->o
        //
        const Vertex v1 = econn_[e].end_;
        const Vertex v = add_vertex(pNew);

        econn_[e].end_ = v;
        vconn_[v1].to_.reset();
        const Edge eNew = new_edge(v, v1);
        if (!is_valid(eNew))
            return {};

        return v;
    }

    Vertex ManifoldCurve2D::split_edge_with_vertex(Vertex v0, Vertex v1, const Point2& pNew)
    {
        // before:
		//
		// v0      e        v1
		//  o--------------->o
		//
		// after:
		//
		// v0  e   v   eNew  v1
		//  o------>o------->o
		//
        const Edge e = vconn_[v0].from_;
        const Vertex v = add_vertex(pNew);

        econn_[e].end_ = v;
        vconn_[v1].to_.reset();
        const Edge eNew = new_edge(v, v1);
        if (!is_valid(eNew))
            return {};

        return v;
    }

} // namespace pmp
