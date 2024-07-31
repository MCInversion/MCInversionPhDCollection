#pragma once

#include <vector>

#include "pmp/Types.h"
#include "pmp/HandleAndPropertyTypes.h"
#include "pmp/BoundingBox.h"

namespace pmp
{
    //! A data structure for discretized manifold curves
    class ManifoldCurve2D
	{
    public:
        //! \name Iterator Types
		//!@{

		//! An iterator class to iterate linearly over all vertices
        class VertexIterator
        {
        public:
            using difference_type = std::ptrdiff_t;
            using value_type = Vertex;
            using reference = Vertex&;
            using pointer = Vertex*;
            using iterator_category = std::bidirectional_iterator_tag;

            //! Default constructor
            VertexIterator(Vertex v = Vertex(), const ManifoldCurve2D* c = nullptr)
                : handle_(v), curve_(c)
            {
                if (curve_ && curve_->has_garbage())
                    while (curve_->is_valid(handle_) && curve_->is_deleted(handle_))
                        ++handle_.idx_;
            }

            //! get the vertex the iterator refers to
            Vertex operator*() const { return handle_; }

            //! are two iterators equal?
            bool operator==(const VertexIterator& rhs) const
            {
                return (handle_ == rhs.handle_);
            }

            //! are two iterators different?
            bool operator!=(const VertexIterator& rhs) const
            {
                return !operator==(rhs);
            }

            //! pre-increment iterator
            VertexIterator& operator++()
            {
                ++handle_.idx_;
                assert(curve_);
                while (curve_->has_garbage() && curve_->is_valid(handle_) &&
                    curve_->is_deleted(handle_))
                    ++handle_.idx_;
                return *this;
            }

            //! post-increment iterator
            VertexIterator operator++(int)
            {
                auto tmp = *this;
                ++(*this);
                return tmp;
            }

            //! pre-decrement iterator
            VertexIterator& operator--()
            {
                --handle_.idx_;
                assert(curve_);
                while (curve_->has_garbage() && curve_->is_valid(handle_) &&
                    curve_->is_deleted(handle_))
                    --handle_.idx_;
                return *this;
            }

            //! post-decrement iterator
            VertexIterator operator--(int)
            {
                auto tmp = *this;
                --(*this);
                return tmp;
            }

            //! add offset to iterator
            VertexIterator operator+(int i) const
            {
                VertexIterator it = *this;
                if (i > 0)
                {
                    while (i > 0)
                    {
                        ++it.handle_.idx_;
                        while (it.curve_->has_garbage() && it.curve_->is_valid(it.handle_) &&
                            it.curve_->is_deleted(it.handle_))
                            ++it.handle_.idx_;
                        --i;
                    }
                }
                else
                {
                    while (i < 0)
                    {
                        --it.handle_.idx_;
                        while (it.curve_->has_garbage() && it.curve_->is_valid(it.handle_) &&
                            it.curve_->is_deleted(it.handle_))
                            --it.handle_.idx_;
                        ++i;
                    }
                }
                return it;
            }

            //! subtract offset from iterator
            VertexIterator operator-(int i) const
            {
                return operator+(-i);
            }

        private:
            Vertex handle_;
            const ManifoldCurve2D* curve_;
        };

        //! this class iterates linearly over all edges
        //! \sa edges_begin(), edges_end()
        //! \sa VertexIterator, HalfedgeIterator, FaceIterator
        class EdgeIterator
        {
        public:
            using difference_type = std::ptrdiff_t;
            using value_type = Edge;
            using reference = Edge&;
            using pointer = Edge*;
            using iterator_category = std::bidirectional_iterator_tag;

            //! Default constructor
            EdgeIterator(Edge e = Edge(), const ManifoldCurve2D* c = nullptr)
                : handle_(e), curve_(c)
            {
                if (curve_ && curve_->has_garbage())
                    while (curve_->is_valid(handle_) && curve_->is_deleted(handle_))
                        ++handle_.idx_;
            }

            //! get the edge the iterator refers to
            Edge operator*() const { return handle_; }

            //! are two iterators equal?
            bool operator==(const EdgeIterator& rhs) const
            {
                return (handle_ == rhs.handle_);
            }

            //! are two iterators different?
            bool operator!=(const EdgeIterator& rhs) const
            {
                return !operator==(rhs);
            }

            //! pre-increment iterator
            EdgeIterator& operator++()
            {
                ++handle_.idx_;
                assert(curve_);
                while (curve_->has_garbage() && curve_->is_valid(handle_) &&
                    curve_->is_deleted(handle_))
                    ++handle_.idx_;
                return *this;
            }

            //! post-increment iterator
            EdgeIterator operator++(int)
            {
                auto tmp = *this;
                ++(*this);
                return tmp;
            }

            //! pre-decrement iterator
            EdgeIterator& operator--()
            {
                --handle_.idx_;
                assert(curve_);
                while (curve_->has_garbage() && curve_->is_valid(handle_) &&
                    curve_->is_deleted(handle_))
                    --handle_.idx_;
                return *this;
            }

            //! post-decrement iterator
            EdgeIterator operator--(int)
            {
                auto tmp = *this;
                --(*this);
                return tmp;
            }

        private:
            Edge handle_;
            const ManifoldCurve2D* curve_;
        };

        //! helper class for iterating through all vertices using range-based
		//! for-loops.
        class VertexContainer
        {
        public:
            VertexContainer(VertexIterator begin, VertexIterator end)
                : begin_(begin), end_(end)
            {
            }
            [[nodiscard]] VertexIterator begin() const { return begin_; }
            [[nodiscard]] VertexIterator end() const { return end_; }

        private:
            VertexIterator begin_;
            VertexIterator end_;
        };

        //! helper class for iterating through all edges using range-based
        //! for-loops. \sa edges()
        class EdgeContainer
        {
        public:
            EdgeContainer(EdgeIterator begin, EdgeIterator end)
                : begin_(begin), end_(end)
            {
            }
            [[nodiscard]] EdgeIterator begin() const { return begin_; }
            [[nodiscard]] EdgeIterator end() const { return end_; }

        private:
            EdgeIterator begin_;
            EdgeIterator end_;
        };

        //!@}
		//! \name Memory Management
		//!@{

    	//! \return number of (deleted and valid) vertices in the curve
    	[[nodiscard]] size_t vertices_size() const { return vprops_.size(); }

        //! \return number of (deleted and valid) edges in the curve
        [[nodiscard]] size_t edges_size() const { return eprops_.size(); }

        //! \return number of vertices in the curve
        [[nodiscard]] size_t n_vertices() const { return vertices_size() - deleted_vertices_; }

        //! \return number of edges in the curve
        [[nodiscard]] size_t n_edges() const { return edges_size() - deleted_edges_; }

        //! \return true if the curve is empty, i.e., has no vertices
        [[nodiscard]] bool is_empty() const { return n_vertices() == 0; }

        //! clear curve: remove all vertices, edges, faces
        virtual void clear();

        //! remove unused memory from vectors
        void free_memory() const;

        //! reserve memory (mainly used in file readers)
        void reserve(size_t nVertices, size_t nEdges) const;

        //! remove deleted elements
        void garbage_collection();

        // are there any deleted entities?
        [[nodiscard]] bool has_garbage() const { return has_garbage_; }

        //
        // ============== Constructors and Destructor =======================
        //
        ManifoldCurve2D();

        virtual ~ManifoldCurve2D() = default;

        //! copy constructor: copies \p rhs to \p *this. performs a deep copy of all
		//! properties.
        ManifoldCurve2D(const ManifoldCurve2D& rhs) { operator=(rhs); }

        // ================ Operators =======================================

        //! assign \p rhs to \p *this. performs a deep copy of all properties.
        ManifoldCurve2D& operator=(const ManifoldCurve2D& rhs);

        //! transform all points of this surface curve with a given 4x4 matrix.
        ManifoldCurve2D& operator*=(const mat3& mat);

        //! assign \p rhs to \p *this. does not copy custom properties.
        ManifoldCurve2D& assign(const ManifoldCurve2D& rhs);

        // ================== Properties ====================================

        //! add a vertex property of type \p T with name \p name and default
		//! value \p t. fails if a property named \p name exists already,
		//! since the name has to be unique. in this case it returns an
		//! invalid property
        template <class T>
        VertexProperty<T> add_vertex_property(const std::string& name, const T t = T())
        {
            return VertexProperty<T>(vprops_.add<T>(name, t));
        }

        //! get the vertex property named \p name of type \p T. returns an
        //! invalid VertexProperty if the property does not exist or if the
        //! type does not match.
        template <class T>
        VertexProperty<T> get_vertex_property(const std::string& name) const
        {
            return VertexProperty<T>(vprops_.get<T>(name));
        }

        //! if a vertex property of type \p T with name \p name exists, it is
		//! returned. otherwise this property is added (with default value \c
		//! t)
        template <class T>
        VertexProperty<T> vertex_property(const std::string& name, const T t = T())
        {
            return VertexProperty<T>(vprops_.get_or_add<T>(name, t));
        }

        //! remove the vertex property \p p
        template <class T>
        void remove_vertex_property(VertexProperty<T>& p)
        {
            vprops_.remove(p);
        }

        //! does the curve have a vertex property with name \p name?
        [[nodiscard]] bool has_vertex_property(const std::string& name) const
        {
            return vprops_.exists(name);
        }

        //! add a edge property of type \p T with name \p name and default
	    //! value \p t.  fails if a property named \p name exists already,
	    //! since the name has to be unique.  in this case it returns an
	    //! invalid property.
        template <class T>
        EdgeProperty<T> add_edge_property(const std::string& name, const T t = T())
        {
            return EdgeProperty<T>(eprops_.add<T>(name, t));
        }

        //! get the edge property named \p name of type \p T. returns an
		//! invalid VertexProperty if the property does not exist or if the
		//! type does not match.
        template <class T>
        EdgeProperty<T> get_edge_property(const std::string& name) const
        {
            return EdgeProperty<T>(eprops_.get<T>(name));
        }

        //! if an edge property of type \p T with name \p name exists, it is
		//! returned.  otherwise this property is added (with default value \c
		//! t)
        template <class T>
        EdgeProperty<T> edge_property(const std::string& name, const T t = T())
        {
            return EdgeProperty<T>(eprops_.get_or_add<T>(name, t));
        }

        //! remove the edge property \p p
        template <class T>
        void remove_edge_property(EdgeProperty<T>& p)
        {
            eprops_.remove(p);
        }

        //! does the curve have an edge property with name \p name?
        [[nodiscard]] bool has_edge_property(const std::string& name) const
        {
            return eprops_.exists(name);
        }

        //!@}
        //! \name Low-level connectivity
        //!@{

        //! set the outgoing edge of vertex \p v to \p e
        void set_edge_from(Vertex v, Edge e) { vconn_[v].from_ = e; }
        //! set the incoming edge of vertex \p v to \p e
        void set_edge_to(Vertex v, Edge e) { vconn_[v].to_ = e; }

        //! \return whether \p v is a boundary vertex
        [[nodiscard]] bool is_boundary(Vertex v) const
        {
            return (!vconn_[v].from_.is_valid() && vconn_[v].to_.is_valid()) ||
                (!vconn_[v].to_.is_valid() && vconn_[v].from_.is_valid());
        }

        //! \return whether \p v is isolated, i.e., not incident to any edge
        [[nodiscard]] bool is_isolated(Vertex v) const
        {
	        return !vconn_[v].from_.is_valid() && !vconn_[v].to_.is_valid();
        }

        //! sets the start vertex the edge \p e to \p v
        void set_start_vertex(Edge e, Vertex v) { econn_[e].start_ = v; }
        //! sets the end vertex the edge \p e to \p v
        void set_end_vertex(Edge e, Vertex v) { econn_[e].end_ = v; }

        //!@}

		//! \name Allocate new elements
		//!@{

		//! \brief Allocate a new vertex, resize vertex properties accordingly.
		//! \throw AllocationException in case of failure to allocate a new vertex.
        Vertex new_vertex()
        {
            if (vertices_size() == PMP_MAX_INDEX - 1)
            {
                auto what =
                    "ManifoldCurve2D: cannot allocate vertex, max. index reached";
                throw AllocationException(what);
            }
            vprops_.push_back();
            return Vertex(static_cast<IndexType>(vertices_size()) - 1);
        }

        //! \brief Allocate a new edge, resize edge property accordingly.
        //! \throw AllocationException in case of failure to allocate a new edge.
        Edge new_edge()
        {
            if (edges_size() == PMP_MAX_INDEX - 1)
            {
                auto what = "ManifoldCurve2D: cannot allocate edge, max. index reached";
                throw AllocationException(what);
            }
            eprops_.push_back();
            return Edge(static_cast<IndexType>(edges_size()) - 1);
        }

        //! \brief Allocate a new edge, resize edge property accordingly.
		//! \throw AllocationException in case of failure to allocate a new edge.
		//! \param start starting Vertex of the new edge
		//! \param end end Vertex of the new edge
        //! \warning this function does not check for non-manifoldness.
        Edge new_edge(Vertex start, Vertex end)
        {
            assert(start != end);

            if (edges_size() == PMP_MAX_INDEX - 1)
            {
                auto what = "ManifoldCurve2D: cannot allocate edge, max. index reached!\n";
                throw AllocationException(what);
            }

            eprops_.push_back();
            const Edge newEdge = Edge(static_cast<IndexType>(edges_size()) - 1);
            set_start_vertex(newEdge, start);
            set_end_vertex(newEdge, end);
            set_edge_from(start, newEdge);
            set_edge_to(end, newEdge);
            return newEdge;
        }

        //! \return whether vertex \p v is deleted
		//! \sa garbage_collection()
        [[nodiscard]] bool is_deleted(Vertex v) const { return vdeleted_[v]; }

        //! \return whether edge \p e is deleted
        //! \sa garbage_collection()
        [[nodiscard]] bool is_deleted(Edge e) const { return edeleted_[e]; }

        //! \return whether vertex \p v is valid.
        [[nodiscard]] bool is_valid(Vertex v) const { return v.idx() < vertices_size(); }

        //! \return whether edge \p e is valid.
        [[nodiscard]] bool is_valid(Edge e) const { return e.idx() < edges_size(); }

        //! deletes the vertex \p v from the curve. If \p reconnect is true, the adjacent edges of v won't be discarded.
        void delete_vertex(Vertex v, bool reconnect = true);

        //! deletes the edge \p e from the curve
        void delete_edge(Edge e);

        //! deletes the edge \p e from the curve
        void remove_edge(Edge e);

        //! deletes the vertex \p v from the curve. If \p reconnect is true, the adjacent edges of v won't be discarded.
        void remove_vertex(Vertex v, bool reconnect = true);

        //! adds a new vertex with position \p p, but does no connectivity adjustments.
        Vertex add_vertex(const Point2& p);

        //! adds a new edge between vertices \p v0 and \p v1.
        //! \throw TopologyException if some of the two endpoints are already connected to two edges.
        Edge add_edge(Vertex start, Vertex end);

        //
        // ========= Vertices and Edges Management (Tessellation changing operations)
        //

        //! Collapses edge \p e to a single vertex (start if \p keepStartVertex == true).
        void collapse_edge(Edge e, bool keepStartVertex = true);

        //! Collapses edge spanning \p v0 and \p v1 to a single vertex (start if \p keepStartVertex == true).
        void collapse_edge(Vertex v0, Vertex v1, bool keepStartVertex = true);

        //! Splits (subdivides) edge \p e with a new vertex lin. interpolated between its endpoints with \p param.
        Vertex split_edge(Edge e, float param = 0.5f);

        //! Splits (subdivides) edge with a new vertex lin. interpolated between its endpoints \p v0 and \p v1 with \p param.
        Vertex split_edge(Vertex v0, Vertex v1, float param = 0.5f);

        //! Splits (subdivides) edge \p e with a new vertex given by position \p pNew.
        Vertex split_edge_with_vertex(Edge e, const Point2& pNew);

        //! Splits (subdivides) edge between vertices \p v0 and \p v1 with a new vertex given by position \p pNew.
        Vertex split_edge_with_vertex(Vertex v0, Vertex v1, const Point2& pNew);

        //!@}
		//! \name Iterators
		//!@{

		//! \return start iterator for vertices
        [[nodiscard]] VertexIterator vertices_begin() const
        {
            return { Vertex(0), this };
        }

        //! \return end iterator for vertices
        [[nodiscard]] VertexIterator vertices_end() const
        {
            return { Vertex(static_cast<IndexType>(vertices_size())),this };
        }

        //! \return vertex container for C++11 range-based for-loops
        [[nodiscard]] VertexContainer vertices() const
        {
            return { vertices_begin(), vertices_end() };
        }

        //! \return start iterator for edges
        [[nodiscard]] EdgeIterator edges_begin() const
        {
	        return { Edge(0), this };
        }

        //! \return end iterator for edges
        [[nodiscard]] EdgeIterator edges_end() const
        {
            return { Edge(static_cast<IndexType>(edges_size())), this };
        }

        //! \return edge container for C++11 range-based for-loops
        [[nodiscard]] EdgeContainer edges() const
        {
            return { edges_begin(), edges_end() };
        }

        // ============ Neighborhood operations =========

        //! \return the edge that contains vertex \p v as a starting vertex
        [[nodiscard]] Edge edge_from(Vertex v) const
    	{
            if (!is_valid(v))
                return {};
            return vconn_[v].from_;
        }

        //! \return the edge that contains vertex \p v as the end vertex
        [[nodiscard]] Edge edge_to(Vertex v) const
        {
            if (!is_valid(v))
                return {};
            return vconn_[v].to_;
        }

        //! \return the vertex that edge \p e starts at
        [[nodiscard]] Vertex from_vertex(Edge e) const
        {
            if (!is_valid(e))
                return {};
            return econn_[e].start_;
        }

        //! \return the vertex that edge \p e ends at
        [[nodiscard]] Vertex to_vertex(Edge e) const
        {
            if (!is_valid(e))
                return {};
            return econn_[e].end_;
        }


        //! \return edges adjacent to vertex \p v.
        [[nodiscard]] std::pair<Edge, Edge> edges(Vertex v) const
        {
            if (!is_valid(v))
                return {};
            return { vconn_[v].to_, vconn_[v].from_ };
        }

        //! \return vertices connected to vertex \p v via edges.
        [[nodiscard]] std::pair<Vertex, Vertex> vertices(Vertex v) const
        {
            if (!is_valid(v))
                return {};
            const auto toEdge = edge_to(v);
            const auto fromEdge = edge_from(v);
            return { from_vertex(toEdge), to_vertex(fromEdge) };
        }

        //! \return start and end vertices on edge \p e.
        [[nodiscard]] std::pair<Vertex, Vertex> vertices(Edge e) const
        {
            if (!is_valid(e))
                return {};
            return { econn_[e].start_, econn_[e].end_ };
        }

        //!@}
        //! \name Geometry-related Functions
        //!@{

        //! position of a vertex (read only)
        [[nodiscard]] const Point2& position(Vertex v) const { return vpoint_[v]; }

    	//! position of a vertex
        Point2& position(Vertex v) { return vpoint_[v]; }

        //! \return vector of point positions
        std::vector<Point2>& positions() { return vpoint_.vector(); }

        //! compute the bounding box of the object
        [[nodiscard]] BoundingBox2 bounds() const
        {
            BoundingBox2 bb;
            for (auto v : vertices())
                bb += position(v);
            return bb;
        }

        //! compute the length of edge \p e.
        [[nodiscard]] Scalar edge_length(Edge e) const
        {
            return norm(vpoint_[from_vertex(e)] - vpoint_[to_vertex(e)]);
        }

        //! compute the squared of edge \p e.
        [[nodiscard]] Scalar edge_length_sq(Edge e) const
        {
            return dot(vpoint_[from_vertex(e)] - vpoint_[to_vertex(e)], vpoint_[from_vertex(e)] - vpoint_[to_vertex(e)]);
        }

    private:

        struct VertexConnectivity
        {
            Edge to_; // incoming edge to the given vertex
            Edge from_; // outgoing edge from the given vertex
        };

        struct EdgeConnectivity
        {
            Vertex start_; // base (first) vertex of the given edge
            Vertex end_; // target (second) vertex of the given edge
        };

        // property containers for each entity type and object
        PropertyContainer vprops_;
        PropertyContainer eprops_;

		// point coordinates
        VertexProperty<Point2> vpoint_;

        // connectivity information
        VertexProperty<VertexConnectivity> vconn_;
        EdgeProperty<EdgeConnectivity> econn_;

        // markers for deleted entities
        VertexProperty<bool> vdeleted_;
        EdgeProperty<bool> edeleted_;

        // numbers of deleted entities
        IndexType deleted_vertices_{ 0 };
        IndexType deleted_edges_{ 0 };

        // indicate garbage present
        bool has_garbage_{ false };

        std::string name_;
    };

    [[nodiscard]] bool write_to_ply(const pmp::ManifoldCurve2D& curve, const std::string& fileName);

} // namespace pmp
