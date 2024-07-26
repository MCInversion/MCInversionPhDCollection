#pragma once

#include <vector>

#include "pmp/Types.h"
#include "pmp/Properties.h"

namespace pmp
{

    //! \addtogroup core
    //!@{

    // Handle Types

    //! Base class for all entity handles types.
    //! \details internally it is basically an index.
    class Handle
    {
    public:
        //! default constructor with invalid index
        explicit Handle(IndexType idx = PMP_MAX_INDEX) : idx_(idx) {}

        //! Get the underlying index of this handle
        [[nodiscard]] IndexType idx() const { return idx_; }

        //! reset handle to be invalid (index=PMP_MAX_INDEX.)
        void reset() { idx_ = PMP_MAX_INDEX; }

        //! \return whether the handle is valid, i.e., the index is not equal to PMP_MAX_INDEX.
        [[nodiscard]] bool is_valid() const { return idx_ != PMP_MAX_INDEX; }

        //! are two handles equal?
        bool operator==(const Handle& rhs) const { return idx_ == rhs.idx_; }

        //! are two handles different?
        bool operator!=(const Handle& rhs) const { return idx_ != rhs.idx_; }

        //! compare operator useful for sorting handles
        bool operator<(const Handle& rhs) const { return idx_ < rhs.idx_; }

    private:
        // =========== !!!!! ===================
        // Supported geometry types
        // -------------------------------------
        friend class SurfaceMesh;
        friend class ManifoldCurve2D;
        // -------------------------------------
        IndexType idx_;
    };

    //! hash function for handles.
    class HandleHashFunction
    {
    public:
        size_t operator()(const Handle& h) const { return std::hash<IndexType>()(h.idx()); }
    };


    //! this type represents a vertex (internally it is basically an index)
    class Vertex : public Handle
    {
        using Handle::Handle;
    };

    //! this type represents a halfedge (internally it is basically an index)
    class Halfedge : public Handle
    {
        using Handle::Handle;
    };

    //! this type represents an edge (internally it is basically an index)
    class Edge : public Handle
    {
        using Handle::Handle;
    };

    //! this type represents a face (internally it is basically an index)
    class Face : public Handle
    {
        using Handle::Handle;
    };

    // Output operators

    //! output a Vertex to a stream
    inline std::ostream& operator<<(std::ostream& os, Vertex v)
    {
        return (os << 'v' << v.idx());
    }

    //! output a Halfedge to a stream
    inline std::ostream& operator<<(std::ostream& os, Halfedge h)
    {
        return (os << 'h' << h.idx());
    }

    //! output an Edge to a stream
    inline std::ostream& operator<<(std::ostream& os, Edge e)
    {
        return (os << 'e' << e.idx());
    }

    //! output a Face to a stream
    inline std::ostream& operator<<(std::ostream& os, Face f)
    {
        return (os << 'f' << f.idx());
    }

    // Property Types

    //! Vertex property of type T
    template <class T>
    class VertexProperty : public Property<T>
    {
    public:
        //! default constructor
        explicit VertexProperty() = default;
        explicit VertexProperty(Property<T> p) : Property<T>(p) {}

        //! access the data stored for vertex \p v
        typename Property<T>::reference operator[](Vertex v)
        {
            return Property<T>::operator[](v.idx());
        }

        //! access the data stored for vertex \p v
        typename Property<T>::const_reference operator[](Vertex v) const
        {
            return Property<T>::operator[](v.idx());
        }
    };

    //! Halfedge property of type T
    template <class T>
    class HalfedgeProperty : public Property<T>
    {
    public:
        //! default constructor
        explicit HalfedgeProperty() = default;
        explicit HalfedgeProperty(Property<T> p) : Property<T>(p) {}

        //! access the data stored for halfedge \p h
        typename Property<T>::reference operator[](Halfedge h)
        {
            return Property<T>::operator[](h.idx());
        }

        //! access the data stored for halfedge \p h
        typename Property<T>::const_reference operator[](Halfedge h) const
        {
            return Property<T>::operator[](h.idx());
        }
    };

    //! Edge property of type T
    template <class T>
    class EdgeProperty : public Property<T>
    {
    public:
        //! default constructor
        explicit EdgeProperty() = default;
        explicit EdgeProperty(Property<T> p) : Property<T>(p) {}

        //! access the data stored for edge \p e
        typename Property<T>::reference operator[](Edge e)
        {
            return Property<T>::operator[](e.idx());
        }

        //! access the data stored for edge \p e
        typename Property<T>::const_reference operator[](Edge e) const
        {
            return Property<T>::operator[](e.idx());
        }
    };

    //! Face property of type T
    template <class T>
    class FaceProperty : public Property<T>
    {
    public:
        //! default constructor
        explicit FaceProperty() = default;
        explicit FaceProperty(Property<T> p) : Property<T>(p) {}

        //! access the data stored for face \p f
        typename Property<T>::reference operator[](Face f)
        {
            return Property<T>::operator[](f.idx());
        }

        //! access the data stored for face \p f
        typename Property<T>::const_reference operator[](Face f) const
        {
            return Property<T>::operator[](f.idx());
        }
    };

    //! Object property of type T
    template <class T>
    class ObjectProperty : public Property<T>
    {
    public:
        //! default constructor
        explicit ObjectProperty() = default;
        explicit ObjectProperty(Property<T> p) : Property<T>(p) {}

        //! access the data stored for the object
        typename Property<T>::reference operator[](IndexType idx)
        {
            return Property<T>::operator[](idx);
        }

        //! access the data stored for the object
        typename Property<T>::const_reference operator[](IndexType idx) const
        {
            return Property<T>::operator[](idx);
        }
    };

} // namespace pmp