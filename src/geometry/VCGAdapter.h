#pragma once

#include <vcg/complex/complex.h>
#include <vcg/simplex/face/component.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/create/ball_pivoting.h>
#include <vcg/complex/algorithms/convex_hull.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>

#include "pmp/Types.h"

// VCG mesh types
typedef pmp::Scalar Scalarm;
typedef vcg::Point2<pmp::Scalar>     Point2m;
typedef vcg::Point3<pmp::Scalar>     Point3m;
typedef vcg::Point4<pmp::Scalar>     Point4m;
typedef vcg::Plane3<pmp::Scalar>     Plane3m;
typedef vcg::Segment2<pmp::Scalar>   Segment2m;
typedef vcg::Segment3<pmp::Scalar>   Segment3m;
typedef vcg::Box3<pmp::Scalar>       Box3m;
typedef vcg::Matrix44<pmp::Scalar>   Matrix44m;
typedef vcg::Matrix33<pmp::Scalar>   Matrix33m;
typedef vcg::Shot<pmp::Scalar>       Shotm;
typedef vcg::Similarity<pmp::Scalar> Similaritym;

namespace vcg
{
	namespace vertex
	{
		template <class T> class Coord3m : public Coord<vcg::Point3<Scalarm>, T> {
		public:	static void Name(std::vector<std::string>& name) { name.push_back(std::string("Coord3m")); T::Name(name); }
		};

		template <class T> class Normal3m : public Normal<vcg::Point3<Scalarm>, T> {
		public:	static void Name(std::vector<std::string>& name) { name.push_back(std::string("Normal3m")); T::Name(name); }
		};

		template <class T> class Qualitym : public Quality<Scalarm, T> {
		public: static void Name(std::vector<std::string>& name) { name.push_back(std::string("Qualitym")); T::Name(name); }
		};

		template <class T> class CurvatureDirmOcf : public CurvatureDirOcf<CurvatureDirTypeOcf<Scalarm>, T> {
		public:	static void Name(std::vector<std::string>& name) { name.push_back(std::string("CurvatureDirmOcf")); T::Name(name); }
		};

		template <class T> class RadiusmOcf : public RadiusOcf<Scalarm, T> {
		public:	static void Name(std::vector<std::string>& name) { name.push_back(std::string("RadiusmOcf")); T::Name(name); }
		};

	}//end namespace vertex

	namespace face
	{
		template <class T> class Normal3m : public NormalAbs<vcg::Point3<Scalarm>, T> {
		public:  static void Name(std::vector<std::string>& name) { name.push_back(std::string("Normal3m")); T::Name(name); }
		};

		template <class T> class QualitymOcf : public QualityOcf<Scalarm, T> {
		public:  static void Name(std::vector<std::string>& name) { name.push_back(std::string("QualitymOcf")); T::Name(name); }
		};

		template <class T> class CurvatureDirmOcf : public CurvatureDirOcf<CurvatureDirOcfBaseType<Scalarm>, T> {
		public:	static void Name(std::vector<std::string>& name) { name.push_back(std::string("CurvatureDirdOcf")); T::Name(name); }
		};

	}//end namespace face
}//end namespace vcg

class VCG_Vertex; class VCG_Edge; class VCG_Face;

struct VCG_UsedTypes : public vcg::UsedTypes<
	vcg::Use<VCG_Vertex>   ::AsVertexType,
	vcg::Use<VCG_Edge>     ::AsEdgeType,
	vcg::Use<VCG_Face>     ::AsFaceType> {};

class VCG_Vertex : public vcg::Vertex< VCG_UsedTypes,
	vcg::vertex::InfoOcf,           /*  4b */
	vcg::vertex::Coord3m,           /* 12b */
	vcg::vertex::BitFlags,          /*  4b */
	vcg::vertex::Normal3m,          /* 12b */
	vcg::vertex::Qualitym,          /*  4b */
	vcg::vertex::Color4b,           /*  4b */
	vcg::vertex::VFAdjOcf,          /*  0b */
	vcg::vertex::MarkOcf,           /*  0b */
	vcg::vertex::TexCoordfOcf,      /*  0b */
	vcg::vertex::CurvatureDirmOcf,  /*  0b */
	vcg::vertex::RadiusmOcf         /*  0b */  > {};
class VCG_Face : public vcg::Face<   VCG_UsedTypes,
	vcg::face::InfoOcf,              /* 4b */
	vcg::face::VertexRef,            /*12b */
	vcg::face::BitFlags,             /* 4b */
	vcg::face::Normal3m,             /*12b */
	vcg::face::QualitymOcf,          /* 0b */
	vcg::face::MarkOcf,              /* 0b */
	vcg::face::Color4bOcf,           /* 0b */
	vcg::face::FFAdjOcf,             /* 0b */
	vcg::face::VFAdjOcf,             /* 0b */
	vcg::face::CurvatureDirmOcf,     /* 0b */
	vcg::face::WedgeTexCoordfOcf     /* 0b */ > {};
class VCG_Edge : public vcg::Edge<   VCG_UsedTypes,
	vcg::edge::BitFlags,          /*  4b */
	vcg::edge::EVAdj,
	vcg::edge::EEAdj> {};

class VCG_Mesh : public vcg::tri::TriMesh< vcg::vertex::vector_ocf<VCG_Vertex>, vcg::face::vector_ocf<VCG_Face> > {};

class MyVertex0 : public vcg::Vertex< VCG_UsedTypes, vcg::vertex::Coord3f, vcg::vertex::BitFlags  > {};
class MyVertex1 : public vcg::Vertex< VCG_UsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::BitFlags  > {};
class MyVertex2 : public vcg::Vertex< VCG_UsedTypes, vcg::vertex::Coord3f, vcg::vertex::Color4b, vcg::vertex::CurvatureDirf,
	vcg::vertex::Qualityf, vcg::vertex::Normal3f, vcg::vertex::BitFlags  > {};