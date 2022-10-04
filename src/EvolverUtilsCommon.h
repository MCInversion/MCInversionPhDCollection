#pragma once

#include "pmp/MatVec.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <functional>

namespace pmp
{
	class SurfaceMesh;
	class Vertex;
}

/// \brief a stats wrapper for co-volume measures affecting the stability of the finite volume method.
struct CoVolumeStats
{
	double Mean{ 0.0 };
	double Max{ -DBL_MAX };
	double Min{ DBL_MAX };
};

const std::string coVolMeasureVertexPropertyName{ "v:coVolumeMeasure" };

/// \brief a co-volume area evaluator.
using AreaFunction = std::function<double(const pmp::SurfaceMesh&, pmp::Vertex)>;

/**
 * \brief Analyzes the stats of co-volumes around each mesh vertex, and creates a vertex property for the measure values.
 * \param mesh             input mesh.
 * \param areaFunction     function to evaluate vertex co-volume area.
 * \return co-volume stats.
 */
[[nodiscard]] CoVolumeStats AnalyzeMeshCoVolumes(pmp::SurfaceMesh& mesh, const AreaFunction& areaFunction);

/*
 * \brief Computes scaling factor for stabilizing the finite volume method on assumed spherical surface meshes based on time step.
 * \param timeStep               time step size.
 * \param icoRadius              radius of an evolving geodesic icosahedron.
 * \param icoSubdiv              subdivision level of an evolving geodesic icosahedron.
 * \param stabilizationFactor    a multiplier for stabilizing mean co-volume area.
 * \return scaling factor for mesh and scalar grid.
 */
[[nodiscard]] float GetStabilizationScalingFactor(const double& timeStep, const float& icoRadius, const unsigned int& icoSubdiv, const float& stabilizationFactor = 1.0f);

/// \brief identifier for sparse matrix.
using SparseMatrix = Eigen::SparseMatrix<double>;

/// \brief A utility for converting Eigen::ComputationInfo to a string message.
[[nodiscard]] std::string InterpretSolverErrorCode(const Eigen::ComputationInfo& cInfo);

/**
 * \brief Computes angle-based tangential update velocity for a mesh vertex with a given weight.
 * \param mesh       a surface mesh to whom the vertex belongs.
 * \param v          vertex handle of a vertex where the tangential velocity is to be computed.
 * \param vNormal    normal to vertex with handle v.
 * \param weight     weight of the velocity vector.
 * \return tangential velocity vector.
 */
[[nodiscard]] pmp::vec3 ComputeTangentialUpdateVelocityAtVertex(const pmp::SurfaceMesh& mesh, const pmp::Vertex& v, const pmp::vec3& vNormal, const float& weight = 1.0f);
