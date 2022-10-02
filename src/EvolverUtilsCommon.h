#pragma once

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
