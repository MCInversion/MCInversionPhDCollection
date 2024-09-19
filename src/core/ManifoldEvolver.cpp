#include "ManifoldEvolver.h"

#include "pmp/algorithms/CurveFactory.h"
#include "pmp/algorithms/Features.h"
#include "pmp/algorithms/Normals.h"

#include "geometry/GridUtil.h"
#include "geometry/IcoSphereBuilder.h"
#include "geometry/MeshAnalysis.h"

#include "sdf/SDF.h"

#include "ConversionUtils.h"

/// \brief A factor by which the radius of any constructed outer/inner sphere is shrunken.
constexpr float SPHERE_RADIUS_FACTOR = 0.8f;

//
// ======================================================================================
//                    The strategy for 1D Curves in 2D space
// --------------------------------------------------------------------------------------
//

void ManifoldCurveEvolutionStrategy::Preprocess()
{
	const auto [minTargetSize, maxTargetSize, targetCenter] = ComputeAmbientFields();
	ConstructInitialManifolds(minTargetSize, maxTargetSize, targetCenter);

	GetFieldCellSize() = m_DistanceField ? m_DistanceField->CellSize() : minTargetSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	ComputeVariableDistanceFields();

	if (GetSettings().UseStabilizationViaScaling)
	{
		StabilizeGeometries();
	}
	AssignRemeshingSettingsToEvolvingManifolds();
	PrepareManifoldProperties();
}

void CustomManifoldCurveEvolutionStrategy::Preprocess()
{
	if (!GetOuterCurve() && GetInnerCurves().empty())
		throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: There's nothing to evolve!\n");
	
	if (!HasValidInnerOuterManifolds())
		throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: Invalid inner /outer manifold geometry! all custom inner curves are contained within the custom outer curve.\n");

	GetEvolBox() = GetOuterCurve()->bounds();
	const auto sizeVec = (GetEvolBox().max() - GetEvolBox().min()) * GetSettings().FieldSettings.FieldExpansionFactor;
	GetEvolBox().expand(sizeVec[0], sizeVec[1]);

	std::tie(std::ignore, std::ignore, std::ignore) = ComputeAmbientFields();

	const auto minTargetSize = std::min(sizeVec[0], sizeVec[1]);
	GetFieldCellSize() = GetDistanceField() ? GetDistanceField()->CellSize() : minTargetSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	ComputeVariableDistanceFields();

	if (GetSettings().UseStabilizationViaScaling)
	{
		const auto [minLength, maxLength] = CalculateCoVolumeRange();
		StabilizeCustomGeometries(minLength, maxLength);
	}
	AssignRemeshingSettingsToEvolvingManifolds();
	PrepareManifoldProperties();
}

void ManifoldCurveEvolutionStrategy::PerformEvolutionStep(unsigned int stepId)
{
	if (stepId > 1)
		ResetManifoldProperties();
	GetIntegrate()(stepId);
}

void ManifoldCurveEvolutionStrategy::Remesh()
{
	for (auto* curveToRemesh : m_RemeshTracker.GetManifoldsToRemesh())
	{
		pmp::CurveRemeshing remesher(*curveToRemesh);
		remesher.adaptive_remeshing(m_RemeshingSettings[curveToRemesh]);
	}
	m_RemeshTracker.Reset();
}

void ManifoldCurveEvolutionStrategy::ResizeRemeshingSettings(float resizeFactor)
{
	GetSettings().TimeStep *= pow(resizeFactor, 2); // we also need to adjust time step to keep the initial stabilization
	m_RemeshingSettings.AdjustAllRemeshingLengths(resizeFactor);
}

void ManifoldCurveEvolutionStrategy::DetectFeatures()
{
	// TODO: implement for pmp::ManifoldCurve2D
}

void ManifoldCurveEvolutionStrategy::ExportCurrentState(unsigned int step, const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_" + std::to_string(step);

	if (m_OuterCurve)
	{
		auto exportedOuterCurve = *m_OuterCurve;
		exportedOuterCurve *= m_TransformToOriginal;
		if (!write_to_ply(exportedOuterCurve, baseOutputFilename + "_Outer" + connectingName + ".ply"))
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Outer" + connectingName + ".ply") << "!\n";
	}

	for (size_t i = 0; const auto & innerCurve : m_InnerCurves)
	{
		auto exportedInnerCurve = *innerCurve;
		exportedInnerCurve *= m_TransformToOriginal;
		if (!write_to_ply(exportedInnerCurve, baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply"))
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply") << "!\n";
	}
}

void ManifoldCurveEvolutionStrategy::ExportFinalResult(const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_Result";

	if (m_OuterCurve)
	{
		auto exportedOuterCurve = *m_OuterCurve;
		exportedOuterCurve *= m_TransformToOriginal;
		if (!write_to_ply(exportedOuterCurve, baseOutputFilename + "_Outer" + connectingName + ".ply"))
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Outer" + connectingName + ".ply") << "!\n";
	}

	for (size_t i = 0; const auto & innerCurve : m_InnerCurves)
	{
		auto exportedInnerCurve = *innerCurve;
		exportedInnerCurve *= m_TransformToOriginal;
		if (!write_to_ply(exportedInnerCurve, baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply"))
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportCurrentState: error writing " << (baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply") << "!\n";
	}
}

void ManifoldCurveEvolutionStrategy::ExportTargetDistanceFieldAsImage(const std::string& baseOutputFilename)
{
	if (!m_DistanceField)
		return; // field not defined

	auto exportedField = *m_DistanceField;
	exportedField *= m_TransformToOriginal;
	exportedField /= static_cast<double>(GetScalingFactor());

	ExportScalarGrid2DToPNG(baseOutputFilename + "_TargetDF.png", exportedField, m_ScalarInterpolate, 
		10, 10, RAINBOW_TO_WHITE_MAP);
	ExportScalarGridDimInfo2D(baseOutputFilename + "_TargetDF.gdim2d", exportedField);
}

std::shared_ptr<pmp::ManifoldCurve2D> ManifoldCurveEvolutionStrategy::GetOuterCurveInOrigScale() const
{
	if (!m_OuterCurve)
		return nullptr;

	auto outerCurveTransformed = std::make_shared<pmp::ManifoldCurve2D>(*m_OuterCurve);
	*outerCurveTransformed *= m_TransformToOriginal;
	return outerCurveTransformed;
}

std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> ManifoldCurveEvolutionStrategy::GetInnerCurvesInOrigScale() const
{
	std::vector<std::shared_ptr<pmp::ManifoldCurve2D>> innerCurvesTransformed{};
	innerCurvesTransformed.reserve(m_InnerCurves.size());
	for (const auto& innerCurve : m_InnerCurves)
	{
		if (!innerCurve)
			continue;

		auto innerCurveTransformed = std::make_shared<pmp::ManifoldCurve2D>(*innerCurve);
		*innerCurveTransformed *= m_TransformToOriginal;
		innerCurvesTransformed.push_back(innerCurveTransformed);
	}
	innerCurvesTransformed.shrink_to_fit();
	return innerCurvesTransformed;
}

// -------------------------------------------------------------------------------------

void ManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep(unsigned int step)
{
	// ================================== Handle m_OuterCurve ==========================================================
	if (m_OuterCurve)
	{
		const auto NVertices = static_cast<unsigned int>(m_OuterCurve->n_vertices());
		SparseMatrix sysMat(NVertices, NVertices);
		Eigen::MatrixXd sysRhs(NVertices, 2);

		auto vDistanceToTarget = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:min_distance");
		const auto tStep = GetSettings().TimeStep;

		pmp::Normals2::compute_vertex_normals(*m_OuterCurve);
		auto vNormalsProp = m_OuterCurve->get_vertex_property<pmp::vec2>("v:normal");

		// prepare matrix & rhs for m_OuterCurve:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 2);

		for (const auto v : m_OuterCurve->vertices())
		{
			const auto vPosToUpdate = m_OuterCurve->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			for (const auto& innerCurveDf : m_InnerCurvesDistanceFields)
			{
				const auto innerDfAtVPos = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *innerCurveDf));
				if (innerDfAtVPos < vMinDistance[v])
					vMinDistance[v] = innerDfAtVPos;
			}

			if (m_OuterCurve->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector2d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);
			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = GetSettings().OuterManifoldEpsilon(static_cast<double>(vMinDistance[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const double advectionDistance = GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v];
			const double etaCtrlWeight = GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal);

			const Eigen::Vector2d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			if (tanRedistWeight > 0.0f)
			{
				// compute tangential velocity
				const auto vTanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_OuterCurve, v, vNormal, tanRedistWeight);
				sysRhs.row(v.idx()) += tStep * Eigen::Vector2d(vTanVelocity);
			}

			const auto laplaceWeightInfo = pmp::laplace_implicit_1D(*m_OuterCurve, v); // Laplacian weights
			tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0 + tStep * epsilonCtrlWeight * static_cast<double>(laplaceWeightInfo.weightSum)));

			for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
			{
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), w.idx(), -1.0 * tStep * epsilonCtrlWeight * static_cast<double>(weight)));
			}
		}

		// After the loop
		sysMat.setFromTriplets(tripletList.begin(), tripletList.end());
		// if (IsRemeshingNecessary(*m_OuterCurve, m_RemeshingSettings[m_OuterCurve.get()]))
			m_RemeshTracker.AddManifold(m_OuterCurve.get());

		// solve
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
		Eigen::MatrixXd x = solver.solve(sysRhs);
		if (solver.info() != Eigen::Success)
		{
			const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(step) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}

		// update vertex positions & verify mesh within bounds
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);
			if (!m_EvolBox.Contains(newPos))
			{
				const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: vertex " + std::to_string(i) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}
			m_OuterCurve->position(pmp::Vertex(i)) = newPos;
		}
	}

	// ================================== Handle m_InnerCurves ==========================================================
	for (const auto& innerCurve : m_InnerCurves)
	{
		const auto NVertices = static_cast<unsigned int>(innerCurve->n_vertices());
		SparseMatrix sysMat(NVertices, NVertices);
		Eigen::MatrixXd sysRhs(NVertices, 2);

		auto vDistanceToTarget = innerCurve->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = innerCurve->get_vertex_property<pmp::Scalar>("v:min_distance");
		const auto tStep = GetSettings().TimeStep;

		pmp::Normals2::compute_vertex_normals(*innerCurve);
		auto vNormalsProp = innerCurve->get_vertex_property<pmp::vec2>("v:normal");

		// prepare matrix & rhs for m_OuterCurve:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 2);  // Assuming 2 entries per vertex for curves

		for (const auto v : innerCurve->vertices())
		{
			const auto vPosToUpdate = innerCurve->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			if (m_OuterCurveDistanceField)
			{
				const auto outerDfAtVPos = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_OuterCurveDistanceField));
				if (outerDfAtVPos < vMinDistance[v])
					vMinDistance[v] = outerDfAtVPos;
			}

			if (innerCurve->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector2d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);
			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = -1.0 * GetSettings().InnerManifoldEpsilon(static_cast<double>(vMinDistance[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const double advectionDistance = GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v];
			const double etaCtrlWeight = GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal);

			const Eigen::Vector2d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			if (tanRedistWeight > 0.0f)
			{
				// compute tangential velocity
				const auto vTanVelocity = ComputeTangentialUpdateVelocityAtVertex(*innerCurve, v, vNormal, tanRedistWeight);
				sysRhs.row(v.idx()) += tStep * Eigen::Vector2d(vTanVelocity);
			}

			const auto laplaceWeightInfo = pmp::laplace_implicit_1D(*innerCurve, v); // Laplacian weights
			tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0 + tStep * epsilonCtrlWeight * static_cast<double>(laplaceWeightInfo.weightSum)));

			for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
			{
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), w.idx(), -1.0 * tStep * epsilonCtrlWeight * static_cast<double>(weight)));
			}
		}

		// After the loop
		sysMat.setFromTriplets(tripletList.begin(), tripletList.end());
		if (IsRemeshingNecessary(*innerCurve, m_RemeshingSettings[innerCurve.get()]))
			m_RemeshTracker.AddManifold(innerCurve.get());

		// solve
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
		Eigen::MatrixXd x = solver.solve(sysRhs);
		if (solver.info() != Eigen::Success)
		{
			const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(step) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}

		// update vertex positions & verify mesh within bounds
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);
			if (!m_EvolBox.Contains(newPos))
			{
				const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: innerCurve vertex " + std::to_string(i) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}
			innerCurve->position(pmp::Vertex(i)) = newPos;
		}
	}
}

void ManifoldCurveEvolutionStrategy::ExplicitIntegrationStep(unsigned int step)
{
	// ================================== Handle m_OuterCurve ==========================================================
	if (m_OuterCurve)
	{
		const auto tStep = GetSettings().TimeStep;
		auto vDistanceToTarget = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vMinDistanceToInner = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:min_distance_to_inner");

		pmp::Normals2::compute_vertex_normals(*m_OuterCurve);
		auto vNormalsProp = m_OuterCurve->get_vertex_property<pmp::vec2>("v:normal");

		for (const auto v : m_OuterCurve->vertices())
		{
			const auto vPosToUpdate = m_OuterCurve->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			for (const auto& innerSurfaceDf : m_InnerCurvesDistanceFields)
			{
				const auto innerDfAtVPos = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *innerSurfaceDf));
				if (innerDfAtVPos < vMinDistance[v])
					vMinDistance[v] = innerDfAtVPos;
				if (innerDfAtVPos < vMinDistanceToInner[v])
					vMinDistanceToInner[v] = innerDfAtVPos;
			}

			if (m_OuterCurve->is_boundary(v))
				continue; // skip boundary vertices

			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);
			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = 
				GetSettings().OuterManifoldEpsilon(static_cast<double>(vMinDistance[v])) +
				GetSettings().OuterManifoldRepulsion(static_cast<double>(vMinDistanceToInner[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const double advectionDistance = GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v];
			const double etaCtrlWeight = GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal);

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * pmp::laplace_1D(*m_OuterCurve, v);

			// Tangential redistribution velocity
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			pmp::vec2 tanVelocity(0.0, 0.0);
			if (tanRedistWeight > 0.0f)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_OuterCurve, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldCurveEvolutionStrategy::ExplicitIntegrationStep: vertex " + std::to_string(v.idx()) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			m_OuterCurve->position(v) = updatedPosition;
		}
	}

	// ================================== Handle m_InnerCurves ==========================================================
	for (const auto& innerCurve : m_InnerCurves)
	{
		const auto tStep = GetSettings().TimeStep;
		auto vDistanceToTarget = innerCurve->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = innerCurve->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vDistanceToOuter = innerCurve->get_vertex_property<pmp::Scalar>("v:distance_to_outer");

		pmp::Normals2::compute_vertex_normals(*innerCurve);
		auto vNormalsProp = innerCurve->get_vertex_property<pmp::vec2>("v:normal");

		for (const auto v : innerCurve->vertices())
		{
			const auto vPosToUpdate = innerCurve->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			if (m_OuterCurveDistanceField)
			{
				vDistanceToOuter[v] = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_OuterCurveDistanceField));
				if (vDistanceToOuter[v] < vMinDistance[v])
					vMinDistance[v] = vDistanceToOuter[v];
			}

			if (innerCurve->is_boundary(v))
				continue; // skip boundary vertices
			
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);
			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = 
				-1.0 * GetSettings().InnerManifoldEpsilon(static_cast<double>(vMinDistance[v])) +
				GetSettings().InnerManifoldRepulsion(static_cast<double>(vDistanceToOuter[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const double advectionDistance = GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v];
			const double etaCtrlWeight = GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal);

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * pmp::laplace_1D(*innerCurve, v);

			// Tangential redistribution velocity
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			pmp::vec2 tanVelocity(0.0, 0.0);
			if (tanRedistWeight > 0.0f)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*innerCurve, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldCurveEvolutionStrategy::ExplicitIntegrationStep: innerCurve vertex " + std::to_string(v.idx()) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			innerCurve->position(v) = updatedPosition;
		}
	}
}

void ManifoldCurveEvolutionStrategy::PrepareManifoldProperties()
{
	if (m_OuterCurve)
	{
		// distance to m_TargetPointCloud
		m_OuterCurve->add_vertex_property<pmp::Scalar>("v:distance_to_target", FLT_MAX);

		// the minimum from the distances to m_TargetPointCloud and all m_InnerCurves
		m_OuterCurve->add_vertex_property<pmp::Scalar>("v:min_distance", FLT_MAX);

		// the minimum from the distances to m_InnerCurves.
		m_OuterCurve->add_vertex_property<pmp::Scalar>("v:min_distance_to_inner", FLT_MAX);
	}

	for (const auto& innerCurve : m_InnerCurves)
	{
		// distance to m_TargetPointCloud
		innerCurve->add_vertex_property<pmp::Scalar>("v:distance_to_target", FLT_MAX);

		// the minimum from the distances to m_TargetPointCloud and m_OuterCurve
		innerCurve->add_vertex_property<pmp::Scalar>("v:min_distance", FLT_MAX);

		// distance to m_OuterCurve
		innerCurve->add_vertex_property<pmp::Scalar>("v:distance_to_outer", FLT_MAX);
	}

	// v:normal will be added during normal computation: pmp::Normals2::compute_vertex_normals
}

void ManifoldCurveEvolutionStrategy::ResetManifoldProperties()
{
	// TODO: find a better way to reset prop values globally, perhaps implement a new property method
	if (m_OuterCurve)
	{
		auto vDistanceToTarget = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vMinDistanceToInner = m_OuterCurve->get_vertex_property<pmp::Scalar>("v:min_distance_to_inner");

		for (const auto v : m_OuterCurve->vertices())
		{
			vDistanceToTarget[v] = FLT_MAX;
			vMinDistance[v] = FLT_MAX;
			vMinDistanceToInner[v] = FLT_MAX;
		}
	}

	for (const auto& innerCurve : m_InnerCurves)
	{
		auto vDistanceToTarget = innerCurve->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = innerCurve->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vDistanceToOuter = innerCurve->get_vertex_property<pmp::Scalar>("v:distance_to_outer");

		for (const auto v : innerCurve->vertices())
		{
			vDistanceToTarget[v] = FLT_MAX;
			vMinDistance[v] = FLT_MAX;
			vDistanceToOuter[v] = FLT_MAX;
		}
	}
}

std::tuple<float, float, pmp::Point2> ManifoldCurveEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldCurveEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		return { FLT_MAX, FLT_MAX, pmp::Point2(0, 0)};
	}

	const pmp::BoundingBox2 ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const float maxSize = std::max(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const float cellSize = minSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	const SDF::PointCloudDistanceField2DSettings dfSettings{
		cellSize,
		GetSettings().FieldSettings.FieldExpansionFactor,
		DBL_MAX
	};
	m_DistanceField = std::make_shared<Geometry::ScalarGrid2D>(
		SDF::PlanarPointCloudDistanceFieldGenerator::Generate(*m_TargetPointCloud, dfSettings));
	m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid2D>(ComputeNormalizedNegativeGradient(*m_DistanceField));
	return { minSize, maxSize, ptCloudBBox.center() };
}

void ManifoldCurveEvolutionStrategy::ComputeVariableDistanceFields()
{
	if (!m_OuterCurve || m_InnerCurves.empty())
	{
		// there's no possibility of interaction between the outer and the inner manifolds
		return;
	}

	const SDF::DistanceField2DSettings curveDFSettings{
		GetFieldCellSize(),
		GetSettings().FieldSettings.FieldExpansionFactor,
		DBL_MAX,
		SDF::KDTreeSplitType::Center,
		SDF::SignComputation2D::None,
		SDF::PreprocessingType2D::Quadtree
	};

	const Geometry::ManifoldCurve2DAdapter outerCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(*m_OuterCurve));
	m_OuterCurveDistanceField = std::make_shared<Geometry::ScalarGrid2D>(
		SDF::PlanarDistanceFieldGenerator::Generate(outerCurveAdapter, curveDFSettings));

	for (const auto& innerCurve : m_InnerCurves)
	{
		const Geometry::ManifoldCurve2DAdapter innerCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(*innerCurve));
		m_InnerCurvesDistanceFields.emplace_back(std::make_shared<Geometry::ScalarGrid2D>(
			SDF::PlanarDistanceFieldGenerator::Generate(innerCurveAdapter, curveDFSettings)));
	}
}

/// \brief the smallest allowed number of vertices in a manifold curve.
constexpr unsigned int N_CIRCLE_VERTS_0{ 5 };

void ManifoldCurveEvolutionStrategy::ConstructInitialManifolds(float minTargetSize, float maxTargetSize, const pmp::Point2& targetBoundsCenter)
{
	if (!GetSettings().UseInnerManifolds && !GetSettings().UseOuterManifolds)
	{
		throw std::invalid_argument("ManifoldCurveEvolutionStrategy::ConstructInitialManifolds: Current setting is: UseInnerManifolds == false && UseOuterManifolds == false. This means there's nothing to evolve!\n");
	}

	const float outerCircleRadius = 0.5f * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5f + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);
	const auto nSegments = static_cast<unsigned int>(pow(2, GetSettings().LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;

	if (GetSettings().UseOuterManifolds)
	{
		m_OuterCurve = std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(pmp::Point2{}, outerCircleRadius, nSegments));
		m_InitialSphereSettings[m_OuterCurve.get()] = Circle2D{ targetBoundsCenter, outerCircleRadius };
	}

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return;

	//const InscribedCircleInputData calcData{
	//	*m_TargetPointCloud,
	//	std::make_shared<Geometry::ScalarGrid2D>(*m_DistanceField) // clone
	//};
	//ParticleSwarmDistanceFieldInscribedCircleCalculator inscribedCircleCalculator;
	//const auto circles = inscribedCircleCalculator.Calculate(calcData);

	// Hardcoded inner curves:
	const auto circles = std::vector{ Circle2D{targetBoundsCenter, 1.85f} };

	m_InnerCurves.reserve(circles.size());

	for (const auto& circle : circles)
	{
		// keep the same vertex density for inner circles
		const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * (circle.Radius) / outerCircleRadius);
		m_InnerCurves.emplace_back(std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(
			circle.Center,
			circle.Radius * SPHERE_RADIUS_FACTOR, 
			nInnerSegments
		)));
		m_InitialSphereSettings[m_InnerCurves.back().get()] = circle;
	}
}

/// \brief The power of the stabilizing scale factor.
constexpr float SCALE_FACTOR_POWER_1D = 1.0f;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr float INV_SHRINK_FACTOR_1D = 20.0f;

void ManifoldCurveEvolutionStrategy::StabilizeGeometries(float stabilizationFactor)
{
	const auto radius = stabilizationFactor * m_InitialSphereSettings.MinRadius() + (1.0f - stabilizationFactor) * m_InitialSphereSettings.MaxRadius();
	if (radius <= 0.0f)
	{
		throw std::invalid_argument("ManifoldCurveEvolutionStrategy::StabilizeGeometries: m_InitialSphereSettings empty!\n");
	}
	const auto expectedVertexCount = static_cast<unsigned int>(pow(2, GetSettings().LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
	const auto expectedMeanCoVolLength = (2.0f * static_cast<float>(M_PI) * radius / static_cast<float>(expectedVertexCount));
	const auto scalingFactor = pow(static_cast<float>(GetSettings().TimeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_1D, SCALE_FACTOR_POWER_1D);
	GetScalingFactor() = scalingFactor;
	GetSettings().FieldSettings.FieldIsoLevel *= scalingFactor;

	const pmp::mat3 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 1.0f
	};
	const pmp::Point2 origin = m_OuterCurve ? m_InitialSphereSettings[m_OuterCurve.get()].Center :
		(!m_InnerCurves.empty() ? m_InitialSphereSettings[m_InnerCurves[0].get()].Center : pmp::Point2{ 0, 0 });
	const pmp::mat3 transfMatrixGeomMove{
		1.0f, 0.0f, -origin[0],
		0.0f, 1.0f, -origin[1],
		0.0f, 0.0f, 1.0f
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	m_TransformToOriginal = inverse(transfMatrixFull);

	// transform geometries
	if (m_OuterCurve)
	{
		(*m_OuterCurve) *= transfMatrixGeomScale;
	}
	for (const auto& innerCurve : m_InnerCurves)
	{
		(*innerCurve) *= transfMatrixGeomScale;
	}

	// test box for geometry validation
	const float evolBoxFactor = 1.2f * scalingFactor;
	m_EvolBox = pmp::BoundingBox2(
		pmp::Point2{ -radius, -radius } * evolBoxFactor,
		pmp::Point2{ radius, radius } * evolBoxFactor);

	if (m_OuterCurveDistanceField)
	{
		(*m_OuterCurveDistanceField) *= transfMatrixFull;
		(*m_OuterCurveDistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerCurveDF : m_InnerCurvesDistanceFields)
	{
		(*innerCurveDF) *= transfMatrixFull;
		(*innerCurveDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}

	if (!m_DistanceField || !m_DFNegNormalizedGradient)
		return; // nothing to scale

	(*m_DistanceField) *= transfMatrixFull;
	(*m_DFNegNormalizedGradient) *= transfMatrixFull;
	(*m_DistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
}

void ManifoldCurveEvolutionStrategy::AssignRemeshingSettingsToEvolvingManifolds()
{
	if (m_OuterCurve)
	{
		const Circle2D& circleSettings = m_InitialSphereSettings[m_OuterCurve.get()];
		//m_RemeshingSettings[m_OuterCurve.get()] = CollectRemeshingSettingsFromCircleCurve(m_OuterCurve, 
		//	circleSettings.Radius * GetScalingFactor(), 
		//	circleSettings.Center);
		m_RemeshingSettings[m_OuterCurve.get()] = CollectRemeshingSettingsFromIcoSphere_OLD(
			GetSettings().LevelOfDetail,
			circleSettings.Radius * GetScalingFactor(),
			GetSettings().RemeshingSettings);
	}

	for (const auto& innerCurve : m_InnerCurves)
	{
		const Circle2D& circleSettings = m_InitialSphereSettings[innerCurve.get()];
		//m_RemeshingSettings[innerCurve.get()] = CollectRemeshingSettingsFromCircleCurve(innerCurve, 
		//	circleSettings.Radius * GetScalingFactor(),
		//	circleSettings.Center);
		m_RemeshingSettings[innerCurve.get()] = CollectRemeshingSettingsFromIcoSphere_OLD(
			GetSettings().LevelOfDetail,
			circleSettings.Radius * GetScalingFactor(),
			GetSettings().RemeshingSettings);
	}
}

void CustomManifoldCurveEvolutionStrategy::AssignRemeshingSettingsToEvolvingManifolds()
{
	if (GetOuterCurve())
	{
		GetRemeshingSettings()[GetOuterCurve().get()] = CollectRemeshingSettingsFromCurve(GetOuterCurve());
	}

	for (const auto& innerCurve : GetInnerCurves())
	{
		GetRemeshingSettings()[innerCurve.get()] = CollectRemeshingSettingsFromCurve(innerCurve);
	}
}

bool CustomManifoldCurveEvolutionStrategy::HasValidInnerOuterManifolds() const
{
	// check self-intersections
	if (GetOuterCurve() && Geometry::PMPManifoldCurve2DHasSelfIntersections(*GetOuterCurve()))
		return false;

	const auto& outerCurve = *GetOuterCurve();
	// check inner curves
	for (const auto& innerCurve : GetInnerCurves())
	{
		// check self-intersections
		if (Geometry::PMPManifoldCurve2DHasSelfIntersections(*innerCurve))
			return false;

		// check whether the inner curve is contained inside the outer curve.
		for (const auto& p : innerCurve->positions())
		{
			if (!Geometry::IsPointInsidePMPManifoldCurve(p, outerCurve))
				return false;
		}
	}

	return true;
}

std::pair<float, float> CustomManifoldCurveEvolutionStrategy::CalculateCoVolumeRange() const
{
	float minCoVolLength = FLT_MAX;
	float maxCoVolLength = -FLT_MAX;

	if (GetOuterCurve())
	{
		const auto& outerCurve = *GetOuterCurve();
		for (const auto v : outerCurve.vertices())
		{
			const auto [eTo, eFrom] = outerCurve.edges(v);
			const auto currentCoVolLength = 0.5f * (outerCurve.edge_length(eTo) + outerCurve.edge_length(eFrom));

			if (currentCoVolLength > maxCoVolLength) maxCoVolLength = currentCoVolLength;
			if (currentCoVolLength < minCoVolLength) minCoVolLength = currentCoVolLength;
		}
	}

	for (const auto& curve : GetInnerCurves())
	{
		const auto& innerCurve = *curve;
		for (const auto v : innerCurve.vertices())
		{
			const auto [eTo, eFrom] = innerCurve.edges(v);
			const auto currentCoVolLength = 0.5f * (innerCurve.edge_length(eTo) + innerCurve.edge_length(eFrom));

			if (currentCoVolLength > maxCoVolLength) maxCoVolLength = currentCoVolLength;
			if (currentCoVolLength < minCoVolLength) minCoVolLength = currentCoVolLength;
		}
	}

	return { minCoVolLength, maxCoVolLength };
}

void CustomManifoldCurveEvolutionStrategy::StabilizeCustomGeometries(float minLength, float maxLength, float stabilizationFactor)
{
	const float expectedMeanCoVolLength = (1.0f - stabilizationFactor) * minLength + stabilizationFactor * maxLength;
	const float scalingFactor = pow(static_cast<float>(GetSettings().TimeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_1D, SCALE_FACTOR_POWER_1D);
	GetScalingFactor() = scalingFactor;
	const pmp::mat3 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 1.0f
	};
	GetTransformToOriginal() = inverse(transfMatrixGeomScale);

	// Transform the geometries
	(*GetOuterCurve()) *= transfMatrixGeomScale;
	for (auto& innerCurve : GetInnerCurves())
	{
		(*innerCurve) *= transfMatrixGeomScale;
	}

	if (GetOuterCurveDistanceField())
	{
		(*GetOuterCurveDistanceField()) *= transfMatrixGeomScale;
		(*GetOuterCurveDistanceField()) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerCurveDF : GetInnerCurvesDistanceFields())
	{
		(*innerCurveDF) *= transfMatrixGeomScale;
		(*innerCurveDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}

	if (!GetDistanceField() || !GetDFNegNormalizedGradient())
		return; // nothing to scale
	(*GetDistanceField()) *= transfMatrixGeomScale;
	(*GetDFNegNormalizedGradient()) *= transfMatrixGeomScale;
	(*GetDistanceField()) *= static_cast<double>(scalingFactor); // Scale the distance values
}

//
// ======================================================================================
//                    The strategy for 2D Surfaces in 3D space
// ---------------------------------------------------------------------------------------
//

void ManifoldSurfaceEvolutionStrategy::Preprocess()
{
	const auto [minTargetSize, maxTargetSize, targetCenter] = ComputeAmbientFields();
	ConstructInitialManifolds(minTargetSize, maxTargetSize, targetCenter);

	GetFieldCellSize() = m_DistanceField ? m_DistanceField->CellSize() : minTargetSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	ComputeVariableDistanceFields();

	if (GetSettings().UseStabilizationViaScaling)
	{
		StabilizeGeometries();
	}
	AssignRemeshingSettingsToEvolvingManifolds();
	PrepareManifoldProperties();
}

void CustomManifoldSurfaceEvolutionStrategy::Preprocess()
{
	if (!GetOuterSurface() && GetInnerSurfaces().empty())
		throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: There's nothing to evolve!\n");

	if (!HasValidInnerOuterManifolds())
		throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: Invalid inner /outer manifold geometry! all custom inner curves are contained within the custom outer curve.\n");

	GetEvolBox() = GetOuterSurface()->bounds();
	const auto sizeVec = (GetEvolBox().max() - GetEvolBox().min()) * GetSettings().FieldSettings.FieldExpansionFactor;
	GetEvolBox().expand(sizeVec[0], sizeVec[1], sizeVec[2]);

	std::tie(std::ignore, std::ignore, std::ignore) = ComputeAmbientFields();

	const auto minTargetSize = std::min({ sizeVec[0], sizeVec[1], sizeVec[2] });
	GetFieldCellSize() = GetDistanceField() ? GetDistanceField()->CellSize() : minTargetSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	ComputeVariableDistanceFields();

	if (GetSettings().UseStabilizationViaScaling)
	{
		const auto [minArea, maxArea] = CalculateCoVolumeRange();
		StabilizeCustomGeometries(minArea, maxArea);
	}
	AssignRemeshingSettingsToEvolvingManifolds();
	PrepareManifoldProperties();
}

void ManifoldSurfaceEvolutionStrategy::PerformEvolutionStep(unsigned int stepId)
{
	if (stepId > 1)
		ResetManifoldProperties();
	GetIntegrate()(stepId);
}

void ManifoldSurfaceEvolutionStrategy::Remesh()
{
	for (auto* surfaceToRemesh : m_RemeshTracker.GetManifoldsToRemesh())
	{
		pmp::Remeshing remesher(*surfaceToRemesh);
		remesher.adaptive_remeshing(m_RemeshingSettings[surfaceToRemesh]);
	}
	m_RemeshTracker.Reset();
}

void ManifoldSurfaceEvolutionStrategy::ResizeRemeshingSettings(float resizeFactor)
{
	GetSettings().TimeStep *= pow(resizeFactor, 2); // we also need to adjust time step to keep the initial stabilization
	m_RemeshingSettings.AdjustAllRemeshingLengths(resizeFactor);
}

void ManifoldSurfaceEvolutionStrategy::DetectFeatures()
{
	for (auto* surfaceToRemesh : m_RemeshTracker.GetManifoldsToRemesh())
	{
		pmp::Features featuresDetector(*surfaceToRemesh);
		const auto curvatureAngle = GetSettings().FeatureSettings.CriticalMeanCurvatureAngle;
		const auto curvatureFactor = GetSettings().FeatureSettings.PrincipalCurvatureFactor;
		const auto nFeature = featuresDetector.detect_vertices_with_high_curvature(curvatureAngle, curvatureFactor, true);
	}
}

void ManifoldSurfaceEvolutionStrategy::ExportCurrentState(unsigned int step, const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_" + std::to_string(step);

	if (m_OuterSurface)
	{
		auto exportedOuterSurface = *m_OuterSurface;
		exportedOuterSurface *= m_TransformToOriginal;
		exportedOuterSurface.write(baseOutputFilename + "_Outer" + connectingName + ".vtk");		
	}

	for (size_t i = 0; const auto& innerSurface : m_InnerSurfaces)
	{
		auto exportedInnerSurface = *innerSurface;
		exportedInnerSurface *= m_TransformToOriginal;
		exportedInnerSurface.write(baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".vtk");
	}
}

void ManifoldSurfaceEvolutionStrategy::ExportFinalResult(const std::string& baseOutputFilename)
{
	const std::string connectingName = "_Evol_Result";

	if (m_OuterSurface)
	{
		auto exportedOuterSurface = *m_OuterSurface;
		exportedOuterSurface *= m_TransformToOriginal;
		exportedOuterSurface.write(baseOutputFilename + "_Outer" + connectingName + ".vtk");
	}

	for (size_t i = 0; const auto & innerSurface : m_InnerSurfaces)
	{
		auto exportedInnerSurface = *innerSurface;
		exportedInnerSurface *= m_TransformToOriginal;
		exportedInnerSurface.write(baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".vtk");
	}
}

void ManifoldSurfaceEvolutionStrategy::ExportTargetDistanceFieldAsImage(const std::string& baseOutputFilename)
{
	if (!m_DistanceField)
		return; // field not defined

	auto exportedField = *m_DistanceField;
	exportedField *= m_TransformToOriginal;
	exportedField /= static_cast<double>(GetScalingFactor());

	ExportToVTI(baseOutputFilename + "_TargetDF", exportedField);
}

std::shared_ptr<pmp::SurfaceMesh> ManifoldSurfaceEvolutionStrategy::GetOuterSurfaceInOrigScale() const
{
	if (!m_OuterSurface)
		return nullptr;

	auto outerSurfaceTransformed = std::make_shared<pmp::SurfaceMesh>(*m_OuterSurface);
	*outerSurfaceTransformed *= m_TransformToOriginal;
	return outerSurfaceTransformed;
}

std::vector<std::shared_ptr<pmp::SurfaceMesh>> ManifoldSurfaceEvolutionStrategy::GetInnerSurfacesInOrigScale() const
{
	std::vector<std::shared_ptr<pmp::SurfaceMesh>> innerSurfacesTransformed{};
	innerSurfacesTransformed.reserve(m_InnerSurfaces.size());
	for (const auto& innerSurface : m_InnerSurfaces)
	{
		if (!innerSurface)
			continue;

		auto innerCurveTransformed = std::make_shared<pmp::SurfaceMesh>(*innerSurface);
		*innerCurveTransformed *= m_TransformToOriginal;
		innerSurfacesTransformed.push_back(innerCurveTransformed);
	}
	innerSurfacesTransformed.shrink_to_fit();
	return innerSurfacesTransformed;
}

// ------------------------------------------------


void ManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep(unsigned int step)
{
	// ================================== Handle m_OuterSurface ==========================================================
	if (m_OuterSurface)
	{
		const auto NVertices = static_cast<unsigned int>(m_OuterSurface->n_vertices());
		SparseMatrix sysMat(NVertices, NVertices);
		Eigen::MatrixXd sysRhs(NVertices, 3);

		auto vDistanceToTarget = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vMinDistanceToInner = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:min_distance_to_inner");
		const auto tStep = GetSettings().TimeStep;

		pmp::Normals::compute_vertex_normals(*m_OuterSurface);
		auto vNormalsProp = m_OuterSurface->get_vertex_property<pmp::vec3>("v:normal");

		// prepare matrix & rhs for m_OuterSurface:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 6);  // Assuming an average of 6 entries per vertex

		for (const auto v : m_OuterSurface->vertices())
		{
			const auto vPosToUpdate = m_OuterSurface->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			for (const auto& innerSurfaceDf : m_InnerSurfacesDistanceFields)
			{
				const auto innerDfAtVPos = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *innerSurfaceDf));
				if (innerDfAtVPos < vMinDistance[v])
					vMinDistance[v] = innerDfAtVPos;
				if (innerDfAtVPos < vMinDistanceToInner[v])
					vMinDistanceToInner[v] = innerDfAtVPos;
			}

			if (m_OuterSurface->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector3d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);
			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = 
				GetSettings().OuterManifoldEpsilon(static_cast<double>(vMinDistance[v])) +
				GetSettings().OuterManifoldRepulsion(static_cast<double>(vMinDistanceToInner[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const double advectionDistance = GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v];
			const double etaCtrlWeight = GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal);

			const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			if (tanRedistWeight > 0.0f)
			{
				// compute tangential velocity
				const auto vTanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_OuterSurface, v, vNormal, tanRedistWeight);
				sysRhs.row(v.idx()) += tStep * Eigen::Vector3d(vTanVelocity);
			}

			const auto laplaceWeightInfo = m_ImplicitLaplacianFunction(*m_OuterSurface, v); // Laplacian weights
			tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0 + tStep * epsilonCtrlWeight * static_cast<double>(laplaceWeightInfo.weightSum)));

			for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
			{
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), w.idx(), -1.0 * tStep * epsilonCtrlWeight * static_cast<double>(weight)));
			}
		}

		// After the loop
		sysMat.setFromTriplets(tripletList.begin(), tripletList.end());
		//if (IsRemeshingNecessary(*m_OuterSurface, m_RemeshingSettings[m_OuterSurface.get()], m_LaplacianAreaFunction))
		if (IsRemeshingNecessary(*m_OuterSurface, GetSettings().QualitySettings.FaceQualityFunc, GetSettings().QualitySettings.Range))
			m_RemeshTracker.AddManifold(m_OuterSurface.get());

		// solve
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
		Eigen::MatrixXd x = solver.solve(sysRhs);
		if (solver.info() != Eigen::Success)
		{
			const std::string msg = "\nManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(step) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}

		// update vertex positions & verify mesh within bounds
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);
			if (!m_EvolBox.Contains(newPos))
			{
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep: vertex " + std::to_string(i) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}
			m_OuterSurface->position(pmp::Vertex(i)) = newPos;
		}
	}

	// ================================== Handle m_InnerSurfaces ==========================================================
	for (const auto& innerSurface : m_InnerSurfaces)
	{
		const auto NVertices = static_cast<unsigned int>(innerSurface->n_vertices());
		SparseMatrix sysMat(NVertices, NVertices);
		Eigen::MatrixXd sysRhs(NVertices, 3);

		auto vDistanceToTarget = innerSurface->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = innerSurface->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vDistanceToOuter = innerSurface->get_vertex_property<pmp::Scalar>("v:distance_to_outer");
		const auto tStep = GetSettings().TimeStep;

		pmp::Normals::compute_vertex_normals(*innerSurface);
		auto vNormalsProp = innerSurface->get_vertex_property<pmp::vec3>("v:normal");

		// prepare matrix & rhs for m_InnerSurface:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 6);  // Assuming an average of 6 entries per vertex

		for (const auto v : innerSurface->vertices())
		{
			const auto vPosToUpdate = innerSurface->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			if (m_OuterSurfaceDistanceField)
			{
				vDistanceToOuter[v] = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_OuterSurfaceDistanceField));
				if (vDistanceToOuter[v] < vMinDistance[v])
					vMinDistance[v] = vDistanceToOuter[v];
			}

			if (innerSurface->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector3d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);
			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = 
				-1.0 * GetSettings().InnerManifoldEpsilon(static_cast<double>(vMinDistance[v])) +
				GetSettings().InnerManifoldRepulsion(static_cast<double>(vDistanceToOuter[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const double advectionDistance = GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v];
			const double etaCtrlWeight = GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal);

			const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			if (tanRedistWeight > 0.0f)
			{
				// compute tangential velocity
				const auto vTanVelocity = ComputeTangentialUpdateVelocityAtVertex(*innerSurface, v, vNormal, tanRedistWeight);
				sysRhs.row(v.idx()) += tStep * Eigen::Vector3d(vTanVelocity);
			}

			const auto laplaceWeightInfo = m_ImplicitLaplacianFunction(*innerSurface, v); // Laplacian weights
			tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0 + tStep * epsilonCtrlWeight * static_cast<double>(laplaceWeightInfo.weightSum)));

			for (const auto& [w, weight] : laplaceWeightInfo.vertexWeights)
			{
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), w.idx(), -1.0 * tStep * epsilonCtrlWeight * static_cast<double>(weight)));
			}
		}

		// After the loop
		sysMat.setFromTriplets(tripletList.begin(), tripletList.end());
		//if (IsRemeshingNecessary(*innerSurface, m_RemeshingSettings[innerSurface.get()], m_LaplacianAreaFunction))
		if (IsRemeshingNecessary(*innerSurface, GetSettings().QualitySettings.FaceQualityFunc, GetSettings().QualitySettings.Range))
			m_RemeshTracker.AddManifold(innerSurface.get());

		// solve
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
		Eigen::MatrixXd x = solver.solve(sysRhs);
		if (solver.info() != Eigen::Success)
		{
			const std::string msg = "\nManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(step) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}

		// update vertex positions & verify mesh within bounds
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);
			if (!m_EvolBox.Contains(newPos))
			{
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep: innerSurface vertex " + std::to_string(i) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}
			innerSurface->position(pmp::Vertex(i)) = newPos;
		}
	}
}

void ManifoldSurfaceEvolutionStrategy::ExplicitIntegrationStep(unsigned int step)
{
	// ================================== Handle m_OuterSurface ==========================================================
	if (m_OuterSurface)
	{
		const auto tStep = GetSettings().TimeStep;
		auto vDistanceToTarget = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:min_distance");

		pmp::Normals::compute_vertex_normals(*m_OuterSurface);
		auto vNormalsProp = m_OuterSurface->get_vertex_property<pmp::vec3>("v:normal");

		for (const auto v : m_OuterSurface->vertices())
		{
			const auto vPosToUpdate = m_OuterSurface->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			for (const auto& innerSurfaceDf : m_InnerSurfacesDistanceFields)
			{
				const auto innerDfAtVPos = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *innerSurfaceDf));
				if (innerDfAtVPos < vMinDistance[v])
					vMinDistance[v] = innerDfAtVPos;
			}

			if (m_OuterSurface->is_boundary(v))
				continue; // skip boundary vertices

			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);
			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = GetSettings().OuterManifoldEpsilon(static_cast<double>(vMinDistance[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const double advectionDistance = GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v];
			const double etaCtrlWeight = GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal);

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * m_ExplicitLaplacianFunction(*m_OuterSurface, v);

			// Tangential redistribution velocity
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			pmp::vec3 tanVelocity(0.0, 0.0, 0.0);
			if (tanRedistWeight > 0.0f)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_OuterSurface, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::ExplicitIntegrationStep: vertex " + std::to_string(v.idx()) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			m_OuterSurface->position(v) = updatedPosition;
		}
	}

	// ================================== Handle m_InnerSurfaces ==========================================================
	for (const auto& innerSurface : m_InnerSurfaces)
	{
		const auto tStep = GetSettings().TimeStep;
		auto vDistanceToTarget = innerSurface->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = innerSurface->get_vertex_property<pmp::Scalar>("v:min_distance");

		pmp::Normals::compute_vertex_normals(*innerSurface);
		auto vNormalsProp = innerSurface->get_vertex_property<pmp::vec3>("v:normal");

		for (const auto v : innerSurface->vertices())
		{
			const auto vPosToUpdate = innerSurface->position(v);
			vDistanceToTarget[v] = m_DistanceField ? static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_DistanceField)) : FLT_MAX;
			vDistanceToTarget[v] -= GetSettings().FieldSettings.FieldIsoLevel;
			vMinDistance[v] = vDistanceToTarget[v];
			if (m_OuterSurfaceDistanceField)
			{
				const auto outerDfAtVPos = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_OuterSurfaceDistanceField));
				if (outerDfAtVPos < vMinDistance[v])
					vMinDistance[v] = outerDfAtVPos;
			}

			if (innerSurface->is_boundary(v))
				continue; // skip boundary vertices

			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);
			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const double epsilonCtrlWeight = -1.0 * GetSettings().InnerManifoldEpsilon(static_cast<double>(vMinDistance[v]));
			const auto negGradDotNormal = pmp::ddot(vNegGradDistanceToTarget, vNormal);
			const auto advectionDistance = static_cast<double>(GetSettings().AdvectionInteractWithOtherManifolds ? vMinDistance[v] : vDistanceToTarget[v]);
			const double etaCtrlWeight = GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal);

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * m_ExplicitLaplacianFunction(*innerSurface, v);

			// Tangential redistribution velocity
			const float tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * epsilonCtrlWeight;
			pmp::vec3 tanVelocity(0.0, 0.0, 0.0);
			if (tanRedistWeight > 0.0f)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*innerSurface, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::ExplicitIntegrationStep: innerSurface vertex " + std::to_string(v.idx()) + " outside m_EvolBox for time step id: "
					+ std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			innerSurface->position(v) = updatedPosition;
		}
	}
}


std::tuple<float, float, pmp::Point> ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		return { FLT_MAX, FLT_MAX, pmp::Point(0, 0, 0)};
	}

	const pmp::BoundingBox ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const float minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const float maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const float cellSize = minSize / static_cast<float>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	const SDF::PointCloudDistanceFieldSettings dfSettings{
		cellSize,
		GetSettings().FieldSettings.FieldExpansionFactor,
		DBL_MAX
	};
	m_DistanceField = std::make_shared<Geometry::ScalarGrid>(
		SDF::PointCloudDistanceFieldGenerator::Generate(*m_TargetPointCloud, dfSettings));
	m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid>(ComputeNormalizedNegativeGradient(*m_DistanceField));
	return { minSize, maxSize, ptCloudBBox.center() };
}

void ManifoldSurfaceEvolutionStrategy::ComputeVariableDistanceFields()
{
	if (!m_OuterSurface || m_InnerSurfaces.empty())
	{
		// there's no possibility of interaction between the outer and the inner manifolds
		return;
	}

	const SDF::DistanceFieldSettings surfaceDFSettings{
		GetFieldCellSize(),
		GetSettings().FieldSettings.FieldExpansionFactor,
		DBL_MAX,
		SDF::KDTreeSplitType::Center,
		SDF::SignComputation::None,
		SDF::BlurPostprocessingType::None,
		SDF::PreprocessingType::Octree
	};
	const Geometry::PMPSurfaceMeshAdapter outerSurfaceAdapter(std::make_shared<pmp::SurfaceMesh>(*m_OuterSurface));
	m_OuterSurfaceDistanceField = std::make_shared<Geometry::ScalarGrid>(
		SDF::DistanceFieldGenerator::Generate(outerSurfaceAdapter, surfaceDFSettings));

	for (const auto& innerSurface : m_InnerSurfaces)
	{
		const Geometry::PMPSurfaceMeshAdapter innerSurfaceAdapter(std::make_shared<pmp::SurfaceMesh>(*innerSurface));
		m_InnerSurfacesDistanceFields.emplace_back(std::make_shared<Geometry::ScalarGrid>(
			SDF::DistanceFieldGenerator::Generate(innerSurfaceAdapter, surfaceDFSettings)));
	}
}

void ManifoldSurfaceEvolutionStrategy::PrepareManifoldProperties()
{
	if (m_OuterSurface)
	{
		// distance to m_TargetPointCloud
		m_OuterSurface->add_vertex_property<pmp::Scalar>("v:distance_to_target", FLT_MAX);

		// the minimum from the distances to m_TargetPointCloud and all m_InnerSurfaces
		m_OuterSurface->add_vertex_property<pmp::Scalar>("v:min_distance", FLT_MAX);

		// the minimum from the distances to m_InnerSurfaces.
		m_OuterSurface->add_vertex_property<pmp::Scalar>("v:min_distance_to_inner", FLT_MAX);
	}

	for (const auto& innerSurface : m_InnerSurfaces)
	{
		// distance to m_TargetPointCloud
		innerSurface->add_vertex_property<pmp::Scalar>("v:distance_to_target", FLT_MAX);

		// the minimum from the distances to m_TargetPointCloud and m_OuterSurface
		innerSurface->add_vertex_property<pmp::Scalar>("v:min_distance", FLT_MAX);

		// distance to m_OuterSurface
		innerSurface->add_vertex_property<pmp::Scalar>("v:distance_to_outer", FLT_MAX);
	}

	// v:normal will be added during normal computation: pmp::Normals::compute_vertex_normals
}

void ManifoldSurfaceEvolutionStrategy::ResetManifoldProperties()
{
	// TODO: find a better way to reset prop values globally, perhaps implement a new property method
	if (m_OuterSurface)
	{
		auto vDistanceToTarget = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vMinDistanceToInner = m_OuterSurface->get_vertex_property<pmp::Scalar>("v:min_distance_to_inner");

		for (const auto v : m_OuterSurface->vertices())
		{
			vDistanceToTarget[v] = FLT_MAX;
			vMinDistance[v] = FLT_MAX;
			vMinDistanceToInner[v] = FLT_MAX;
		}
	}

	for (const auto& innerSurface : m_InnerSurfaces)
	{
		auto vDistanceToTarget = innerSurface->get_vertex_property<pmp::Scalar>("v:distance_to_target");
		auto vMinDistance = innerSurface->get_vertex_property<pmp::Scalar>("v:min_distance");
		auto vDistanceToOuter = innerSurface->get_vertex_property<pmp::Scalar>("v:distance_to_outer");

		for (const auto v : innerSurface->vertices())
		{
			vDistanceToTarget[v] = FLT_MAX;
			vMinDistance[v] = FLT_MAX;
			vDistanceToOuter[v] = FLT_MAX;
		}
	}
}

void ManifoldSurfaceEvolutionStrategy::ConstructInitialManifolds(float minTargetSize, float maxTargetSize, const pmp::Point& targetBoundsCenter)
{
	if (!GetSettings().UseInnerManifolds && !GetSettings().UseOuterManifolds)
	{
		throw std::invalid_argument("ManifoldSurfaceEvolutionStrategy::ConstructInitialManifolds: Current setting is: UseInnerManifolds == false && UseOuterManifolds == false. This means there's nothing to evolve!\n");
	}

	//const pmp::mat4 transfMatrixGeomMove{
	//	1.0f, 0.0f, 0.0f, targetBoundsCenter[0],
	//	0.0f, 1.0f, 0.0f, targetBoundsCenter[1],
	//	0.0f, 0.0f, 1.0f, targetBoundsCenter[2],
	//	0.0f, 0.0f, 0.0f, 1.0f
	//};
	const float outerSphereRadius = 0.5f * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5f + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);

	if (GetSettings().UseOuterManifolds)
	{
		Geometry::IcoSphereBuilder icoBuilder({ GetSettings().LevelOfDetail, outerSphereRadius });
		icoBuilder.BuildBaseData();
		icoBuilder.BuildPMPSurfaceMesh();
		m_OuterSurface = std::make_shared<pmp::SurfaceMesh>(icoBuilder.GetPMPSurfaceMeshResult());
		//(*m_OuterSurface) *= transfMatrixGeomMove; // center to target bounds
		m_InitialSphereSettings[m_OuterSurface.get()] = Sphere3D{ targetBoundsCenter, outerSphereRadius };
	}

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return;

	// TODO: This is hardcoded, implement 3D version of ParticleSwarmDistanceFieldInscribedSphereCalculator
	//const InscribedSphereInputData calcData{
	//	*m_TargetPointCloud,
	//	std::make_shared<Geometry::ScalarGrid>(*m_DistanceField) // clone
	//};	
	//ParticleSwarmDistanceFieldInscribedSphereCalculator inscribedSphereCalculator;
	//const auto spheres = inscribedSphereCalculator.Calculate(calcData);
	const std::vector spheres = {Sphere3D{pmp::Point(0, 0, 0), 0.7f}};

	m_InnerSurfaces.reserve(spheres.size());

	for (const auto& sphere : spheres)
	{
		// keep the same vertex density for inner spheres
		const auto innerSubdiv = static_cast<unsigned int>(static_cast<pmp::Scalar>(GetSettings().LevelOfDetail) * (sphere.Radius * SPHERE_RADIUS_FACTOR) / outerSphereRadius);
		Geometry::IcoSphereBuilder innerIcoBuilder({ 
			innerSubdiv,
			sphere.Radius * SPHERE_RADIUS_FACTOR
		});
		innerIcoBuilder.BuildBaseData();
		innerIcoBuilder.BuildPMPSurfaceMesh();
		auto mesh = innerIcoBuilder.GetPMPSurfaceMeshResult();
		if (sphere.Center != pmp::Point(0, 0, 0))
		{
			const auto translationMatrix = translation_matrix(sphere.Center);
			mesh *= translationMatrix;
		}		
		m_InnerSurfaces.push_back(std::make_shared<pmp::SurfaceMesh>(mesh));
		m_InitialSphereSettings[m_InnerSurfaces.back().get()] = sphere;
	}
}

/// \brief The power of the stabilizing scale factor.
constexpr float SCALE_FACTOR_POWER_2D = 1.0f / 2.0f;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr float INV_SHRINK_FACTOR_2D = 5.0f;

void ManifoldSurfaceEvolutionStrategy::StabilizeGeometries(float stabilizationFactor)
{
	const auto radius = stabilizationFactor * m_InitialSphereSettings.MinRadius() + (1.0f - stabilizationFactor) * m_InitialSphereSettings.MaxRadius();
	if (radius <= 0.0f)
	{
		throw std::invalid_argument("ManifoldSurfaceEvolutionStrategy::StabilizeGeometries: m_InitialSphereSettings empty!\n");
	}
	const unsigned int expectedVertexCount = (N_ICO_EDGES_0 * static_cast<unsigned int>(pow(4, GetSettings().LevelOfDetail) - 1) + 3 * N_ICO_VERTS_0) / 3;
	const float expectedMeanCoVolArea = (4.0f * static_cast<float>(M_PI) * radius * radius / static_cast<float>(expectedVertexCount));
	const auto scalingFactor = pow(static_cast<float>(GetSettings().TimeStep) / expectedMeanCoVolArea * INV_SHRINK_FACTOR_2D, SCALE_FACTOR_POWER_2D);
	GetScalingFactor() = scalingFactor;
	GetSettings().FieldSettings.FieldIsoLevel *= scalingFactor;

	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f, 0.0f,
		0.0f, 0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	const pmp::Point origin = m_OuterSurface ? m_InitialSphereSettings[m_OuterSurface.get()].Center :
		(!m_InnerSurfaces.empty() ? m_InitialSphereSettings[m_InnerSurfaces[0].get()].Center : pmp::Point{0, 0, 0});
	const pmp::mat4 transfMatrixGeomMove{
		1.0f, 0.0f, 0.0f, -origin[0],
		0.0f, 1.0f, 0.0f, -origin[1],
		0.0f, 0.0f, 1.0f, -origin[2],
		0.0f, 0.0f, 0.0f, 1.0f
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	m_TransformToOriginal = inverse(transfMatrixFull);

	// transform geometries
	if (m_OuterSurface)
	{
		(*m_OuterSurface) *= transfMatrixGeomScale;		
	}
	for (const auto& innerSurface : m_InnerSurfaces)
	{
		(*innerSurface) *= transfMatrixGeomScale;
	}

	// test box for geometry validation
	const float evolBoxFactor = 1.2f * scalingFactor;
	m_EvolBox = pmp::BoundingBox(
		pmp::Point{ -radius, -radius, -radius } * evolBoxFactor,
		pmp::Point{ radius, radius, radius } * evolBoxFactor);

	if (m_OuterSurfaceDistanceField)
	{
		(*m_OuterSurfaceDistanceField) *= transfMatrixFull;
		(*m_OuterSurfaceDistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerSurfaceDF : m_InnerSurfacesDistanceFields)
	{
		(*innerSurfaceDF) *= transfMatrixFull;
		(*innerSurfaceDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}

	if (!m_DistanceField || !m_DFNegNormalizedGradient)
		return; // nothing to scale

	(*m_DistanceField) *= transfMatrixFull;
	(*m_DFNegNormalizedGradient) *= transfMatrixFull;
	(*m_DistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
}

void ManifoldSurfaceEvolutionStrategy::AssignRemeshingSettingsToEvolvingManifolds()
{
	if (m_OuterSurface)
	{
		const Sphere3D& sphereSettings = m_InitialSphereSettings[m_OuterSurface.get()];
		m_RemeshingSettings[m_OuterSurface.get()] = CollectRemeshingSettingsFromIcoSphere_OLD(
			GetSettings().LevelOfDetail,
			sphereSettings.Radius * GetScalingFactor(),
			GetSettings().RemeshingSettings);
	}

	for (const auto& innerSurface : m_InnerSurfaces)
	{
		const Sphere3D& sphereSettings = m_InitialSphereSettings[innerSurface.get()];
		m_RemeshingSettings[innerSurface.get()] = CollectRemeshingSettingsFromIcoSphere_OLD(
			GetSettings().LevelOfDetail,
			sphereSettings.Radius * GetScalingFactor(),
			GetSettings().RemeshingSettings);
	}
}

void CustomManifoldSurfaceEvolutionStrategy::AssignRemeshingSettingsToEvolvingManifolds()
{
	if (GetOuterSurface())
	{
		GetRemeshingSettings()[GetOuterSurface().get()] = CollectRemeshingSettingsFromMesh(GetOuterSurface());
	}

	for (const auto& innerSurface : GetInnerSurfaces())
	{
		GetRemeshingSettings()[innerSurface.get()] = CollectRemeshingSettingsFromMesh(innerSurface);
	}
}

bool CustomManifoldSurfaceEvolutionStrategy::HasValidInnerOuterManifolds() const
{
	// check self-intersections
	if (GetOuterSurface() && Geometry::PMPSurfaceMeshHasSelfIntersections(*GetOuterSurface()))
		return false;

	const auto& outerSurface = *GetOuterSurface();
	const Geometry::PMPSurfaceMeshAdapter outerSurfaceAdapter(std::make_shared<pmp::SurfaceMesh>(outerSurface));
	const auto outerSurfaceKdTree = std::make_shared<Geometry::CollisionKdTree>(outerSurfaceAdapter, Geometry::CenterSplitFunction);
	// check inner surfaces
	for (const auto& innerSurface : GetInnerSurfaces())
	{
		// check self-intersections
		if (Geometry::PMPSurfaceMeshHasSelfIntersections(*innerSurface))
			return false;

		// check whether the inner surface is contained inside the outer surface.
		for (const auto& p : innerSurface->positions())
		{
			if (!IsPointInsidePMPSurfaceMesh(p, outerSurfaceKdTree))
				return false;
		}
	}

	return true;
}

std::pair<float, float> CustomManifoldSurfaceEvolutionStrategy::CalculateCoVolumeRange() const
{
	float minCoVolArea = FLT_MAX;
	float maxCoVolArea = -FLT_MAX;

	if (GetOuterSurface())
	{
		const auto& outerSurface = *GetOuterSurface();
		for (const auto v : outerSurface.vertices())
		{
			const auto currentCoVolArea = GetLaplacianAreaFunction()(outerSurface, v);

			if (currentCoVolArea > maxCoVolArea) maxCoVolArea = currentCoVolArea;
			if (currentCoVolArea < minCoVolArea) minCoVolArea = currentCoVolArea;
		}
	}

	for (const auto& surface : GetInnerSurfaces())
	{
		const auto& innerSurface = *surface;
		for (const auto v : innerSurface.vertices())
		{
			const auto currentCoVolArea = GetLaplacianAreaFunction()(innerSurface, v);

			if (currentCoVolArea > maxCoVolArea) maxCoVolArea = currentCoVolArea;
			if (currentCoVolArea < minCoVolArea) minCoVolArea = currentCoVolArea;
		}
	}

	return { minCoVolArea, maxCoVolArea };
}

void CustomManifoldSurfaceEvolutionStrategy::StabilizeCustomGeometries(float minArea, float maxArea, float stabilizationFactor)
{
	const float expectedMeanCoVolLength = (1.0f - stabilizationFactor) * minArea + stabilizationFactor * maxArea;
	const float scalingFactor = pow(static_cast<float>(GetSettings().TimeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_1D, SCALE_FACTOR_POWER_1D);
	GetScalingFactor() = scalingFactor;
	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0f, 0.0f, 0.0f,
		0.0f, scalingFactor, 0.0f, 0.0f,
		0.0f, 0.0f, scalingFactor, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	};
	GetTransformToOriginal() = inverse(transfMatrixGeomScale);

	// Transform the geometries
	(*GetOuterSurface()) *= transfMatrixGeomScale;
	for (auto& innerSurface : GetInnerSurfaces())
	{
		(*innerSurface) *= transfMatrixGeomScale;
	}

	if (GetOuterSurfraceDistanceField())
	{
		(*GetOuterSurfraceDistanceField()) *= transfMatrixGeomScale;
		(*GetOuterSurfraceDistanceField()) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerSurfaceDF : GetInnerSurfacesDistanceFields())
	{
		(*innerSurfaceDF) *= transfMatrixGeomScale;
		(*innerSurfaceDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}

	if (!GetDistanceField() || !GetDFNegNormalizedGradient())
		return; // nothing to scale

	(*GetDistanceField()) *= transfMatrixGeomScale;
	GetEvolBox() = GetDistanceField()->Box(); // test box for geometry validation
	(*GetDFNegNormalizedGradient()) *= transfMatrixGeomScale;
	(*GetDistanceField()) *= static_cast<double>(scalingFactor); // Scale the distance values
}

// ================================================
