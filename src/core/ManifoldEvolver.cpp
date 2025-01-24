#include "ManifoldEvolver.h"

#include "pmp/algorithms/CurveFactory.h"
#include "pmp/algorithms/Features.h"
#include "pmp/algorithms/Normals.h"

#include "utils/NumericalUtils.h"

#include "geometry/GridUtil.h"
#include "geometry/IcoSphereBuilder.h"
#include "geometry/MeshAnalysis.h"
#include "geometry/GeometryIOUtils.h"

#include "sdf/SDF.h"

#include "ConversionUtils.h"

//
// ======================================================================================
//                    The strategy for 1D Curves in 2D space
// --------------------------------------------------------------------------------------
//

void ManifoldCurveEvolutionStrategy::InitLogger(const std::string& baseOutputFileName)
{
	if (!LogManifoldValues())
		return;

	m_Logger.Init(baseOutputFileName + "_log.json");
}

void ManifoldCurveEvolutionStrategy::InitNewTimeStepLog(unsigned int stepId)
{
	if (!LogManifoldValues())
		return;

	m_Logger.StartNewTimeStep(stepId); // Explicitly start a new logging step
}

void ManifoldCurveEvolutionStrategy::SaveLog(bool omitLastTimeStep)
{
	if (!LogManifoldValues())
		return;

	m_Logger.Save(omitLastTimeStep);
}

void ManifoldCurveEvolutionStrategy::Preprocess()
{
	const auto [minTargetSize, maxTargetSize, targetCenter] = ComputeAmbientFields();
	ConstructInitialManifolds(minTargetSize, maxTargetSize, targetCenter);

	GetFieldCellSize() = m_DistanceField ? m_DistanceField->CellSize() : minTargetSize / static_cast<pmp::Scalar>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	ComputeVariableDistanceFields();

	if (GetSettings().UseStabilizationViaScaling)
	{
		StabilizeGeometries();
	}
	ComputeControlFunctionsLowerBounds();

	AssignRemeshingSettingsToEvolvingManifolds();

	if (LogManifoldValues())
		InitializeArcLengthCalculation();
}

void ManifoldCurveEvolutionStrategy::Postprocess()
{
}

void CustomManifoldCurveEvolutionStrategy::Preprocess()
{
	if (!GetOuterCurve() && GetInnerCurves().empty())
		throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: There's nothing to evolve!\n");
	
	//if (!HasValidInnerOuterManifolds())
	//	throw std::invalid_argument("CustomManifoldCurveEvolutionStrategy::Preprocess: Invalid inner /outer manifold geometry! Not all custom inner curves are contained within the custom outer curve.\n");

	if (LogManifoldValues())
	{
		if (GetOuterCurve())
			GetLogger().AddManifold(GetOuterCurve().get());
		for (auto& innerCurve : GetInnerCurves())
		{
			if (!innerCurve)
				continue;

			GetLogger().AddManifold(innerCurve.get());
		}
	}

	if (NeedsFieldsCalculation())
	{
		GetEvolBox() = GetOuterCurve() ? GetOuterCurve()->bounds() : GetInnerCurves()[0]->bounds();
		const auto sizeVec = (GetEvolBox().max() - GetEvolBox().min()) * GetSettings().FieldSettings.FieldExpansionFactor;
		GetEvolBox().expand(sizeVec[0], sizeVec[1]);

		std::tie(std::ignore, std::ignore, std::ignore) = ComputeAmbientFields();

		const auto minTargetSize = std::min(sizeVec[0], sizeVec[1]);
		GetFieldCellSize() = GetDistanceField() ? GetDistanceField()->CellSize() : minTargetSize / static_cast<pmp::Scalar>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
		ComputeVariableDistanceFields();
	}

	//if (GetSettings().UseStabilizationViaScaling)
	{
		const auto [minLength, maxLength] = CalculateCoVolumeRange();

		std::cout << "Before stabilization: { minLength: " << minLength << ", maxLength: " << maxLength << "} ... vs ... timeStep: " << GetSettings().TimeStep << "\n";

		StabilizeCustomGeometries(minLength, maxLength);

		const auto [minLengthAfter, maxLengthAfter] = CalculateCoVolumeRange();
		std::cout << "After stabilization: { minLengthAfter: " << minLengthAfter << ", maxLengthAfter: " << maxLengthAfter << "} ... vs ... timeStep: " << GetSettings().TimeStep << "\n";
	}
	ComputeControlFunctionsLowerBounds();

	AssignRemeshingSettingsToEvolvingManifolds();
	if (LogManifoldValues())
		InitializeArcLengthCalculation();
}

void ManifoldCurveEvolutionStrategy::PerformEvolutionStep(unsigned int stepId)
{
	if (LogManifoldValues())
		InitNewTimeStepLog(stepId);

	GetIntegrate()(stepId);
}

void ManifoldCurveEvolutionStrategy::Remesh()
{
	for (auto* curveToRemesh : m_RemeshTracker.GetManifoldsToRemesh())
	{
		pmp::CurveRemeshing remesher(*curveToRemesh, m_ArcLengthCalculators[curveToRemesh]);
		remesher.adaptive_remeshing(m_RemeshingSettings[curveToRemesh]);

		//if (LogManifoldValues())
		//{
		//	m_Logger.ReserveBuffers(curveToRemesh);
		//}
	}
	m_RemeshTracker.Reset();
}

void ManifoldCurveEvolutionStrategy::ResizeRemeshingSettings(pmp::Scalar resizeFactor)
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
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportFinalResult: error writing " << (baseOutputFilename + "_Outer" + connectingName + ".ply") << "!\n";
	}

	for (size_t i = 0; const auto & innerCurve : m_InnerCurves)
	{
		auto exportedInnerCurve = *innerCurve;
		exportedInnerCurve *= m_TransformToOriginal;
		if (!write_to_ply(exportedInnerCurve, baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply"))
			std::cerr << "ManifoldCurveEvolutionStrategy::ExportFinalResult: error writing " << (baseOutputFilename + "_Inner" + std::to_string(i++) + connectingName + ".ply") << "!\n";
	}
}

void ManifoldCurveEvolutionStrategy::ExportTargetDistanceFieldAsImage(const std::string& baseOutputFilename)
{
	constexpr double colorMapPlotScaleFactor = 0.5; // scale the distance field color map down to show more detail

	// ----------------------------------------------------------------
	if (GetSettings().ExportVariableScalarFieldsDimInfo && m_OuterCurveDistanceField)
	{
		auto exportedOuterCurveDF = *m_OuterCurveDistanceField;
		exportedOuterCurveDF *= m_TransformToOriginal;
		exportedOuterCurveDF /= static_cast<double>(GetScalingFactor());
		ExportScalarGridDimInfo2D(baseOutputFilename + "_OuterDF.gdim2d", exportedOuterCurveDF);
		ExportScalarGrid2DToPNG(baseOutputFilename + "_OuterDF.png", exportedOuterCurveDF, m_ScalarInterpolate, 
			10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);
	}
	if (GetSettings().ExportVariableVectorFieldsDimInfo && m_OuterCurveDFNegNormalizedGradient)
	{
		auto exportedOuterCurveDFNegGrad = *m_OuterCurveDFNegNormalizedGradient;
		exportedOuterCurveDFNegGrad *= m_TransformToOriginal;
		ExportVectorGridDimInfo2D(baseOutputFilename + "_OuterDFNegGrad.gdim2d", exportedOuterCurveDFNegGrad);
	}

	if (GetSettings().ExportVariableScalarFieldsDimInfo)
	{
		for (size_t i = 0; i < m_InnerCurvesDistanceFields.size(); ++i)
		{
			if (!m_InnerCurvesDistanceFields[i])
				continue;

			auto exportedInnerCurveDF = *m_InnerCurvesDistanceFields[i];
			exportedInnerCurveDF *= m_TransformToOriginal;
			exportedInnerCurveDF /= static_cast<double>(GetScalingFactor());
			ExportScalarGridDimInfo2D(baseOutputFilename + "_InnerDF" + std::to_string(i) + ".gdim2d", exportedInnerCurveDF);
			ExportScalarGrid2DToPNG(baseOutputFilename + "_InnerDF" + std::to_string(i) + ".png", exportedInnerCurveDF, m_ScalarInterpolate, 
				10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);
		}
	}
	if (GetSettings().ExportVariableVectorFieldsDimInfo)
	{
		for (size_t i = 0; i < m_InnerCurvesDFNegNormalizedGradients.size(); ++i)
		{
			if (!m_InnerCurvesDFNegNormalizedGradients[i])
				continue;

			auto exportedInnerCurveDFNegGrad = *m_InnerCurvesDFNegNormalizedGradients[i];
			exportedInnerCurveDFNegGrad *= m_TransformToOriginal;
			ExportVectorGridDimInfo2D(baseOutputFilename + "_InnerDFNegGrad" + std::to_string(i) + ".gdim2d", exportedInnerCurveDFNegGrad);
		}
	}
	// ----------------------------------------------------------------

	if (!m_DistanceField)
		return; // field not defined

	auto exportedField = *m_DistanceField;
	exportedField *= m_TransformToOriginal;
	exportedField /= static_cast<double>(GetScalingFactor());

	ExportScalarGrid2DToPNG(baseOutputFilename + "_TargetDF.png", exportedField, m_ScalarInterpolate, 
		10, 10, RAINBOW_TO_WHITE_MAP * colorMapPlotScaleFactor);
	const auto [dfMin, dfMax] = std::ranges::minmax_element(m_DistanceField->Values());
	constexpr double colorLegendScaleFactor = 2.0; // stretch color map keys (ratios) to fit the gradient to the full legend bar
	constexpr double dfValueClampFactor = colorMapPlotScaleFactor / colorLegendScaleFactor; // value cutoff for the chosen color legend.
	constexpr bool verticalLegend = true;
	constexpr unsigned int resolutionFactor = 4;
	constexpr unsigned int legendPxHeight = (verticalLegend ? 400 : 100) * resolutionFactor;
	constexpr unsigned int legendPxWidth = (verticalLegend ? 100 : 600) * resolutionFactor;
	ExportColorScaleToPNG(
		baseOutputFilename + "_TargetDF_Scale.png",
		(*dfMin) / GetScalingFactor(),
		dfValueClampFactor * (*dfMax) / GetScalingFactor(),
		RAINBOW_TO_WHITE_MAP * colorLegendScaleFactor, 
		legendPxHeight, legendPxWidth);
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

		const auto tStep = GetSettings().TimeStep;

		pmp::Normals2::compute_vertex_normals(*m_OuterCurve);
		auto vNormalsProp = m_OuterCurve->get_vertex_property<pmp::vec2>("v:normal");

		// prepare matrix & rhs for m_OuterCurve:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 2);

		std::vector<pmp::Scalar> arcLengths;
		if (LogOuterManifoldValues() && m_ArcLengthCalculators[m_OuterCurve.get()])
		{
			arcLengths = m_ArcLengthCalculators[m_OuterCurve.get()]->CalculateArcLengths();
		}

		for (const auto v : m_OuterCurve->vertices())
		{
			const auto vPosToUpdate = m_OuterCurve->position(v);

			InteractionDistanceInfo<pmp::dvec2> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec2>{vDistanceToTarget, vNegGradDistanceToTarget};

			double vMinDistanceToInner = DBL_MAX;
			for (unsigned int i = 0; i < m_InnerCurvesDistanceFields.size(); ++i)
			{
				if (!m_InnerCurvesDistanceFields[i] || !m_InnerCurvesDFNegNormalizedGradients[i])
					continue;

				auto innerDfAtVPos = m_ScalarInterpolate(vPosToUpdate, *m_InnerCurvesDistanceFields[i]);
				innerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToInner = m_VectorInterpolate(vPosToUpdate, *m_InnerCurvesDFNegNormalizedGradients[i]);

				interaction << InteractionDistanceInfo<pmp::dvec2>{innerDfAtVPos, vNegGradDistanceToInner};

				if (innerDfAtVPos < vMinDistanceToInner)
					vMinDistanceToInner = innerDfAtVPos;
			}

			if (m_OuterCurve->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector2d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const double epsilonCtrlWeight =
				GetSettings().OuterManifoldEpsilon(static_cast<double>(interaction.Distance));

			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				(GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget), vNormal), -1.0, 1.0);
			const double advectionDistance = 
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight = 
				((m_DistanceField || !m_InnerCurvesDistanceFields.empty()) ? GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal) : 0.0) +
				GetSettings().OuterManifoldRepulsion(static_cast<double>(vMinDistanceToInner));

			if (LogOuterManifoldValues())
			{
				if (GetSettings().DiagSettings.LogOuterManifoldEpsilon)
					m_Logger.LogValue(m_OuterCurve.get(), "epsilonCtrlWeight", v.idx(), epsilonCtrlWeight);
				if (GetSettings().DiagSettings.LogOuterManifoldEta)
					m_Logger.LogValue(m_OuterCurve.get(), "etaCtrlWeight", v.idx(), etaCtrlWeight);
				if (!arcLengths.empty())
					m_Logger.LogValue(m_OuterCurve.get(), "arcLength", v.idx(), arcLengths[v.idx()]);
			}

			const Eigen::Vector2d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;

			//const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * 0.0; // * (1.0 - exp(-interaction.Distance * interaction.Distance / 1.0));
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight); // (1.0 - exp(-interaction.Distance * interaction.Distance));

			if (tanRedistWeight > 0.0)
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
		if (IsRemeshingNecessary(*m_OuterCurve, m_RemeshingSettings[m_OuterCurve.get()]))
			m_RemeshTracker.AddManifold(m_OuterCurve.get());

		// solve
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
		Eigen::MatrixXd x = solver.solve(sysRhs);
		if (solver.info() != Eigen::Success)
		{
			PrintSparseMatrixAndRHS(sysMat, sysRhs, "outerCurve Sys");

			const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(step) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}

		// update vertex positions & verify mesh within bounds
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);

			if (LogOuterManifoldValues())
			{
				if (GetSettings().DiagSettings.LogOuterManifoldErrors)
				{
					const pmp::vec2 remainderVec{
						Utils::Get32BitRemainder(newPos[0]),
						Utils::Get32BitRemainder(newPos[1])
					};
					m_Logger.LogValue(m_OuterCurve.get(), "outerCurve32RemainderMagnitude", i, norm(remainderVec));
				}
				else if (GetSettings().DiagSettings.LogOuterManifoldXErrors)
				{
					m_Logger.LogValue(m_OuterCurve.get(), "outerCurve32RemainderX", i, Utils::Get32BitRemainder(newPos[0]));
				}
				else if (GetSettings().DiagSettings.LogOuterManifoldYErrors)
				{
					m_Logger.LogValue(m_OuterCurve.get(), "outerCurve32RemainderY", i, Utils::Get32BitRemainder(newPos[1]));
				}				
			}

			if (!m_EvolBox.Contains(newPos))
			{
				const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: vertex " + std::to_string(i) 
					+ ": (" + std::to_string(newPos[0]) + ", " + std::to_string(newPos[1]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
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

		const auto tStep = GetSettings().TimeStep;

		pmp::Normals2::compute_vertex_normals(*innerCurve);
		auto vNormalsProp = innerCurve->get_vertex_property<pmp::vec2>("v:normal");

		// prepare matrix & rhs for m_OuterCurve:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 2);  // Assuming 2 entries per vertex for curves

		std::vector<pmp::Scalar> arcLengths;
		if (LogInnerManifoldValues() && m_ArcLengthCalculators[innerCurve.get()])
		{
			arcLengths = m_ArcLengthCalculators[innerCurve.get()]->CalculateArcLengths();
		}

		for (const auto v : innerCurve->vertices())
		{
			const auto vPosToUpdate = innerCurve->position(v);

			InteractionDistanceInfo<pmp::dvec2> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec2>{vDistanceToTarget, vNegGradDistanceToTarget};

			double outerDfAtVPos = DBL_MAX;
			if (m_OuterCurveDistanceField && m_OuterCurveDFNegNormalizedGradient)
			{
				outerDfAtVPos = m_ScalarInterpolate(vPosToUpdate, *m_OuterCurveDistanceField);
				outerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToOuter = m_VectorInterpolate(vPosToUpdate, *m_OuterCurveDFNegNormalizedGradient);

				interaction << InteractionDistanceInfo<pmp::dvec2>{outerDfAtVPos, vNegGradDistanceToOuter};
			}

			if (innerCurve->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector2d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const double epsilonCtrlWeight =
				GetSettings().InnerManifoldEpsilon(static_cast<double>(interaction.Distance));

			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget, vNormal), -1.0, 1.0);
			const double advectionDistance = 
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight = 
				((m_DistanceField || m_OuterCurveDistanceField) ? GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal) : 0.0) +
				GetSettings().InnerManifoldRepulsion(static_cast<double>(outerDfAtVPos));

			if (LogInnerManifoldValues())
			{
				if (GetSettings().DiagSettings.LogInnerManifoldsEpsilon)
					m_Logger.LogValue(innerCurve.get(), "epsilonCtrlWeight", v.idx(), epsilonCtrlWeight);
				if (GetSettings().DiagSettings.LogInnerManifoldsEta)
					m_Logger.LogValue(innerCurve.get(), "etaCtrlWeight", v.idx(), etaCtrlWeight);
				if (!arcLengths.empty())
					m_Logger.LogValue(innerCurve.get(), "arcLength", v.idx(), arcLengths[v.idx()]);
			}

			const Eigen::Vector2d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;

			//const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * 0.0; // (1.0 - exp(-interaction.Distance * interaction.Distance / 1.0));
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight); //(1.0 - exp(-interaction.Distance * interaction.Distance));

			if (tanRedistWeight > 0.0)
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

		////PrintSparseMatrix(sysMat);
		//if (!IsDiagonallyDominant(sysMat))
		//{
		//	constexpr double lambda = 0.01; // Example regularization parameter
		//	RegularizeMatrixInPlace(sysMat, lambda);
		//	//std::cout << "step " << step << ": sysMat(" << NVertices << ", " << NVertices << ") is NOT diagonally-dominant!\n";
		//}

		// solve
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> solver(sysMat);
		Eigen::MatrixXd x = solver.solve(sysRhs);
		if (solver.info() != Eigen::Success)
		{
			PrintSparseMatrixAndRHS(sysMat, sysRhs, "innerCurve Sys");

			const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: solver.info() != Eigen::Success for time step id: "
				+ std::to_string(step) + ", Error code: " + InterpretSolverErrorCode(solver.info()) + "\n";
			std::cerr << msg;
			throw std::runtime_error(msg);
		}

		// update vertex positions & verify mesh within bounds
		//bool hasProblematicVertices = false;
		for (unsigned int i = 0; i < NVertices; i++)
		{
			const auto newPos = x.row(i);

			if (LogInnerManifoldValues())
			{
				if (GetSettings().DiagSettings.LogInnerManifoldsErrors)
				{
					const pmp::vec2 remainderVec{
						Utils::Get32BitRemainder(newPos[0]),
						Utils::Get32BitRemainder(newPos[1])
					};
					m_Logger.LogValue(innerCurve.get(), "innerCurve32RemainderMagnitude", i, norm(remainderVec));
				}
				else if (GetSettings().DiagSettings.LogInnerManifoldsXErrors)
				{
					m_Logger.LogValue(innerCurve.get(), "innerCurve32RemainderX", i, Utils::Get32BitRemainder(newPos[0]));
				}
				else if (GetSettings().DiagSettings.LogInnerManifoldsYErrors)
				{
					m_Logger.LogValue(innerCurve.get(), "innerCurve32RemainderY", i, Utils::Get32BitRemainder(newPos[1]));
				}				
			}

			if (!m_EvolBox.Contains(newPos))
			{
				//const auto newPosInOrigScale = affine_transform(m_TransformToOriginal, pmp::Point2(newPos));
				//std::cout << "\nv" << i << " = (" + std::to_string(newPosInOrigScale[0]) + ", " + std::to_string(newPosInOrigScale[1]) + ") outside bounds!\n";
				//hasProblematicVertices = true;
				const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: innerCurve vertex " + std::to_string(i)
					+ ": (" + std::to_string(newPos[0]) + ", " + std::to_string(newPos[1]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			//if (!hasProblematicVertices)
			//	innerCurve->position(pmp::Vertex(i)) = newPos;
			innerCurve->position(pmp::Vertex(i)) = newPos;
		}

		//if (hasProblematicVertices)
		//{
		//	const std::string msg = "\nManifoldCurveEvolutionStrategy::SemiImplicitIntegrationStep: found problematic vertices for time step id: " + std::to_string(step) + "!\n";
		//	std::cerr << msg;
		//	throw std::runtime_error(msg);
		//}
	}
}

void ManifoldCurveEvolutionStrategy::ExplicitIntegrationStep(unsigned int step)
{
	// ================================== Handle m_OuterCurve ==========================================================
	if (m_OuterCurve)
	{
		const auto tStep = GetSettings().TimeStep;

		pmp::Normals2::compute_vertex_normals(*m_OuterCurve);
		auto vNormalsProp = m_OuterCurve->get_vertex_property<pmp::vec2>("v:normal");

		for (const auto v : m_OuterCurve->vertices())
		{
			const auto vPosToUpdate = m_OuterCurve->position(v);

			InteractionDistanceInfo<pmp::dvec2> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec2>{vDistanceToTarget, vNegGradDistanceToTarget};

			double vMinDistanceToInner = DBL_MAX;
			for (unsigned int i = 0; i < m_InnerCurvesDistanceFields.size(); ++i)
			{
				if (!m_InnerCurvesDistanceFields[i] || !m_InnerCurvesDFNegNormalizedGradients[i])
					continue;

				auto innerDfAtVPos = static_cast<pmp::Scalar>(m_ScalarInterpolate(vPosToUpdate, *m_InnerCurvesDistanceFields[i]));
				innerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToInner = m_VectorInterpolate(vPosToUpdate, *m_InnerCurvesDFNegNormalizedGradients[i]);

				interaction << InteractionDistanceInfo<pmp::dvec2>{innerDfAtVPos, vNegGradDistanceToInner};

				if (innerDfAtVPos < vMinDistanceToInner)
					vMinDistanceToInner = innerDfAtVPos;
			}

			if (m_OuterCurve->is_boundary(v))
				continue; // skip boundary vertices

			const double epsilonCtrlWeight =
				GetSettings().OuterManifoldEpsilon(static_cast<double>(interaction.Distance)) +
				GetSettings().OuterManifoldRepulsion(static_cast<double>(vMinDistanceToInner));

			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				(GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget), vNormal), -1.0, 1.0);
			const double advectionDistance =
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight =
				(m_DistanceField || !m_InnerCurvesDistanceFields.empty()) ? GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal) : 0.0;

			if (LogOuterManifoldValues())
			{
				if (GetSettings().DiagSettings.LogOuterManifoldEpsilon)
					m_Logger.LogValue(m_OuterCurve.get(), "epsilonCtrlWeight", v.idx(), epsilonCtrlWeight);
				if (GetSettings().DiagSettings.LogOuterManifoldEta)
					m_Logger.LogValue(m_OuterCurve.get(), "etaCtrlWeight", v.idx(), etaCtrlWeight);
			}

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * pmp::laplace_1D(*m_OuterCurve, v);

			// Tangential redistribution velocity
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight);
			pmp::vec2 tanVelocity(0.0, 0.0);
			if (tanRedistWeight > 0.0)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_OuterCurve, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldCurveEvolutionStrategy::ExplicitIntegrationStep: vertex " + std::to_string(v.idx()) 
					+ ": (" + std::to_string(updatedPosition[0]) + ", " + std::to_string(updatedPosition[1]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
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

		pmp::Normals2::compute_vertex_normals(*innerCurve);
		auto vNormalsProp = innerCurve->get_vertex_property<pmp::vec2>("v:normal");

		for (const auto v : innerCurve->vertices())
		{
			const auto vPosToUpdate = innerCurve->position(v);

			InteractionDistanceInfo<pmp::dvec2> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec2(0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec2>{vDistanceToTarget, vNegGradDistanceToTarget};

			double outerDfAtVPos = DBL_MAX;
			if (m_OuterCurveDistanceField && m_OuterCurveDFNegNormalizedGradient)
			{
				outerDfAtVPos = m_ScalarInterpolate(vPosToUpdate, *m_OuterCurveDistanceField);
				outerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToOuter = m_VectorInterpolate(vPosToUpdate, *m_OuterCurveDFNegNormalizedGradient);

				interaction << InteractionDistanceInfo<pmp::dvec2>{outerDfAtVPos, vNegGradDistanceToOuter};
			}

			if (innerCurve->is_boundary(v))
				continue; // skip boundary vertices
			
			const double epsilonCtrlWeight =
				GetSettings().InnerManifoldEpsilon(static_cast<double>(interaction.Distance)) +
				GetSettings().InnerManifoldRepulsion(static_cast<double>(outerDfAtVPos));

			const auto vNormal = static_cast<pmp::vec2>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget, vNormal), -1.0, 1.0);
			const double advectionDistance =
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight =
				(m_DistanceField || m_OuterCurveDistanceField) ? GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal) : 0.0;

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * pmp::laplace_1D(*innerCurve, v);

			// Tangential redistribution velocity
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight);
			pmp::vec2 tanVelocity(0.0, 0.0);
			if (tanRedistWeight > 0.0)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*innerCurve, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldCurveEvolutionStrategy::ExplicitIntegrationStep: innerCurve vertex " + std::to_string(v.idx()) 
					+ ": (" + std::to_string(updatedPosition[0]) + ", " + std::to_string(updatedPosition[1]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			innerCurve->position(v) = updatedPosition;
		}
	}
}

std::tuple<pmp::Scalar, pmp::Scalar, pmp::Point2> ManifoldCurveEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldCurveEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		return { FLT_MAX, FLT_MAX, pmp::Point2(0, 0)};
	}

	const pmp::BoundingBox2 ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const pmp::Scalar minSize = std::min(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const pmp::Scalar maxSize = std::max(ptCloudBBoxSize[0], ptCloudBBoxSize[1]);
	const auto cellSize = minSize / static_cast<pmp::Scalar>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	const SDF::PointCloudDistanceField2DSettings dfSettings{
		cellSize,
		GetSettings().FieldSettings.FieldExpansionFactor,
		DBL_MAX
	};
	m_DistanceField = std::make_shared<Geometry::ScalarGrid2D>(
		SDF::PlanarPointCloudDistanceFieldGenerator::Generate(*m_TargetPointCloud, dfSettings));
	RepairScalarGrid2D(*m_DistanceField);
	m_DFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid2D>(ComputeNormalizedNegativeGradient(*m_DistanceField));
	return { minSize, maxSize, ptCloudBBox.center() };
}

void ManifoldCurveEvolutionStrategy::ComputeVariableDistanceFields()
{
	if (!NeedsVariableFieldsCalculation())
	{
		// there's no possibility of interaction between the outer and the inner manifolds
		return;
	}

	// clear fields
	m_OuterCurveDistanceField.reset();
	m_OuterCurveDFNegNormalizedGradient.reset();
	m_InnerCurvesDistanceFields.clear();
	m_InnerCurvesDFNegNormalizedGradients.clear();

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
	RepairScalarGrid2D(*m_OuterCurveDistanceField);
	m_OuterCurveDFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid2D>(ComputeNormalizedNegativeGradient(*m_OuterCurveDistanceField));

	for (const auto& innerCurve : m_InnerCurves)
	{
		const Geometry::ManifoldCurve2DAdapter innerCurveAdapter(std::make_shared<pmp::ManifoldCurve2D>(*innerCurve));
		m_InnerCurvesDistanceFields.emplace_back(std::make_shared<Geometry::ScalarGrid2D>(
			SDF::PlanarDistanceFieldGenerator::Generate(innerCurveAdapter, curveDFSettings, m_OuterCurveDistanceField->Box())));
		RepairScalarGrid2D(*m_InnerCurvesDistanceFields.back());
		m_InnerCurvesDFNegNormalizedGradients.emplace_back(
			std::make_shared<Geometry::VectorGrid2D>(ComputeNormalizedNegativeGradient(*m_InnerCurvesDistanceFields.back())));
	}
}

void ManifoldCurveEvolutionStrategy::ConstructInitialManifolds(pmp::Scalar minTargetSize, pmp::Scalar maxTargetSize, const pmp::Point2& targetBoundsCenter)
{
	if (!GetSettings().UseInnerManifolds && !GetSettings().UseOuterManifolds)
	{
		throw std::invalid_argument("ManifoldCurveEvolutionStrategy::ConstructInitialManifolds: Current setting is: UseInnerManifolds == false && UseOuterManifolds == false. This means there's nothing to evolve!\n");
	}

	const pmp::Scalar outerCircleRadius = 0.5 * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5 + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);
	const auto nSegments = static_cast<unsigned int>(pow(2, GetSettings().LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;

	if (GetSettings().UseOuterManifolds)
	{
		m_OuterCurve = std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(pmp::Point2{}, outerCircleRadius, nSegments));
		// DISCLAIMER: Since we want to evolve manifolds centered at the origin, we will not move the outer manifold into its "true" position.
		// The "true" position will be stored in m_InitialSphereSettings:
		m_InitialSphereSettings[m_OuterCurve.get()] = Circle2D{ targetBoundsCenter, outerCircleRadius };

		if (LogManifoldValues())
			m_Logger.AddManifold(m_OuterCurve.get());
	}

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return;

	throw std::invalid_argument(
		"ManifoldCurveEvolutionStrategy::ConstructInitialManifolds: Automatic calculation of inscribed manifolds is not reliable yet! Please use CustomManifoldCurveEvolutionStrategy instead.");

	//const InscribedCircleInputData calcData{
	//	*m_TargetPointCloud,
	//	std::make_shared<Geometry::ScalarGrid2D>(*m_DistanceField) // clone
	//};
	//ParticleSwarmDistanceFieldInscribedCircleCalculator inscribedCircleCalculator;
	//const auto circles = inscribedCircleCalculator.Calculate(calcData);

	// Hardcoded inner curves:
	//const auto circles = std::vector{ Circle2D{targetBoundsCenter, 1.85} };

	//m_InnerCurves.reserve(circles.size());

	//for (const auto& circle : circles)
	//{
	//	// keep the same vertex density for inner circles
	//	const auto nInnerSegments = static_cast<unsigned int>(static_cast<pmp::Scalar>(nSegments) * (circle.Radius) / outerCircleRadius);
	//	m_InnerCurves.emplace_back(std::make_shared<pmp::ManifoldCurve2D>(pmp::CurveFactory::circle(
	//		circle.Center,
	//		circle.Radius * SPHERE_RADIUS_FACTOR, 
	//		nInnerSegments
	//	)));
	//	m_InitialSphereSettings[m_InnerCurves.back().get()] = circle;
	//}
}

/// \brief The power of the stabilizing scale factor.
constexpr pmp::Scalar SCALE_FACTOR_POWER_1D = 1.0;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr pmp::Scalar INV_SHRINK_FACTOR_1D = 20.0; // TODO: seems like an overkill. Investigate

void ManifoldCurveEvolutionStrategy::StabilizeGeometries(pmp::Scalar stabilizationFactor)
{
	const pmp::Scalar radius = stabilizationFactor * m_InitialSphereSettings.MinRadius() + (1.0 - stabilizationFactor) * m_InitialSphereSettings.MaxRadius();
	if (radius <= 0.0)
	{
		throw std::invalid_argument("ManifoldCurveEvolutionStrategy::StabilizeGeometries: m_InitialSphereSettings empty!\n");
	}
	const auto expectedVertexCount = static_cast<unsigned int>(pow(2, GetSettings().LevelOfDetail - 1)) * N_CIRCLE_VERTS_0;
	const auto expectedMeanCoVolLength = (2.0 * static_cast<pmp::Scalar>(M_PI) * radius / static_cast<pmp::Scalar>(expectedVertexCount));
	const auto scalingFactor = pow(static_cast<pmp::Scalar>(GetSettings().TimeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_1D, SCALE_FACTOR_POWER_1D);
	GetScalingFactor() = scalingFactor;

	// -----------------------------------------------------------------------------------------------
	// All geometric quantities need to be scaled by the scalingFactor to ensure numerical stability.
	// For spatial conveniences, geometries must also be centered, i.e.: translated by -origin.
	// Geometries that are already centered at (0, 0) are not translated.
	// -----------------------------------------------------------------------------------------------
	GetSettings().FieldSettings.FieldIsoLevel *= scalingFactor;
	GetFieldCellSize() *= scalingFactor;

	const pmp::mat3 transfMatrixGeomScale{
		(pmp::Scalar)scalingFactor, 0.0, 0.0,
		0.0, (pmp::Scalar)scalingFactor, 0.0,
		0.0, 0.0, 1.0
	};
	const pmp::Point2 origin = m_OuterCurve ? m_InitialSphereSettings[m_OuterCurve.get()].Center :
		(!m_InnerCurves.empty() ? m_InitialSphereSettings[m_InnerCurves[0].get()].Center : pmp::Point2{ 0, 0 });
	const pmp::mat3 transfMatrixGeomMove{
		1.0, 0.0, -origin[0],
		0.0, 1.0, -origin[1],
		0.0, 0.0, 1.0
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	m_TransformToOriginal = inverse(transfMatrixFull);
	m_ScaleFieldToOriginal = inverse(transfMatrixGeomScale);

	// transform geometries
	if (m_OuterCurve)
	{
		(*m_OuterCurve) *= transfMatrixGeomScale;
	}
	for (const auto& innerCurve : m_InnerCurves)
	{
		(*innerCurve) *= transfMatrixFull;
	}

	// test box for geometry validation
	const pmp::Scalar evolBoxFactor = 5.0 * scalingFactor;
	m_EvolBox = pmp::BoundingBox2(
		pmp::Point2{ -radius, -radius } * evolBoxFactor,
		pmp::Point2{ radius, radius } * evolBoxFactor);

	if (m_OuterCurveDistanceField)
	{
		(*m_OuterCurveDistanceField) *= transfMatrixGeomScale;
		(*m_OuterCurveDistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	if (m_OuterCurveDFNegNormalizedGradient)
	{
		(*m_OuterCurveDFNegNormalizedGradient) *= transfMatrixGeomScale;
		// values are supposed to be unit vectors regardless of scaling
	}
	for (const auto& innerCurveDF : m_InnerCurvesDistanceFields)
	{
		(*innerCurveDF) *= transfMatrixFull;
		(*innerCurveDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerCurveDFGradient : m_InnerCurvesDFNegNormalizedGradients)
	{
		(*innerCurveDFGradient) *= transfMatrixFull;
		// values are supposed to be unit vectors regardless of scaling
	}

	if (!m_DistanceField || !m_DFNegNormalizedGradient)
		return; // nothing to scale

	(*m_DistanceField) *= transfMatrixFull;
	(*m_DistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
	(*m_DFNegNormalizedGradient) *= transfMatrixFull;
	// values are supposed to be unit vectors regardless of scaling
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

void ManifoldCurveEvolutionStrategy::InitializeArcLengthCalculation()
{
	if (m_OuterCurve)
	{
		m_ArcLengthCalculators[m_OuterCurve.get()] = std::make_shared<pmp::EvolvingArcLengthCalculator>(*m_OuterCurve);
	}

	for (const auto& innerCurve : m_InnerCurves)
	{
		m_ArcLengthCalculators[innerCurve.get()] = std::make_shared<pmp::EvolvingArcLengthCalculator>(*innerCurve);
	}
}

bool ManifoldCurveEvolutionStrategy::NeedsVariableFieldsCalculation()
{
	return m_OuterCurve && !m_InnerCurves.empty();
}

bool ManifoldCurveEvolutionStrategy::NeedsFieldsCalculation()
{
	return NeedsVariableFieldsCalculation() || m_TargetPointCloud;
}

void ManifoldCurveEvolutionStrategy::ComputeControlFunctionsLowerBounds()
{
	// find reasonable bounds for distance
	double maxDistance{ 0.0 };
	if (m_DistanceField)
	{
		maxDistance = *std::ranges::max_element(m_DistanceField->Values());
	}
	else if (!m_DistanceField && m_OuterCurveDistanceField)
	{
		maxDistance = *std::ranges::max_element(m_OuterCurveDistanceField->Values());
	}
	else if (!m_DistanceField && !m_OuterCurveDistanceField && !m_InnerCurvesDistanceFields.empty())
	{
		maxDistance = *std::ranges::max_element(m_InnerCurvesDistanceFields.front()->Values());
	}

	// normalize the max distance to expansion factor
	const auto expFactor = GetSettings().FieldSettings.FieldExpansionFactor;
	if (expFactor < FLT_EPSILON)
	{
		throw std::invalid_argument("ManifoldCurveEvolutionStrategy::ComputeControlFunctionsLowerBounds: expFactor < FLT_EPSILON!\n");
	}
	maxDistance /= (0.5 + expFactor);

	// apply the limit
	GetSettings().OuterManifoldEpsilon.ComputeLimit(maxDistance);
	GetSettings().OuterManifoldEta.ComputeLimit(maxDistance);
	GetSettings().InnerManifoldEpsilon.ComputeLimit(maxDistance);
	GetSettings().InnerManifoldEta.ComputeLimit(maxDistance);
}

void CustomManifoldCurveEvolutionStrategy::AssignRemeshingSettingsToEvolvingManifolds()
{
	if (GetOuterCurve())
	{
		GetRemeshingSettings()[GetOuterCurve().get()] = 
			CollectRemeshingSettingsForManifold_EXPERIMENTAL(GetOuterCurve(), GetSettings().RemeshingSettings);
	}

	for (const auto& innerCurve : GetInnerCurves())
	{
		GetRemeshingSettings()[innerCurve.get()] = 
			CollectRemeshingSettingsForManifold_EXPERIMENTAL(innerCurve, GetSettings().RemeshingSettings);
	}
}

bool CustomManifoldCurveEvolutionStrategy::HasValidInnerOuterManifolds() const
{
	// check self-intersections
	if (GetOuterCurve() && Geometry::PMPManifoldCurve2DHasSelfIntersections(*GetOuterCurve()))
		return false;

	// check inner curves
	for (const auto& innerCurve : GetInnerCurves())
	{
		// check self-intersections
		if (Geometry::PMPManifoldCurve2DHasSelfIntersections(*innerCurve))
			return false;

		if (!GetOuterCurve())
			continue;

		// check whether the inner curve is contained inside the outer curve.
		//const auto& outerCurve = *GetOuterCurve();
		//for (const auto& p : innerCurve->positions())
		//{
		//	if (!Geometry::IsPointInsidePMPManifoldCurve(p, outerCurve))
		//		return false;
		//}
	}

	return true;
}

std::pair<pmp::Scalar, pmp::Scalar> CustomManifoldCurveEvolutionStrategy::CalculateCoVolumeRange() const
{
	pmp::Scalar minCoVolLength = FLT_MAX;
	pmp::Scalar maxCoVolLength = -FLT_MAX;

	if (GetOuterCurve())
	{
		const auto& outerCurve = *GetOuterCurve();
		for (const auto v : outerCurve.vertices())
		{
			const auto [eTo, eFrom] = outerCurve.edges(v);
			const auto currentCoVolLength = 0.5 * (outerCurve.edge_length(eTo) + outerCurve.edge_length(eFrom));

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
			const auto currentCoVolLength = 0.5 * (innerCurve.edge_length(eTo) + innerCurve.edge_length(eFrom));

			if (currentCoVolLength > maxCoVolLength) maxCoVolLength = currentCoVolLength;
			if (currentCoVolLength < minCoVolLength) minCoVolLength = currentCoVolLength;
		}
	}

	return { minCoVolLength, maxCoVolLength };
}

void CustomManifoldCurveEvolutionStrategy::StabilizeCustomGeometries(pmp::Scalar minLength, pmp::Scalar maxLength, pmp::Scalar stabilizationFactor)
{
	pmp::Scalar scalingFactor{ 1.0 };
	if (GetSettings().UseStabilizationViaScaling)
	{
		const pmp::Scalar expectedMeanCoVolLength = (1.0 - stabilizationFactor) * minLength + stabilizationFactor * maxLength;
		scalingFactor = pow(static_cast<pmp::Scalar>(GetSettings().TimeStep) / expectedMeanCoVolLength * 1.0, SCALE_FACTOR_POWER_1D);
		std::cout << "StabilizeCustomGeometries: Calculated scaling factor: " << scalingFactor << "\n";
	}
	GetScalingFactor() = scalingFactor;		

	// -----------------------------------------------------------------------------------------------
	// All geometric quantities need to be scaled by the scalingFactor to ensure numerical stability.
	// For spatial conveniences, geometries must also be centered, i.e.: translated by -origin.
	// Geometries that are already centered at (0, 0) are not translated.
	// -----------------------------------------------------------------------------------------------
	GetSettings().FieldSettings.FieldIsoLevel *= scalingFactor;
	GetFieldCellSize() *= scalingFactor;

	const pmp::mat3 transfMatrixGeomScale{
		scalingFactor, 0.0, 0.0,
		0.0, scalingFactor, 0.0,
		0.0, 0.0, 1.0
	};

	// extract origin and radius for m_EvolBox from the custom geometries
	pmp::Point2 origin{};
	pmp::Scalar radius = 1.0;
	if (GetOuterCurve())
	{
		const auto outerBounds = GetOuterCurve()->bounds();
		const auto outerBoundsSize = outerBounds.max() - outerBounds.min();
		radius = std::max(outerBoundsSize[0], outerBoundsSize[1]) * 0.5;
		origin = outerBounds.center();
	}
	else
	{
		if (GetDistanceField())
		{
			// non-empty target distance field
			const auto& dfBox = GetDistanceField()->Box();
			const auto dfBoxSize = dfBox.max() - dfBox.min();
			radius = std::max(std::max(dfBoxSize[0], dfBoxSize[1]), radius);
			origin = dfBox.center();
		}
		else
		{
			for (const auto& innerCurve : GetInnerCurves())
			{
				const auto innerBounds = innerCurve->bounds();
				const auto innerBoundsSize = innerBounds.max() - innerBounds.min();
				radius = std::max(std::max(innerBoundsSize[0], innerBoundsSize[1]) * (pmp::Scalar)0.8, radius);
				origin += innerBounds.center();
			}
			// mean bounds center position for all initial curves
			if (!GetInnerCurves().empty())
				origin /= static_cast<pmp::Scalar>(GetInnerCurves().size());
		}
	}

	const pmp::mat3 transfMatrixGeomMove{
		1.0, 0.0, -origin[0],
		0.0, 1.0, -origin[1],
		0.0, 0.0, 1.0
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	GetTransformToOriginal() = inverse(transfMatrixFull);
	GetScaleFieldToOriginal() = inverse(transfMatrixGeomScale);

	// Transform the geometries so that they're centered at (0, 0) and in proper scale
	if (GetOuterCurve())
	{
		(*GetOuterCurve()) *= transfMatrixFull;
	}
	for (auto& innerCurve : GetInnerCurves())
	{
		(*innerCurve) *= transfMatrixFull;
	}

	// test box for geometry validation
	const pmp::Scalar evolBoxFactor = 5.0 * scalingFactor;
	GetEvolBox() = pmp::BoundingBox2(
		pmp::Point2{ -radius, -radius } * evolBoxFactor,
		pmp::Point2{ radius, radius } * evolBoxFactor);

	if (GetOuterCurveDistanceField())
	{
		(*GetOuterCurveDistanceField()) *= transfMatrixFull;
		(*GetOuterCurveDistanceField()) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	if (GetOuterCurveDFNegNormalizedGradient())
	{
		(*GetOuterCurveDFNegNormalizedGradient()) *= transfMatrixFull;
		// values are supposed to be unit vectors regardless of scaling
	}
	for (const auto& innerCurveDF : GetInnerCurvesDistanceFields())
	{
		(*innerCurveDF) *= transfMatrixFull;
		(*innerCurveDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerCurveDFGradient : GetInnerCurvesDFNegNormalizedGradients())
	{
		(*innerCurveDFGradient) *= transfMatrixFull;
		// values are supposed to be unit vectors regardless of scaling
	}

	if (!GetDistanceField() || !GetDFNegNormalizedGradient())
		return; // nothing to scale

	(*GetDistanceField()) *= transfMatrixFull;
	(*GetDistanceField()) *= static_cast<double>(scalingFactor); // Scale the distance values
	(*GetDFNegNormalizedGradient()) *= transfMatrixFull;
	// values are supposed to be unit vectors regardless of scaling
}

//
// ======================================================================================
//                    The strategy for 2D Surfaces in 3D space
// ---------------------------------------------------------------------------------------
//

void ManifoldSurfaceEvolutionStrategy::InitLogger(const std::string& baseOutputFileName)
{
	if (!LogManifoldValues())
		return;

	m_Logger.Init(baseOutputFileName + "_log.json");
}

void ManifoldSurfaceEvolutionStrategy::InitNewTimeStepLog(unsigned int stepId)
{
	if (!LogManifoldValues())
		return;

	m_Logger.StartNewTimeStep(stepId); // Explicitly start a new logging step
}

void ManifoldSurfaceEvolutionStrategy::SaveLog(bool omitLastTimeStep)
{
	if (!LogManifoldValues())
		return;

	m_Logger.Save(omitLastTimeStep);
}

void ManifoldSurfaceEvolutionStrategy::Preprocess()
{
	const auto [minTargetSize, maxTargetSize, targetCenter] = ComputeAmbientFields();
	ConstructInitialManifolds(minTargetSize, maxTargetSize, targetCenter);

	GetFieldCellSize() = m_DistanceField ? m_DistanceField->CellSize() : minTargetSize / static_cast<pmp::Scalar>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
	ComputeVariableDistanceFields();

	if (GetSettings().UseStabilizationViaScaling)
	{
		StabilizeGeometries();
	}
	ComputeControlFunctionsLowerBounds();

	AssignRemeshingSettingsToEvolvingManifolds();
}

void ManifoldSurfaceEvolutionStrategy::Postprocess()
{
}

void CustomManifoldSurfaceEvolutionStrategy::Preprocess()
{
	if (!GetOuterSurface() && GetInnerSurfaces().empty())
		throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: There's nothing to evolve!\n");

	//if (!HasValidInnerOuterManifolds())
	//	throw std::invalid_argument("CustomManifoldSurfaceEvolutionStrategy::Preprocess: Invalid inner /outer manifold geometry! Not all custom inner surfaces are contained within the custom outer surface.\n");

	if (LogManifoldValues())
	{
		if (GetOuterSurface())
			GetLogger().AddManifold(GetOuterSurface().get());
		for (auto& innerSurface : GetInnerSurfaces())
		{
			if (!innerSurface)
				continue;

			GetLogger().AddManifold(innerSurface.get());
		}
	}

	if (NeedsFieldsCalculation())
	{
		GetEvolBox() = GetOuterSurface()->bounds();
		const auto sizeVec = (GetEvolBox().max() - GetEvolBox().min()) * GetSettings().FieldSettings.FieldExpansionFactor;
		GetEvolBox().expand(sizeVec[0], sizeVec[1], sizeVec[2]);

		std::tie(std::ignore, std::ignore, std::ignore) = ComputeAmbientFields();

		const auto minTargetSize = std::min({ sizeVec[0], sizeVec[1], sizeVec[2] });
		GetFieldCellSize() = GetDistanceField() ? GetDistanceField()->CellSize() : minTargetSize / static_cast<pmp::Scalar>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
		ComputeVariableDistanceFields();
	}

	//if (GetSettings().UseStabilizationViaScaling)
	{
		const auto [minArea, maxArea] = CalculateCoVolumeRange();
		std::cout << "Before stabilization: { minArea: " << minArea << ", maxArea: " << maxArea << "} ... vs ... timeStep: " << GetSettings().TimeStep << "\n";

		StabilizeCustomGeometries(minArea, maxArea);
		const auto [minAreaAfter, maxAreaAfter] = CalculateCoVolumeRange();
		std::cout << "After stabilization: { minAreaAfter: " << minAreaAfter << ", maxAreaAfter: " << maxAreaAfter << "} ... vs ... timeStep: " << GetSettings().TimeStep << "\n";
	}
	ComputeControlFunctionsLowerBounds();

	AssignRemeshingSettingsToEvolvingManifolds();
}

void ManifoldSurfaceEvolutionStrategy::PerformEvolutionStep(unsigned int stepId)
{
	if (LogManifoldValues())
		InitNewTimeStepLog(stepId);

	GetIntegrate()(stepId);
}

void ManifoldSurfaceEvolutionStrategy::Remesh()
{
	for (auto* surfaceToRemesh : m_RemeshTracker.GetManifoldsToRemesh())
	{
		pmp::Remeshing remesher(*surfaceToRemesh);
		remesher.adaptive_remeshing(m_RemeshingSettings[surfaceToRemesh]);

		//if (LogManifoldValues())
		//{
		//	m_Logger.ReserveBuffers(surfaceToRemesh);
		//}
	}
	m_RemeshTracker.Reset();
}

void ManifoldSurfaceEvolutionStrategy::ResizeRemeshingSettings(pmp::Scalar resizeFactor)
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
	if (GetSettings().ExportVariableScalarFieldsDimInfo && m_OuterSurfaceDistanceField)
	{
		auto exportedOuterSurfaceDF = *m_OuterSurfaceDistanceField;
		exportedOuterSurfaceDF *= m_TransformToOriginal;
		exportedOuterSurfaceDF /= static_cast<double>(GetScalingFactor());
		ExportToVTI(baseOutputFilename + "_OuterDF", exportedOuterSurfaceDF);
	}
	if (GetSettings().ExportVariableVectorFieldsDimInfo && m_OuterSurfaceDFNegNormalizedGradient)
	{
		auto exportedOuterSurfaceDFNegGrad = *m_OuterSurfaceDFNegNormalizedGradient;
		exportedOuterSurfaceDFNegGrad *= m_TransformToOriginal;
		ExportToVTK(baseOutputFilename + "_OuterDFNegGrad", exportedOuterSurfaceDFNegGrad);
	}
	if (GetSettings().ExportVariableScalarFieldsDimInfo)
	{
		for (size_t i = 0; i < m_InnerSurfacesDistanceFields.size(); ++i)
		{
			if (!m_InnerSurfacesDistanceFields[i])
				continue;

			auto exportedInnerSurfaceDF = *m_InnerSurfacesDistanceFields[i];
			exportedInnerSurfaceDF *= m_TransformToOriginal;
			exportedInnerSurfaceDF /= static_cast<double>(GetScalingFactor());
			ExportToVTI(baseOutputFilename + "_InnerDF" + std::to_string(i), exportedInnerSurfaceDF);
		}
	}
	if (GetSettings().ExportVariableVectorFieldsDimInfo)
	{
		for (size_t i = 0; i < m_InnerSurfacesDFNegNormalizedGradients.size(); ++i)
		{
			if (!m_InnerSurfacesDFNegNormalizedGradients[i])
				continue;

			auto exportedInnerSurfaceDFNegGrad = *m_InnerSurfacesDFNegNormalizedGradients[i];
			exportedInnerSurfaceDFNegGrad *= m_TransformToOriginal;
			ExportToVTK(baseOutputFilename + "_InnerDFNegGrad" + std::to_string(i), exportedInnerSurfaceDFNegGrad);
		}
	}
	// ----------------------------------------------------------------

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

		const auto tStep = GetSettings().TimeStep;

		pmp::Normals::compute_vertex_normals(*m_OuterSurface);
		auto vNormalsProp = m_OuterSurface->get_vertex_property<pmp::vec3>("v:normal");

		// prepare matrix & rhs for m_OuterSurface:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 6);  // Assuming an average of 6 entries per vertex

		for (const auto v : m_OuterSurface->vertices())
		{
			const auto vPosToUpdate = m_OuterSurface->position(v);

			InteractionDistanceInfo<pmp::dvec3> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec3>{vDistanceToTarget, vNegGradDistanceToTarget};

			double vMinDistanceToInner = DBL_MAX;
			for (unsigned int i = 0; i < m_InnerSurfacesDistanceFields.size(); ++i)
			{
				if (!m_InnerSurfacesDistanceFields[i] || !m_InnerSurfacesDFNegNormalizedGradients[i])
					continue;

				auto innerDfAtVPos = m_ScalarInterpolate(vPosToUpdate, *m_InnerSurfacesDistanceFields[i]);
				innerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToInner = m_VectorInterpolate(vPosToUpdate, *m_InnerSurfacesDFNegNormalizedGradients[i]);

				interaction << InteractionDistanceInfo<pmp::dvec3>{innerDfAtVPos, vNegGradDistanceToInner};

				if (innerDfAtVPos < vMinDistanceToInner)
					vMinDistanceToInner = innerDfAtVPos;
			}

			if (m_OuterSurface->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector3d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const double epsilonCtrlWeight =
				GetSettings().OuterManifoldEpsilon(static_cast<double>(interaction.Distance));

			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				(GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget), vNormal), -1.0, 1.0);
			const double advectionDistance =
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight =
				(m_DistanceField || !m_InnerSurfacesDistanceFields.empty()) ? GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal) : 0.0;

			if (LogOuterManifoldValues())
			{
				if (GetSettings().DiagSettings.LogOuterManifoldEpsilon)
					m_Logger.LogValue(m_OuterSurface.get(), "epsilonCtrlWeight", v.idx(), epsilonCtrlWeight);
				if (GetSettings().DiagSettings.LogOuterManifoldEta)
					m_Logger.LogValue(m_OuterSurface.get(), "etaCtrlWeight", v.idx(), etaCtrlWeight);
			}

			const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight);
			if (tanRedistWeight > 0.0)
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
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep: vertex " + std::to_string(i) 
					+ ": (" + std::to_string(newPos[0]) + ", " + std::to_string(newPos[1]) + ", " + std::to_string(newPos[2]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1]) + ", " + std::to_string(m_EvolBox.min()[2])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1]) + ", " + std::to_string(m_EvolBox.max()[2])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
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

		const auto tStep = GetSettings().TimeStep;

		pmp::Normals::compute_vertex_normals(*innerSurface);
		auto vNormalsProp = innerSurface->get_vertex_property<pmp::vec3>("v:normal");

		// prepare matrix & rhs for m_InnerSurface:
		std::vector<Eigen::Triplet<double>> tripletList;
		tripletList.reserve(static_cast<size_t>(NVertices) * 6);  // Assuming an average of 6 entries per vertex

		for (const auto v : innerSurface->vertices())
		{
			const auto vPosToUpdate = innerSurface->position(v);

			InteractionDistanceInfo<pmp::dvec3> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec3>{vDistanceToTarget, vNegGradDistanceToTarget};

			double outerDfAtVPos = DBL_MAX;
			if (m_OuterSurfaceDistanceField && m_OuterSurfaceDFNegNormalizedGradient)
			{
				outerDfAtVPos = m_ScalarInterpolate(vPosToUpdate, *m_OuterSurfaceDistanceField);
				outerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToOuter = m_VectorInterpolate(vPosToUpdate, *m_OuterSurfaceDFNegNormalizedGradient);

				interaction << InteractionDistanceInfo<pmp::dvec3>{outerDfAtVPos, vNegGradDistanceToOuter};
			}

			if (innerSurface->is_boundary(v))
			{
				// freeze boundary/feature vertices
				const Eigen::Vector3d vertexRhs = vPosToUpdate;
				sysRhs.row(v.idx()) = vertexRhs;
				tripletList.emplace_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
				continue;
			}

			const double epsilonCtrlWeight =
				GetSettings().InnerManifoldEpsilon(static_cast<double>(interaction.Distance));

			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget, vNormal), -1.0, 1.0);
			const double advectionDistance =
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight =
				(m_DistanceField || m_OuterSurfaceDistanceField) ? GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal) : 0.0;

			if (LogInnerManifoldValues())
			{
				if (GetSettings().DiagSettings.LogInnerManifoldsEpsilon)
					m_Logger.LogValue(innerSurface.get(), "epsilonCtrlWeight", v.idx(), epsilonCtrlWeight);
				if (GetSettings().DiagSettings.LogInnerManifoldsEta)
					m_Logger.LogValue(innerSurface.get(), "etaCtrlWeight", v.idx(), etaCtrlWeight);
			}

			const Eigen::Vector3d vertexRhs = vPosToUpdate + tStep * etaCtrlWeight * vNormal;
			sysRhs.row(v.idx()) = vertexRhs;
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight);
			if (tanRedistWeight > 0.0)
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
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::SemiImplicitIntegrationStep: innerSurface vertex " + std::to_string(i)
					+ ": (" + std::to_string(newPos[0]) + ", " + std::to_string(newPos[1]) + ", " + std::to_string(newPos[2]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1]) + ", " + std::to_string(m_EvolBox.min()[2])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1]) + ", " + std::to_string(m_EvolBox.max()[2])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
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

		pmp::Normals::compute_vertex_normals(*m_OuterSurface);
		auto vNormalsProp = m_OuterSurface->get_vertex_property<pmp::vec3>("v:normal");

		for (const auto v : m_OuterSurface->vertices())
		{
			const auto vPosToUpdate = m_OuterSurface->position(v);

			InteractionDistanceInfo<pmp::dvec3> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec3>{vDistanceToTarget, vNegGradDistanceToTarget};

			double vMinDistanceToInner = DBL_MAX;
			for (unsigned int i = 0; i < m_InnerSurfacesDistanceFields.size(); ++i)
			{
				if (!m_InnerSurfacesDistanceFields[i] || !m_InnerSurfacesDFNegNormalizedGradients[i])
					continue;

				auto innerDfAtVPos = m_ScalarInterpolate(vPosToUpdate, *m_InnerSurfacesDistanceFields[i]);
				innerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToInner = m_VectorInterpolate(vPosToUpdate, *m_InnerSurfacesDFNegNormalizedGradients[i]);

				interaction << InteractionDistanceInfo<pmp::dvec3>{innerDfAtVPos, vNegGradDistanceToInner};

				if (innerDfAtVPos < vMinDistanceToInner)
					vMinDistanceToInner = innerDfAtVPos;
			}

			if (m_OuterSurface->is_boundary(v))
				continue; // skip boundary vertices

			const double epsilonCtrlWeight =
				GetSettings().OuterManifoldEpsilon(static_cast<double>(interaction.Distance)) +
				GetSettings().OuterManifoldRepulsion(static_cast<double>(vMinDistanceToInner));

			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				(GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget), vNormal), -1.0, 1.0);
			const double advectionDistance =
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight =
				(m_DistanceField || !m_InnerSurfacesDistanceFields.empty()) ? GetSettings().OuterManifoldEta(advectionDistance, negGradDotNormal) : 0.0;

			if (LogOuterManifoldValues())
			{
				if (GetSettings().DiagSettings.LogOuterManifoldEpsilon)
					m_Logger.LogValue(m_OuterSurface.get(), "epsilonCtrlWeight", v.idx(), epsilonCtrlWeight);
				if (GetSettings().DiagSettings.LogOuterManifoldEta)
					m_Logger.LogValue(m_OuterSurface.get(), "etaCtrlWeight", v.idx(), etaCtrlWeight);
			}

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * m_ExplicitLaplacianFunction(*m_OuterSurface, v);

			// Tangential redistribution velocity
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight);
			pmp::vec3 tanVelocity(0.0, 0.0, 0.0);
			if (tanRedistWeight > 0.0)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*m_OuterSurface, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::ExplicitIntegrationStep: vertex " + std::to_string(v.idx())
					+ ": (" + std::to_string(updatedPosition[0]) + ", " + std::to_string(updatedPosition[1]) + ", " + std::to_string(updatedPosition[2]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1]) + ", " + std::to_string(m_EvolBox.min()[2])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1]) + ", " + std::to_string(m_EvolBox.max()[2])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
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

		pmp::Normals::compute_vertex_normals(*innerSurface);
		auto vNormalsProp = innerSurface->get_vertex_property<pmp::vec3>("v:normal");

		for (const auto v : innerSurface->vertices())
		{
			const auto vPosToUpdate = innerSurface->position(v);

			InteractionDistanceInfo<pmp::dvec3> interaction{};

			double vDistanceToTarget = m_DistanceField ? m_ScalarInterpolate(vPosToUpdate, *m_DistanceField) : DBL_MAX;
			vDistanceToTarget -= GetSettings().FieldSettings.FieldIsoLevel;
			const auto vNegGradDistanceToTarget = m_DFNegNormalizedGradient ? m_VectorInterpolate(vPosToUpdate, *m_DFNegNormalizedGradient) : pmp::dvec3(0, 0, 0);

			interaction << InteractionDistanceInfo<pmp::dvec3>{vDistanceToTarget, vNegGradDistanceToTarget};

			double outerDfAtVPos = DBL_MAX;
			if (m_OuterSurfaceDistanceField && m_OuterSurfaceDFNegNormalizedGradient)
			{
				outerDfAtVPos = m_ScalarInterpolate(vPosToUpdate, *m_OuterSurfaceDistanceField);
				outerDfAtVPos -= GetSettings().FieldSettings.FieldIsoLevel;
				const auto vNegGradDistanceToOuter = m_VectorInterpolate(vPosToUpdate, *m_OuterSurfaceDFNegNormalizedGradient);

				interaction << InteractionDistanceInfo<pmp::dvec3>{outerDfAtVPos, vNegGradDistanceToOuter};
			}

			if (innerSurface->is_boundary(v))
				continue; // skip boundary vertices

			const double epsilonCtrlWeight =
				GetSettings().InnerManifoldEpsilon(static_cast<double>(interaction.Distance)) +
				GetSettings().InnerManifoldRepulsion(static_cast<double>(outerDfAtVPos));

			const auto vNormal = static_cast<pmp::vec3>(vNormalsProp[v]); // vertex unit normal

			const auto negGradDotNormal = std::clamp(pmp::ddot(
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.NegGradient : vNegGradDistanceToTarget, vNormal), -1.0, 1.0);
			const double advectionDistance =
				GetSettings().AdvectionInteractWithOtherManifolds ? interaction.Distance : vDistanceToTarget;
			const double etaCtrlWeight =
				(m_DistanceField || m_OuterSurfaceDistanceField) ? GetSettings().InnerManifoldEta(advectionDistance, negGradDotNormal) : 0.0;

			if (LogInnerManifoldValues())
			{
				if (GetSettings().DiagSettings.LogInnerManifoldsEpsilon)
					m_Logger.LogValue(innerSurface.get(), "epsilonCtrlWeight", v.idx(), epsilonCtrlWeight);
				if (GetSettings().DiagSettings.LogInnerManifoldsEta)
					m_Logger.LogValue(innerSurface.get(), "etaCtrlWeight", v.idx(), etaCtrlWeight);
			}

			// Laplacian term (already weighted by epsilon and area)
			const auto laplacianTerm = epsilonCtrlWeight * m_ExplicitLaplacianFunction(*innerSurface, v);

			// Tangential redistribution velocity
			const auto tanRedistWeight = static_cast<double>(GetSettings().TangentialVelocityWeight) * std::abs(epsilonCtrlWeight);
			pmp::vec3 tanVelocity(0.0, 0.0, 0.0);
			if (tanRedistWeight > 0.0)
			{
				tanVelocity = ComputeTangentialUpdateVelocityAtVertex(*innerSurface, v, vNormal, tanRedistWeight);
			}

			// Update the vertex position explicitly
			auto updatedPosition = vPosToUpdate + tStep * (laplacianTerm + etaCtrlWeight * vNormal + tanVelocity);

			// Check if the updated position is within bounds
			if (!m_EvolBox.Contains(updatedPosition))
			{
				const std::string msg = "\nManifoldSurfaceEvolutionStrategy::ExplicitIntegrationStep: innerSurface vertex " + std::to_string(v.idx()) 
					+ ": (" + std::to_string(updatedPosition[0]) + ", " + std::to_string(updatedPosition[1]) + ", " + std::to_string(updatedPosition[2]) + ") outside m_EvolBox: {"
					+ "min_=(" + std::to_string(m_EvolBox.min()[0]) + ", " + std::to_string(m_EvolBox.min()[1]) + ", " + std::to_string(m_EvolBox.min()[2])
					+ "), max_=(" + std::to_string(m_EvolBox.max()[0]) + ", " + std::to_string(m_EvolBox.max()[1]) + ", " + std::to_string(m_EvolBox.max()[2])
					+ ")} for time step id: " + std::to_string(step) + "!\n";
				std::cerr << msg;
				throw std::runtime_error(msg);
			}

			innerSurface->position(v) = updatedPosition;
		}
	}
}


std::tuple<pmp::Scalar, pmp::Scalar, pmp::Point> ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields()
{
	if (!m_TargetPointCloud)
	{
		std::cerr << "ManifoldSurfaceEvolutionStrategy::ComputeAmbientFields: No m_TargetPointCloud found! Initializing empty fields: m_DistanceField and m_DFNegNormalizedGradient.\n";
		return { FLT_MAX, FLT_MAX, pmp::Point(0, 0, 0)};
	}

	const pmp::BoundingBox ptCloudBBox(*m_TargetPointCloud);
	const auto ptCloudBBoxSize = ptCloudBBox.max() - ptCloudBBox.min();
	const pmp::Scalar minSize = std::min({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const pmp::Scalar maxSize = std::max({ ptCloudBBoxSize[0], ptCloudBBoxSize[1], ptCloudBBoxSize[2] });
	const auto cellSize = minSize / static_cast<pmp::Scalar>(GetSettings().FieldSettings.NVoxelsPerMinDimension);
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
	if (!NeedsVariableFieldsCalculation())
	{
		// there's no possibility of interaction between the outer and the inner manifolds
		return;
	}

	// clear fields
	m_OuterSurfaceDistanceField.reset();
	m_OuterSurfaceDFNegNormalizedGradient.reset();
	m_InnerSurfacesDistanceFields.clear();
	m_InnerSurfacesDFNegNormalizedGradients.clear();

	SDF::DistanceFieldSettings surfaceDFSettings{
		GetFieldCellSize(),
		GetSettings().FieldSettings.FieldExpansionFactor,
		DBL_MAX,
		SDF::KDTreeSplitType::Center,
		SDF::SignComputation::None,
		SDF::BlurPostprocessingType::None,
		SDF::PreprocessingType::Octree
	};
	const std::optional<pmp::BoundingBox> outerFieldBox = m_DistanceField ? std::make_optional(m_DistanceField->Box()) : std::nullopt;
	const Geometry::PMPSurfaceMeshAdapter outerSurfaceAdapter(std::make_shared<pmp::SurfaceMesh>(*m_OuterSurface));
	m_OuterSurfaceDistanceField = std::make_shared<Geometry::ScalarGrid>(
		SDF::DistanceFieldGenerator::Generate(outerSurfaceAdapter, surfaceDFSettings, outerFieldBox));
	RepairScalarGrid(*m_OuterSurfaceDistanceField);
	m_OuterSurfaceDFNegNormalizedGradient = std::make_shared<Geometry::VectorGrid>(ComputeNormalizedNegativeGradient(*m_OuterSurfaceDistanceField));

	for (const auto& innerSurface : m_InnerSurfaces)
	{
		const Geometry::PMPSurfaceMeshAdapter innerSurfaceAdapter(std::make_shared<pmp::SurfaceMesh>(*innerSurface));
		m_InnerSurfacesDistanceFields.emplace_back(std::make_shared<Geometry::ScalarGrid>(
			SDF::DistanceFieldGenerator::Generate(innerSurfaceAdapter, surfaceDFSettings, m_OuterSurfaceDistanceField->Box())));
		m_InnerSurfacesDFNegNormalizedGradients.emplace_back(
			std::make_shared<Geometry::VectorGrid>(ComputeNormalizedNegativeGradient(*m_InnerSurfacesDistanceFields.back())));
	}
}

void ManifoldSurfaceEvolutionStrategy::ConstructInitialManifolds(pmp::Scalar minTargetSize, pmp::Scalar maxTargetSize, const pmp::Point& targetBoundsCenter)
{
	if (!GetSettings().UseInnerManifolds && !GetSettings().UseOuterManifolds)
	{
		throw std::invalid_argument("ManifoldSurfaceEvolutionStrategy::ConstructInitialManifolds: Current setting is: UseInnerManifolds == false && UseOuterManifolds == false. This means there's nothing to evolve!\n");
	}

	const pmp::Scalar outerSphereRadius = 0.5 * SPHERE_RADIUS_FACTOR *
		(minTargetSize + (0.5 + GetSettings().FieldSettings.FieldExpansionFactor) * maxTargetSize);

	if (GetSettings().UseOuterManifolds)
	{
		Geometry::IcoSphereBuilder icoBuilder({ GetSettings().LevelOfDetail, outerSphereRadius });
		icoBuilder.BuildBaseData();
		icoBuilder.BuildPMPSurfaceMesh();
		m_OuterSurface = std::make_shared<pmp::SurfaceMesh>(icoBuilder.GetPMPSurfaceMeshResult());
		// DISCLAIMER: Since we want to evolve manifolds centered at the origin, we will not move the outer manifold into its "true" position.
		// The "true" position will be stored in m_InitialSphereSettings:
		m_InitialSphereSettings[m_OuterSurface.get()] = Sphere3D{ targetBoundsCenter, outerSphereRadius };

		if (LogManifoldValues())
			m_Logger.AddManifold(m_OuterSurface.get());
	}

	if (!GetSettings().UseInnerManifolds || !m_TargetPointCloud || !m_DistanceField)
		return;

	throw std::invalid_argument(
		"ManifoldSurfaceEvolutionStrategy::ConstructInitialManifolds: Automatic calculation of inscribed manifolds is not reliable yet! Please use CustomManifoldSurfaceEvolutionStrategy instead.");

	// TODO: This is hardcoded, implement 3D version of ParticleSwarmDistanceFieldInscribedSphereCalculator
	//const InscribedSphereInputData calcData{
	//	*m_TargetPointCloud,
	//	std::make_shared<Geometry::ScalarGrid>(*m_DistanceField) // clone
	//};	
	//ParticleSwarmDistanceFieldInscribedSphereCalculator inscribedSphereCalculator;
	//const auto spheres = inscribedSphereCalculator.Calculate(calcData);
	//const std::vector spheres = {Sphere3D{pmp::Point(0, 0, 0), 0.7}};

	//m_InnerSurfaces.reserve(spheres.size());

	//for (const auto& sphere : spheres)
	//{
	//	// keep the same vertex density for inner spheres
	//	const auto innerSubdiv = static_cast<unsigned int>(static_cast<pmp::Scalar>(GetSettings().LevelOfDetail) * (sphere.Radius * SPHERE_RADIUS_FACTOR) / outerSphereRadius);
	//	Geometry::IcoSphereBuilder innerIcoBuilder({ 
	//		innerSubdiv,
	//		sphere.Radius * SPHERE_RADIUS_FACTOR
	//	});
	//	innerIcoBuilder.BuildBaseData();
	//	innerIcoBuilder.BuildPMPSurfaceMesh();
	//	auto mesh = innerIcoBuilder.GetPMPSurfaceMeshResult();
	//	if (sphere.Center != pmp::Point(0, 0, 0))
	//	{
	//		const auto translationMatrix = translation_matrix(sphere.Center);
	//		mesh *= translationMatrix;
	//	}		
	//	m_InnerSurfaces.push_back(std::make_shared<pmp::SurfaceMesh>(mesh));
	//	m_InitialSphereSettings[m_InnerSurfaces.back().get()] = sphere;
	//}
}

/// \brief The power of the stabilizing scale factor.
constexpr pmp::Scalar SCALE_FACTOR_POWER_2D = 1.0 / 2.0;
/// \brief the reciprocal value of how many times the surface area element shrinks during evolution.
constexpr pmp::Scalar INV_SHRINK_FACTOR_2D = 5.0;

void ManifoldSurfaceEvolutionStrategy::StabilizeGeometries(pmp::Scalar stabilizationFactor)
{
	const pmp::Scalar radius = stabilizationFactor * m_InitialSphereSettings.MinRadius() + (1.0 - stabilizationFactor) * m_InitialSphereSettings.MaxRadius();
	if (radius <= 0.0)
	{
		throw std::invalid_argument("ManifoldSurfaceEvolutionStrategy::StabilizeGeometries: m_InitialSphereSettings empty!\n");
	}
	const unsigned int expectedVertexCount = (N_ICO_EDGES_0 * static_cast<unsigned int>(pow(4, GetSettings().LevelOfDetail) - 1) + 3 * N_ICO_VERTS_0) / 3;
	const pmp::Scalar expectedMeanCoVolArea = (4.0 * static_cast<pmp::Scalar>(M_PI) * radius * radius / static_cast<pmp::Scalar>(expectedVertexCount));
	const auto scalingFactor = pow(static_cast<pmp::Scalar>(GetSettings().TimeStep) / expectedMeanCoVolArea * INV_SHRINK_FACTOR_2D, SCALE_FACTOR_POWER_2D);
	GetScalingFactor() = scalingFactor;
	
	// -----------------------------------------------------------------------------------------------
	// All geometric quantities need to be scaled by the scalingFactor to ensure numerical stability.
	// For spatial conveniences, geometries must also be centered, i.e.: translated by -origin.
	// Geometries that are already centered at (0, 0, 0) are not translated.
	// -----------------------------------------------------------------------------------------------
	GetSettings().FieldSettings.FieldIsoLevel *= scalingFactor;
	GetFieldCellSize() *= scalingFactor;

	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0, 0.0, 0.0,
		0.0, scalingFactor, 0.0, 0.0,
		0.0, 0.0, scalingFactor, 0.0,
		0.0, 0.0, 0.0, 1.0
	};
	const pmp::Point origin = m_OuterSurface ? m_InitialSphereSettings[m_OuterSurface.get()].Center :
		(!m_InnerSurfaces.empty() ? m_InitialSphereSettings[m_InnerSurfaces[0].get()].Center : pmp::Point{0, 0, 0});
	const pmp::mat4 transfMatrixGeomMove{
		1.0, 0.0, 0.0, -origin[0],
		0.0, 1.0, 0.0, -origin[1],
		0.0, 0.0, 1.0, -origin[2],
		0.0, 0.0, 0.0, 1.0
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	m_TransformToOriginal = inverse(transfMatrixFull);
	m_ScaleFieldToOriginal = inverse(transfMatrixGeomScale);

	// transform geometries
	if (m_OuterSurface)
	{
		(*m_OuterSurface) *= transfMatrixGeomScale;		
	}
	for (const auto& innerSurface : m_InnerSurfaces)
	{
		(*innerSurface) *= transfMatrixFull;
	}

	// test box for geometry validation
	const pmp::Scalar evolBoxFactor = 5.0 * scalingFactor;
	m_EvolBox = pmp::BoundingBox(
		pmp::Point{ -radius, -radius, -radius } * evolBoxFactor,
		pmp::Point{ radius, radius, radius } * evolBoxFactor);

	if (m_OuterSurfaceDistanceField)
	{
		(*m_OuterSurfaceDistanceField) *= transfMatrixGeomScale;
		(*m_OuterSurfaceDistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	if (m_OuterSurfaceDFNegNormalizedGradient)
	{
		(*m_OuterSurfaceDFNegNormalizedGradient) *= transfMatrixGeomScale;
		// values are supposed to be unit vectors regardless of scaling
	}
	for (const auto& innerSurfaceDF : m_InnerSurfacesDistanceFields)
	{
		(*innerSurfaceDF) *= transfMatrixFull;
		(*innerSurfaceDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerSurfaceDFGradient : m_InnerSurfacesDFNegNormalizedGradients)
	{
		(*innerSurfaceDFGradient) *= transfMatrixFull;
		// values are supposed to be unit vectors regardless of scaling
	}

	if (!m_DistanceField || !m_DFNegNormalizedGradient)
		return; // nothing to scale

	(*m_DistanceField) *= transfMatrixFull;
	(*m_DistanceField) *= static_cast<double>(scalingFactor); // scale also the distance values.
	(*m_DFNegNormalizedGradient) *= transfMatrixFull;
	// values are supposed to be unit vectors regardless of scaling
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

bool ManifoldSurfaceEvolutionStrategy::NeedsVariableFieldsCalculation()
{
	return m_OuterSurface && !m_InnerSurfaces.empty();
}

bool ManifoldSurfaceEvolutionStrategy::NeedsFieldsCalculation()
{
	return NeedsVariableFieldsCalculation() || m_TargetPointCloud;
}

void ManifoldSurfaceEvolutionStrategy::ComputeControlFunctionsLowerBounds()
{
	// find reasonable bounds for distance
	double maxDistance{ 0.0 };
	if (m_DistanceField)
	{
		maxDistance = *std::ranges::max_element(m_DistanceField->Values());
	}
	else if (!m_DistanceField && m_OuterSurfaceDistanceField)
	{
		maxDistance = *std::ranges::max_element(m_OuterSurfaceDistanceField->Values());
	}
	else if (!m_DistanceField && !m_OuterSurfaceDistanceField && !m_InnerSurfacesDistanceFields.empty())
	{
		maxDistance = *std::ranges::max_element(m_InnerSurfacesDistanceFields.front()->Values());
	}

	// normalize the max distance to expansion factor
	const auto expFactor = GetSettings().FieldSettings.FieldExpansionFactor;
	if (expFactor < FLT_EPSILON)
	{
		throw std::invalid_argument("ManifoldSurfaceEvolutionStrategy::ComputeControlFunctionsLowerBounds: expFactor < FLT_EPSILON!\n");
	}
	maxDistance /= (0.5 + expFactor);

	// apply the limit
	GetSettings().OuterManifoldEpsilon.ComputeLimit(maxDistance);
	GetSettings().OuterManifoldEta.ComputeLimit(maxDistance);
	GetSettings().InnerManifoldEpsilon.ComputeLimit(maxDistance);
	GetSettings().InnerManifoldEta.ComputeLimit(maxDistance);
}

void CustomManifoldSurfaceEvolutionStrategy::AssignRemeshingSettingsToEvolvingManifolds()
{
	if (GetOuterSurface())
	{
		GetRemeshingSettings()[GetOuterSurface().get()] = 
			CollectRemeshingSettingsForManifold_EXPERIMENTAL(GetOuterSurface(), GetSettings().RemeshingSettings);
	}

	for (const auto& innerSurface : GetInnerSurfaces())
	{
		GetRemeshingSettings()[innerSurface.get()] = 
			CollectRemeshingSettingsForManifold_EXPERIMENTAL(innerSurface, GetSettings().RemeshingSettings);
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

std::pair<pmp::Scalar, pmp::Scalar> CustomManifoldSurfaceEvolutionStrategy::CalculateCoVolumeRange() const
{
	pmp::Scalar minCoVolArea = FLT_MAX;
	pmp::Scalar maxCoVolArea = -FLT_MAX;

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

void CustomManifoldSurfaceEvolutionStrategy::StabilizeCustomGeometries(pmp::Scalar minArea, pmp::Scalar maxArea, pmp::Scalar stabilizationFactor)
{
	pmp::Scalar scalingFactor{ 1.0 };
	if (GetSettings().UseStabilizationViaScaling)
	{
		const pmp::Scalar expectedMeanCoVolLength = (1.0 - stabilizationFactor) * minArea + stabilizationFactor * maxArea;
		scalingFactor = pow(static_cast<pmp::Scalar>(GetSettings().TimeStep) / expectedMeanCoVolLength * INV_SHRINK_FACTOR_2D, SCALE_FACTOR_POWER_2D);
		std::cout << "StabilizeCustomGeometries: Calculated scaling factor: " << scalingFactor << "\n";
	}
	GetScalingFactor() = scalingFactor;		

	// -----------------------------------------------------------------------------------------------
	// All geometric quantities need to be scaled by the scalingFactor to ensure numerical stability.
	// For spatial conveniences, geometries must also be centered, i.e.: translated by -origin.
	// Geometries that are already centered at (0, 0, 0) are not translated.
	// -----------------------------------------------------------------------------------------------
	GetSettings().FieldSettings.FieldIsoLevel *= scalingFactor;
	GetFieldCellSize() *= scalingFactor;

	const pmp::mat4 transfMatrixGeomScale{
		scalingFactor, 0.0, 0.0, 0.0,
		0.0, scalingFactor, 0.0, 0.0,
		0.0, 0.0, scalingFactor, 0.0,
		0.0, 0.0, 0.0, 1.0
	};

	// extract origin and radius for m_EvolBox from the custom geometries
	pmp::Point origin{};
	pmp::Scalar radius = 1.0;
	if (GetOuterSurface())
	{
		const auto outerBounds = GetOuterSurface()->bounds();
		const auto outerBoundsSize = outerBounds.max() - outerBounds.min();
		radius = std::max({ outerBoundsSize[0], outerBoundsSize[1], outerBoundsSize[2] }) * 0.5;
		origin = outerBounds.center();
	}
	else
	{
		for (const auto& innerSurface : GetInnerSurfaces())
		{
			const auto innerBounds = innerSurface->bounds();
			const auto innerBoundsSize = innerBounds.max() - innerBounds.min();
			radius = std::max(std::max({ innerBoundsSize[0], innerBoundsSize[1], innerBoundsSize[2] }) * (pmp::Scalar)0.5, radius);
			origin += innerBounds.center();
		}
		// mean bounds center position for all initial curves
		if (!GetInnerSurfaces().empty())
			origin /= static_cast<pmp::Scalar>(GetInnerSurfaces().size());
	}
	const pmp::mat4 transfMatrixGeomMove{
		1.0, 0.0, 0.0, -origin[0],
		0.0, 1.0, 0.0, -origin[1],
		0.0, 0.0, 1.0, -origin[2],
		0.0, 0.0, 0.0, 1.0
	};
	const auto transfMatrixFull = transfMatrixGeomScale * transfMatrixGeomMove;
	GetTransformToOriginal() = inverse(transfMatrixFull);
	GetScaleFieldToOriginal() = inverse(transfMatrixGeomScale);

	// Transform the geometries
	if (GetOuterSurface())
	{
		(*GetOuterSurface()) *= transfMatrixFull;
	}
	for (auto& innerSurface : GetInnerSurfaces())
	{
		(*innerSurface) *= transfMatrixFull;
	}

	// test box for geometry validation
	const pmp::Scalar evolBoxFactor = 5.0 * scalingFactor;
	GetEvolBox() = pmp::BoundingBox(
		pmp::Point{ -radius, -radius, -radius } * evolBoxFactor,
		pmp::Point{ radius, radius, radius } * evolBoxFactor);

	if (GetOuterSurfraceDistanceField())
	{
		(*GetOuterSurfraceDistanceField()) *= transfMatrixFull; //transfMatrixGeomScale; //transfMatrixFull;
		(*GetOuterSurfraceDistanceField()) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	if (GetOuterSurfaceDFNegNormalizedGradient())
	{
		(*GetOuterSurfaceDFNegNormalizedGradient()) *= transfMatrixFull;
		// values are supposed to be unit vectors regardless of scaling
	}
	for (const auto& innerSurfaceDF : GetInnerSurfacesDistanceFields())
	{
		(*innerSurfaceDF) *= transfMatrixFull; //transfMatrixGeomScale; //transfMatrixFull;
		(*innerSurfaceDF) *= static_cast<double>(scalingFactor); // scale also the distance values.
	}
	for (const auto& innerSurfaceDFGradient : GetInnerSurfacesDFNegNormalizedGradients())
	{
		(*innerSurfaceDFGradient) *= transfMatrixFull; //transfMatrixGeomScale; //transfMatrixFull;
		// values are supposed to be unit vectors regardless of scaling
	}

	if (!GetDistanceField() || !GetDFNegNormalizedGradient())
		return; // nothing to scale

	(*GetDistanceField()) *= transfMatrixFull;
	(*GetDistanceField()) *= static_cast<double>(scalingFactor); // Scale the distance values
	(*GetDFNegNormalizedGradient()) *= transfMatrixFull;
	// values are supposed to be unit vectors regardless of scaling
}

// ================================================

bool ManifoldEvolutionStrategy::LogOuterManifoldValues() const
{
	return m_Settings.DiagSettings.LogOuterManifoldEpsilon ||
		m_Settings.DiagSettings.LogOuterManifoldEta ||
		m_Settings.DiagSettings.LogOuterManifoldErrors ||
		m_Settings.DiagSettings.LogOuterManifoldXErrors ||
		m_Settings.DiagSettings.LogInnerManifoldsXErrors;
}

bool ManifoldEvolutionStrategy::LogInnerManifoldValues() const
{
	return m_Settings.DiagSettings.LogInnerManifoldsEpsilon ||
		m_Settings.DiagSettings.LogInnerManifoldsEta ||
		m_Settings.DiagSettings.LogInnerManifoldsErrors ||
		m_Settings.DiagSettings.LogOuterManifoldYErrors ||
		m_Settings.DiagSettings.LogInnerManifoldsYErrors;
}

// ================================================