#include "pmp/algorithms/EdgeKdTree.h"

#include "ArcLengthCalculator.h"

using namespace pmp;

void pmp::EvolvingArcLengthCalculator::RecordPrevRefEdgePositions()
{
	if (!RefEdgeValid(true))
	{
		// nothing to do
		return;
	}

	const auto [v0, v1] = m_Curve.vertices(m_RefEdge.Edge);
	m_PrevRefEdge = std::make_shared<PrevRefEdgeRecord>(PrevRefEdgeRecord{ m_Curve.position(v0), m_Curve.position(v1), m_RefEdge.Param });
}

void EvolvingArcLengthCalculator::UpdateRefEdge()
{
	if (RefEdgeValid(true) || !m_PrevRefEdge)
	{
		// nothing to do
		return;
	}

	// project the original edge to the new remeshed curve
	const auto refPt = m_PrevRefEdge->StartPt * m_PrevRefEdge->Param + m_PrevRefEdge->EndPt * (1.0f - m_PrevRefEdge->Param);
	const auto kdTree = std::make_unique<EdgeKdTree>(std::make_shared<ManifoldCurve2D>(m_Curve), 0);
	const auto nearest = kdTree->nearest(refPt);
	m_RefEdge.Edge = nearest.edge;
	const auto newRefPt = nearest.nearest;
	const auto totalParamLength = m_Curve.edge_length(m_RefEdge.Edge);

	if (totalParamLength < FLT_EPSILON)
	{
		std::cerr << "\nEvolvingArcLengthCalculator::UpdateRefEdge: totalParamLength < FLT_EPSILON (this shouldn't happen)! \n";
		m_RefEdge.Param = 1.0f;
		m_PrevRefEdge = nullptr;
		return;
	}

	const auto v0New = m_Curve.from_vertex(m_RefEdge.Edge);
	const auto newParamDist = norm(newRefPt - m_Curve.position(v0New));
	if (newParamDist > totalParamLength)
	{
		std::cerr << "\nEvolvingArcLengthCalculator::UpdateRefEdge: newParamDist > totalParamLength (this shouldn't happen)! \n";
		m_RefEdge.Param = 1.0f;
		m_PrevRefEdge = nullptr;
		return;
	}
	m_RefEdge.Param = newParamDist / totalParamLength;

	m_PrevRefEdge = nullptr;
}

std::vector<Scalar> EvolvingArcLengthCalculator::CalculateArcLengths()
{
	if (!RefEdgeValid())
	{
		// nothing to do
		return {};
	}

	std::vector<Scalar> result(m_Curve.n_vertices(), -1.0f);

	Vertex vNext = m_Curve.to_vertex(m_RefEdge.Edge);
	Scalar cumulativeLength = m_Curve.edge_length(m_RefEdge.Edge) * (1.0f - m_RefEdge.Param);
	result[vNext.idx()] = cumulativeLength;
	Edge eCurrent{ m_Curve.edge_from(vNext) };

	unsigned int safetyCounter = 0;
	while (eCurrent != m_RefEdge.Edge || safetyCounter >= m_Curve.n_edges())
	{
		vNext = m_Curve.to_vertex(eCurrent);
		cumulativeLength += m_Curve.edge_length(eCurrent);
		result[vNext.idx()] = cumulativeLength;
		eCurrent = m_Curve.edge_from(vNext);
		++safetyCounter;
	}

	return result;
}

bool EvolvingArcLengthCalculator::RefEdgeValid(const bool& checkPosition)
{
	if (!m_RefEdge.Edge.is_valid())
	{
		// invalid indices
		return false;
	}

	if (m_Curve.is_deleted(m_RefEdge.Edge))
	{
		// indices referring to deleted vertices
		return false;
	}

	if (!checkPosition)
		return true;

	// position check for interpolation validity. 
	// By the triangle inequality, the Position must lie on the edge between PrevVertex and NextVertex
	const auto totalAffectedDistance = m_Curve.edge_length(m_RefEdge.Edge);
	if (totalAffectedDistance < FLT_EPSILON)
	{
		std::cerr << "\nEvolvingArcLengthCalculator::RefPointValid: totalAffectedDistance < FLT_EPSILON (this shouldn't happen)!\n";
		return false;
	}
	return m_RefEdge.Param >= 0.0f && m_RefEdge.Param <= totalAffectedDistance;
}
