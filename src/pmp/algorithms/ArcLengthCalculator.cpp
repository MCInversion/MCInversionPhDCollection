#include "ArcLengthCalculator.h"

using namespace pmp;

void EvolvingArcLengthCalculator::UpdateRefEdge()
{
	if (RefEdgeValid(true))
	{
		// nothing to do
		return;
	}

	// TODO: no idea
}

std::vector<Scalar> EvolvingArcLengthCalculator::CalculateArcLengths()
{
	//if (!RefEdgeValid())
	//{
	//	// nothing to do
	//	return {};
	//}

	std::vector<Scalar> result(m_Curve.n_vertices(), -1.0f);

	Vertex vNext = m_Curve.to_vertex(m_RefEdge.Edge);
	Scalar cumulativeLength = m_Curve.edge_length(m_RefEdge.Edge) * (1.0f - m_RefEdge.Param);
	result[vNext.idx()] = cumulativeLength;
	Edge eCurrent{ m_Curve.edge_from(vNext) };

	unsigned int safetyCounter = 0;
	while (eCurrent != m_RefEdge.Edge || safetyCounter >= m_Curve.n_edges())
	{
		cumulativeLength = m_Curve.edge_length(eCurrent);
		result[vNext.idx()] = cumulativeLength;

		vNext = m_Curve.to_vertex(eCurrent);
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

	if (!m_Curve.is_deleted(m_RefEdge.Edge))
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
		//std::cerr << "EvolvingArcLengthCalculator::RefPointValid: totalAffectedDistance < FLT_EPSILON!\n";
		return false;
	}
	return m_RefEdge.Param >= 0.0f && m_RefEdge.Param <= totalAffectedDistance;
}
