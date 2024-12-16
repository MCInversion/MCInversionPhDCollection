#pragma once

#include "pmp/ManifoldCurve2D.h"

namespace pmp
{
	class EvolvingArcLengthCalculator
	{
	public:

		EvolvingArcLengthCalculator(ManifoldCurve2D& curve)
			: m_Curve(curve)
		{
			const auto vFirst = Vertex(0);
			m_RefEdge = RefEdge{ m_Curve.edge_to(vFirst), 1.0f };
		}

		void UpdateRefEdge();

		[[nodiscard]] std::vector<Scalar> CalculateArcLengths();

	private:

		[[nodiscard]] bool RefEdgeValid(const bool& checkPosition = false);

		struct RefEdge
		{
			Edge Edge;
			Scalar Param;
		};

		RefEdge m_RefEdge;
		ManifoldCurve2D& m_Curve;
	};

} // namespace pmp