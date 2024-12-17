#pragma once

#include "pmp/ManifoldCurve2D.h"

namespace pmp
{
	class EvolvingArcLengthCalculator
	{
	public:

		explicit EvolvingArcLengthCalculator(ManifoldCurve2D& curve)
			: m_Curve(curve)
		{
			const auto vFirst = Vertex(0);
			m_RefEdge = RefEdge{ m_Curve.edge_to(vFirst), 1.0f };
		}

		void RecordPrevRefEdgePositions();

		void UpdateRefEdge();

		[[nodiscard]] std::vector<Scalar> CalculateArcLengths();

	protected:

		[[nodiscard]] Edge GetRefEdgeId()
		{
			return m_RefEdge.Edge;
		}

	private:

		[[nodiscard]] bool RefEdgeValid(const bool& checkPosition = false);

		struct RefEdge
		{
			Edge Edge;
			Scalar Param;
		};

		struct PrevRefEdgeRecord
		{
			Point2 StartPt;
			Point2 EndPt;
			Scalar Param;
		};

		RefEdge m_RefEdge;
		ManifoldCurve2D& m_Curve;

		std::shared_ptr<PrevRefEdgeRecord> m_PrevRefEdge{nullptr};
	};

} // namespace pmp