#pragma once

#include "pmp/Types.h"

namespace Geometry
{
	//! \brief A recursive builder for 2D convex hull using the q-hull algorithm.
	class QHull2D
	{
	public:
		/// \brief Construct with input points.
		explicit QHull2D(const std::vector<pmp::Point2>& points)
			: m_SourcePoints(points)
		{
		}

		/// \brief Main builder functionality
		[[nodiscard]] std::vector<pmp::Point2> Build();

	private:

		void BuildRecurse(const std::pair<pmp::Point2, pmp::Point2>& edge);

		const std::vector<pmp::Point2>& m_SourcePoints;
		std::vector<pmp::Point2> m_ResultPoints;
	};

} // namespace Geometry