#include "GeometryUtil.h"

#include "QHull2D.h"

#include <ranges>

namespace
{
	[[nodiscard]] pmp::Point GetFarthestPointFromLine2D(const std::pair<pmp::Point2, pmp::Point2>& line)
	{

	}


} // anonymous namespace

namespace Geometry
{
	std::vector<pmp::Point2> QHull2D::Build()
	{
		std::vector<pmp::Point2> resultPts;

		const auto [minPtIt, maxPtIt] = std::ranges::minmax_element(m_SourcePoints, [](const auto& a, const auto& b) {
			return (a[0] < b[0]);
		});
		m_ResultPoints.push_back(*minPtIt);
		m_ResultPoints.push_back(*maxPtIt);

		BuildRecurse({ *minPtIt, *maxPtIt });

		return std::vector<pmp::Point2>();
	}

	void QHull2D::BuildRecurse(const std::pair<pmp::Point2, pmp::Point2>& edge)
	{
		if (m_SourcePoints.empty())
		{
			return;
		}

		std::vector<pmp::Point2> leftBucket;
		std::vector<pmp::Point2> rightBucket;
		leftBucket.reserve(m_SourcePoints.size());
		rightBucket.reserve(m_SourcePoints.size());

		for (const auto& pt : m_SourcePoints)
		{
			if (IsPointLeftOfLine2D(pt, edge))
			{
				leftBucket.push_back(pt);
				continue;
			}

			rightBucket.push_back(pt);
		}

		leftBucket.shrink_to_fit();
		rightBucket.shrink_to_fit();

		if (!leftBucket.empty())
		{

		}

		if (!rightBucket.empty())
		{

		}
	}


} // namespace Geometry