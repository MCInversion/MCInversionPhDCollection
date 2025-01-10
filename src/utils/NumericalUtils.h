#pragma once

#include <cmath>
#include <ranges>
#include <stdexcept>

namespace Utils
{
	/// \brief A utility for extracting the 32-bit remainder of a 34-bit value
	[[nodiscard]] float Get32BitRemainder(const double& value);

    /// \brief A utility for generating an exponential range of values
    template <typename T>
    [[nodiscard]] auto GetExponentialValueRange(const T& start, const T& end, const T& factor)
	{
        if (factor <= 0 || start <= 0 || end <= 0 || (factor == 1 && start < end)) 
        {
            throw std::invalid_argument("Invalid range parameters");
        }

        return std::ranges::views::iota(0)
            | std::ranges::views::transform([=](int i) { return start * std::pow(factor, i); })
            | std::ranges::views::take_while([=](T value) { return value > end; });
    }

    /// \brief A utility for generating a linear range of values
    template <typename T>
    [[nodiscard]] auto GetLinearValueRange(const T& start, const T& end, const T& step) {
        if (step == 0 || (start < end && step < 0) || (start > end && step > 0)) {
            throw std::invalid_argument("Invalid range parameters");
        }

        return std::ranges::views::iota(0)
            | std::ranges::views::transform([=](int i) { return start + i * step; })
            | std::ranges::views::take_while([=](T value) {
            return (step > 0) ? value < end : value > end;
                });
    }


} // namespace Utils