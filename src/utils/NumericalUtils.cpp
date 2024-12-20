#include "NumericalUtils.h"

#include <cmath>

namespace Utils
{
	float Get32BitRemainder(const double& value)
	{
        // Define the divisor, which is 2^32
        constexpr auto divisor = static_cast<double>(1ULL << 32);

        // Calculate the remainder using std::fmod
        double remainder = std::fmod(value, divisor);

        // If the remainder is negative, adjust it to be positive
        if (remainder < 0)
        {
            remainder += divisor;
        }

        // Cast the result to a float and return
        return static_cast<float>(remainder);
	}

} // namespace Utils