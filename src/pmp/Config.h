#pragma once

#include <cstdint>

// --------------------------------------------------------------
/// Since the original macros like PMP_SCALAR_TYPE_64 in Types.h didn't work
/// because it collided with MatVec.h types, we need to extract them to a 
/// separate header. I doubt the original PMP implementation ever used these
/// definitions.
// --------------------------------------------------------------

// Scalar type configuration
#define PMP_SCALAR_TYPE_64 false
// Index type configuration
#define PMP_INDEX_TYPE_64 false

//! The pmp-library namespace
namespace pmp {

#if PMP_SCALAR_TYPE_64
	using Scalar = double;
#else
	using Scalar = float;
#endif

	// define index type to be used
#if PMP_INDEX_TYPE_64
	using IndexType = std::uint_least64_t;
#define PMP_MAX_INDEX UINT_LEAST64_MAX
#else
	using IndexType = std::uint_least32_t;
#define PMP_MAX_INDEX UINT_LEAST32_MAX
#endif

} // namespace pmp