#pragma once

// ======================================================================
// WARNING: Remove if you want IncrementalMeshBuilder to perform.
// ======================================================================
#define DEBUG_PRINT true

#if DEBUG_PRINT
#include "TimingUtils.h"

#include <iomanip>
#include <sstream>
#include <cstdint>
#endif

#define DBG_OUT std::cout << "[" << GetCurrentTimestamp() << "] "

#if DEBUG_PRINT
inline [[nodiscard]] std::string FormatAddresses(const char* start, const char* end)
{
    std::ostringstream ss;
    ss << "0x" << std::hex << std::setw(16) << std::setfill('0') << reinterpret_cast<uintptr_t>(start)
        << ", 0x" << std::hex << std::setw(16) << std::setfill('0') << reinterpret_cast<uintptr_t>(end);
    return ss.str();
}
#endif
