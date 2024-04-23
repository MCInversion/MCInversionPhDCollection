#pragma once

// ======================================================================
// WARNING: Remove if you want IncrementalMeshBuilder to perform.
// ======================================================================
#define DEBUG_PRINT true

#if DEBUG_PRINT
#include "TimingUtils.h"
#define DBG_OUT std::cout << "[" << GetCurrentTimestamp() << "] "
#endif

#include <cstdint>
#include <sstream>
#include <iomanip>

inline [[nodiscard]] std::string FormatAddresses(const char* start, const char* end)
{
    std::ostringstream ss;
    ss << "0x" << std::hex << std::setw(16) << std::setfill('0') << reinterpret_cast<uintptr_t>(start)
        << ", 0x" << std::hex << std::setw(16) << std::setfill('0') << reinterpret_cast<uintptr_t>(end);
    return ss.str();
}

constexpr unsigned int DEFAULT_MAX_VERTEX_CAPACITY = 1'000'000;