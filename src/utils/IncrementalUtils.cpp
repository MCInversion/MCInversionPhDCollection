
#include "IncrementalUtils.h"

#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>   // for sysconf on POSIX
#include <limits>     // for std::numeric_limits
#endif

bool HasEnoughMemoryForFisherYates(size_t totalCount)
{
    // Check for overflow when computing neededBytes = totalCount * sizeof(size_t)
    if (totalCount > (std::numeric_limits<std::size_t>::max() / sizeof(std::size_t))) {
        return false;
    }
    std::size_t neededBytes = totalCount * sizeof(std::size_t);

#if defined(_WIN32)
    // On Windows, use GlobalMemoryStatusEx to get available physical memory in bytes
    MEMORYSTATUSEX memInfo;
    memInfo.dwLength = sizeof(MEMORYSTATUSEX);
    if (!GlobalMemoryStatusEx(&memInfo)) 
    {
        // If the call fails, conservatively return false
        return false;
    }
    // ullAvailPhys is the amount of physical memory currently available, in bytes
    unsigned long long availableBytes = memInfo.ullAvailPhys;

    // Leave safety margin: require neededBytes <= 50% of available physical RAM
    unsigned long long threshold = availableBytes / 2ULL;
    return (static_cast<unsigned long long>(neededBytes) <= threshold);

#else
    // On POSIX, use sysconf to get number of available pages and page size
    long pages = sysconf(_SC_AVPHYS_PAGES);
    long pageSize = sysconf(_SC_PAGESIZE);
    if (pages < 0 || pageSize < 0) 
    {
        // sysconf failed; be conservative
        return false;
    }
    unsigned long long availableBytes = static_cast<unsigned long long>(pages) * static_cast<unsigned long long>(pageSize);

    // Leave safety margin: require neededBytes <= 50% of available RAM
    unsigned long long threshold = availableBytes / 2ULL;
    return (static_cast<unsigned long long>(neededBytes) <= threshold);
#endif
}
