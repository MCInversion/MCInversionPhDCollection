#include "StringUtils.h"

#include <iostream>
#include <algorithm>

std::string Utils::ExtractLowercaseFileExtensionFromPath(const std::string& pathStr)
{
    // Find the last dot after the last directory separator
    const size_t lastSlash = pathStr.find_last_of("/\\");
    const size_t dot = pathStr.rfind('.');

    if (dot == std::string::npos || (lastSlash != std::string::npos && dot < lastSlash)) 
    {
        std::cerr << "Utils::ExtractFileExtensionFromPath [ERROR]: dot not found or in the wrong place!\n";
        return "";
    }

    // Extract the extension
    std::string ext = pathStr.substr(dot + 1);
    std::string lowerExt;
    std::ranges::transform(ext, std::back_inserter(lowerExt),
        [](unsigned char c) -> char { return std::tolower(c); });

    return lowerExt;
}
