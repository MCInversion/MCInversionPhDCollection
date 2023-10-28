#include "StringUtils.h"

#include <iostream>
#include <algorithm>

std::string Utils::ExtractLowercaseFileExtensionFromPath(const std::string& pathStr)
{
    const std::string::size_type dot(pathStr.rfind('.'));
    if (dot == std::string::npos)
    {
        std::cerr << "Utils::ExtractFileExtensionFromPath [ERROR]: dot not found!\n";
        return "";
    }
    std::string ext = pathStr.substr(dot + 1, pathStr.length() - dot - 1);
    std::ranges::transform(ext, std::back_inserter(ext), tolower);

	return pathStr.substr(pathStr.find_last_of("/\\") + 1);
}
