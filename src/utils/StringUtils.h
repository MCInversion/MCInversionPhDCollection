#pragma once

#include <string>

namespace Utils
{
	/**
	 * \brief Extracts lowercase extension of a given file path string.
	 * \param pathStr    input string.
	 * \return lowercase extension of a given file path string.
	 */
	std::string ExtractLowercaseFileExtensionFromPath(const std::string& pathStr);
	
} // namespace Utils