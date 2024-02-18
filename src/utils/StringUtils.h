#pragma once

#include <string>
#include <unordered_map>

namespace Utils
{
	/**
	 * \brief Extracts lowercase extension of a given file path string.
	 * \param pathStr    input string.
	 * \return lowercase extension of a given file path string.
	 */
	std::string ExtractLowercaseFileExtensionFromPath(const std::string& pathStr);

	/// \brief prints the face intersection ids multimap.
	void PrintFaceIntersectionsMultimap(const std::unordered_multimap<unsigned int, unsigned int>& faceIntersections);
	
} // namespace Utils