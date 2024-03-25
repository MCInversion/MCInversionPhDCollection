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
	[[nodiscard]] std::string ExtractLowercaseFileExtensionFromPath(const std::string& pathStr);

	/// \brief prints the face intersection ids multimap.
	void PrintFaceIntersectionsMultimap(const std::unordered_multimap<unsigned int, unsigned int>& faceIntersections);

	/// \brief outputs a string from a given index.
	inline [[nodiscard]] std::string FormatIndexSimple(const unsigned int& i)
	{
		return std::to_string(i);
	}

	/// \brief outputs a string from a given index, such that its length is always 4, and the missing digits are filled with zeros.
	[[nodiscard]] std::string FormatIndex4DigitFill(const unsigned int& i);

} // namespace Utils