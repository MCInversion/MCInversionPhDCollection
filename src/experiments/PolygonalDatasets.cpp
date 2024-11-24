#include "PolygonalDatasets.h"

#include <vector>
#include <sstream>
#include <cctype>
#include <iterator>
#include <regex>

namespace
{
    [[nodiscard]] std::string ExtractDAttribute(const std::string& svgPath) 
    {
        std::regex dAttrRegex(R"d(d="([^"]*)")d"); // Match d="...".
        std::smatch match;
        if (std::regex_search(svgPath, match, dAttrRegex))
        {
            return match[1].str();
        }
        throw std::runtime_error("Failed to extract 'd' attribute from SVG path.");
    }

    [[nodiscard]] // Tokenizer: Split a string by whitespace or commas
    std::vector<std::string> Tokenize(const std::string& str) {
        std::vector<std::string> tokens;
        std::string current;

        for (char ch : str) {
            if (std::isspace(ch) || ch == ',') {
                if (!current.empty()) {
                    tokens.push_back(current);
                    current.clear();
                }
            }
            else {
                current += ch;
            }
        }

        if (!current.empty()) {
            tokens.push_back(current);
        }

        return tokens;
    }

    [[nodiscard]] std::vector<pmp::Point2> ParseDAttribute(const std::string& dData)
    {
        std::vector<pmp::Point2> points;
        float currentX = 0.0f, currentY = 0.0f;
        float startX = 0.0f, startY = 0.0f;

        auto tokens = Tokenize(dData);
        char currentCommand = 0;

        size_t i = 0;
        while (i < tokens.size()) {
            const std::string& token = tokens[i];

            if (std::isalpha(token[0])) {
                // If the token is a command (e.g., 'M', 'L')
                currentCommand = token[0];
                ++i;
            }

            switch (currentCommand) {
            case 'M': { // Move to (absolute)
                float x = std::stof(tokens[i++]);
                float y = std::stof(tokens[i++]);
                currentX = x;
                currentY = y;
                startX = currentX;
                startY = currentY;
                points.emplace_back(currentX, currentY);
                currentCommand = 'L'; // Treat subsequent coordinates as line-to
                break;
            }
            case 'm': { // Move to (relative)
                float dx = std::stof(tokens[i++]);
                float dy = std::stof(tokens[i++]);
                currentX += dx;
                currentY += dy;
                startX = currentX;
                startY = currentY;
                points.emplace_back(currentX, currentY);
                currentCommand = 'l'; // Treat subsequent coordinates as relative line-to
                break;
            }
            case 'L': { // Line to (absolute)
                while (i < tokens.size() && !std::isalpha(tokens[i][0])) {
                    float x = std::stof(tokens[i++]);
                    float y = std::stof(tokens[i++]);
                    currentX = x;
                    currentY = y;
                    points.emplace_back(currentX, currentY);
                }
                break;
            }
            case 'l': { // Line to (relative)
                while (i < tokens.size() && !std::isalpha(tokens[i][0])) {
                    float dx = std::stof(tokens[i++]);
                    float dy = std::stof(tokens[i++]);
                    currentX += dx;
                    currentY += dy;
                    points.emplace_back(currentX, currentY);
                }
                break;
            }
            case 'V': { // Vertical line to (absolute)
                while (i < tokens.size() && !std::isalpha(tokens[i][0])) {
                    float y = std::stof(tokens[i++]);
                    currentY = y;
                    points.emplace_back(currentX, currentY);
                }
                break;
            }
            case 'v': { // Vertical line to (relative)
                while (i < tokens.size() && !std::isalpha(tokens[i][0])) {
                    float dy = std::stof(tokens[i++]);
                    currentY += dy;
                    points.emplace_back(currentX, currentY);
                }
                break;
            }
            case 'H': { // Horizontal line to (absolute)
                while (i < tokens.size() && !std::isalpha(tokens[i][0])) {
                    float x = std::stof(tokens[i++]);
                    currentX = x;
                    points.emplace_back(currentX, currentY);
                }
                break;
            }
            case 'h': { // Horizontal line to (relative)
                while (i < tokens.size() && !std::isalpha(tokens[i][0])) {
                    float dx = std::stof(tokens[i++]);
                    currentX += dx;
                    points.emplace_back(currentX, currentY);
                }
                break;
            }
            case 'Z': // Close path (absolute)
            case 'z': // Close path (relative)
                points.emplace_back(startX, startY);
                ++i;
                break;
            default:
                throw std::runtime_error("Unsupported SVG command: " + std::string(1, currentCommand));
            }
        }

        return points;
    }
} // anonymous namespace

std::vector<pmp::Point2> ParsePolygonalSVGPath(const std::string& svgPath)
{
    // Extract the d attribute
    std::string dData = ExtractDAttribute(svgPath);

    // Parse the d attribute data into points
    // TODO: fix potential duplicate points
    return ParseDAttribute(dData);
}
