#include "FileMappingWrapper.h"

#include <iostream>

namespace Utils
{
	FileMappingWrapper::FileMappingWrapper(const std::string& filePath)
	{
	    OpenFile(filePath);
	}

	FileMappingWrapper::~FileMappingWrapper()
	{
	    if (m_FileMemory) UnmapViewOfFile(m_FileMemory);
	    if (m_FileMapping) CloseHandle(m_FileMapping);
	    if (m_FileHandle != INVALID_HANDLE_VALUE) CloseHandle(m_FileHandle);
	}

	char* FileMappingWrapper::GetFileMemory() const
	{
	    return static_cast<char*>(m_FileMemory);
	}

	size_t FileMappingWrapper::GetFileSize() const
	{
	    return static_cast<size_t>(m_FileSize.QuadPart);
	}

	size_t FileMappingWrapper::GetLineCount() const
	{
		char* fileStart = GetFileMemory();
		char* fileEnd = fileStart + GetFileSize();
		return std::count(fileStart, fileEnd, '\n');
	}

	void FileMappingWrapper::OpenFile(const std::string& filePath)
	{
	    // Open the file with GENERIC_READ access and FILE_SHARE_READ mode.
		m_FileHandle = CreateFile(filePath.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
	    if (m_FileHandle == INVALID_HANDLE_VALUE)
	    {
	        std::cerr << "FileMappingWrapper::OpenFile: Failed to open file.";
	        return;
	    }

	    // Retrieve the size of the file.
	    if (!GetFileSizeEx(m_FileHandle, &m_FileSize))
	    {
	        CloseHandle(m_FileHandle);
	        std::cerr << "FileMappingWrapper::OpenFile: Failed to get file size.";
	        return;
	    }

	    // Create the file mapping object.
		m_FileMapping = CreateFileMapping(m_FileHandle, nullptr, PAGE_READONLY, m_FileSize.HighPart, m_FileSize.LowPart, nullptr);
	    if (!m_FileMapping)
	    {
	        CloseHandle(m_FileHandle);
	        std::cerr << "FileMappingWrapper::OpenFile: Failed to create file mapping.";
	        return;
	    }

	    // Map a view of the file into the address space of the calling process.
		m_FileMemory = MapViewOfFile(m_FileMapping, FILE_MAP_READ, 0, 0, 0);
	    if (!m_FileMemory)
	    {
	        CloseHandle(m_FileMapping);
	        CloseHandle(m_FileHandle);
	        std::cerr << "FileMappingWrapper::OpenFile: Failed to map view of file.";
	    }
	}
	
} // namespace utils
