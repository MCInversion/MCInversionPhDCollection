#include "FileMappingWrapper.h"

#include <iostream>

using namespace Utils;

FileMappingWrapper::FileMappingWrapper(const std::string& filePath)
{
    OpenFile(filePath);
}

FileMappingWrapper::~FileMappingWrapper() 
{
    if (FileMemory) UnmapViewOfFile(FileMemory);
    if (FileMapping) CloseHandle(FileMapping);
    if (FileHandle != INVALID_HANDLE_VALUE) CloseHandle(FileHandle);
}

const void* FileMappingWrapper::GetFileMemory() const
{
    return FileMemory;
}

size_t FileMappingWrapper::GetFileSize() const
{
    return static_cast<size_t>(FileSize.QuadPart);
}

void FileMappingWrapper::OpenFile(const std::string& filePath)
{
    // Open the file with GENERIC_READ access and FILE_SHARE_READ mode.
    FileHandle = CreateFile(filePath.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (FileHandle == INVALID_HANDLE_VALUE)
    {
        std::cerr << "FileMappingWrapper::OpenFile: Failed to open file.";
        return;
    }

    // Retrieve the size of the file.
    if (!GetFileSizeEx(FileHandle, &FileSize)) 
    {
        CloseHandle(FileHandle);
        std::cerr << "FileMappingWrapper::OpenFile: Failed to get file size.";
        return;
    }

    // Create the file mapping object.
    FileMapping = CreateFileMapping(FileHandle, nullptr, PAGE_READONLY, FileSize.HighPart, FileSize.LowPart, nullptr);
    if (!FileMapping)
    {
        CloseHandle(FileHandle);
        std::cerr << "FileMappingWrapper::OpenFile: Failed to create file mapping.";
        return;
    }

    // Map a view of the file into the address space of the calling process.
    FileMemory = MapViewOfFile(FileMapping, FILE_MAP_READ, 0, 0, 0);
    if (!FileMemory) 
    {
        CloseHandle(FileMapping);
        CloseHandle(FileHandle);
        std::cerr << "FileMappingWrapper::OpenFile: Failed to map view of file.";
    }
}