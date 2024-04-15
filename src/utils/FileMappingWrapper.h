#pragma once

#include "IFileMappingWrapper.h"

#ifdef _WINDOWS
// Windows-specific headers
#include <windows.h>
#include <fcntl.h>
#else
// Unsupported platform
#error "Unsupported platform"
#endif

#include <string>

namespace Utils
{
    /// ==================================================================
    /// \brief Implementation of the file mapping wrapper.
    /// \class FileMappingWrapper
    /// ================================================================== 
    class FileMappingWrapper : public IFileMappingWrapper
    {
    public:
        FileMappingWrapper(const std::string& filePath);
        virtual ~FileMappingWrapper();

        const void* GetFileMemory() const override;
        size_t GetFileSize() const override;

    private:
        HANDLE FileHandle = INVALID_HANDLE_VALUE;
        HANDLE FileMapping = nullptr;
        LPVOID FileMemory = nullptr;
        LARGE_INTEGER FileSize;

        void OpenFile(const std::string& filePath);
    };

} // namespace Utils