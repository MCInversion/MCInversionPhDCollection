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
        explicit FileMappingWrapper(const std::string& filePath);

        ~FileMappingWrapper() override;

        [[nodiscard]] char* GetFileMemory() const override;

        [[nodiscard]] size_t GetFileSize() const override;

    private:
        void OpenFile(const std::string& filePath);

        HANDLE m_FileHandle = INVALID_HANDLE_VALUE;
        HANDLE m_FileMapping = nullptr;
        LPVOID m_FileMemory = nullptr;
        LARGE_INTEGER m_FileSize;
    };

} // namespace Utils