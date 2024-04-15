#pragma once

namespace Utils
{
    /// ==================================================================
    /// \brief Interface for file mapping wrapper.
    /// \class IFileMappingWrapper
    /// ================================================================== 
    class IFileMappingWrapper
    {
    public:
        virtual ~IFileMappingWrapper() = default;

        virtual const void* GetFileMemory() const = 0;
        virtual size_t GetFileSize() const = 0;
    };

} // namespace Utils