#pragma once

#include "pmp/SurfaceMesh.h"

#include "utils/IFileMappingWrapper.h"
#include "geometry/GeometryConversionUtils.h"

#include "IncrementalProgressUtils.h"
#include "PointCloudMeshingStrategies.h"

namespace IMB
{
    /// \brief A functor for rendering a mesh.
    using MeshRenderFunction = std::function<void(const Geometry::BaseMeshGeometryData&)>;

    /// ==================================================================
    /// \brief The main singleton class for incremental mesh building.
    /// \class IncrementalMeshBuilder
    /// ==================================================================
    class IncrementalMeshBuilder
    {
    public:
        static IncrementalMeshBuilder& GetInstance() 
        {
            static IncrementalMeshBuilder instance;
            return instance;
        }

        IncrementalMeshBuilder(const IncrementalMeshBuilder&) = delete;
        IncrementalMeshBuilder& operator=(const IncrementalMeshBuilder&) = delete;

        /// ==================================================================
        /// \brief Initializes the mesh builder with the given parameters.
        /// \param fileName                the file name of the mesh to be loaded.
        /// \param completionFrequency     the frequency of mesh completion.
        /// \param reconstructType         the type of mesh reconstruction function.
        /// \param vertSelType             the type of vertex selection function.
        /// \param renderCallback          gets called when the updated mesh is ready.
        /// ==================================================================
        void Init(const std::string& fileName, const unsigned int& completionFrequency, 
            const ReconstructionFunctionType& reconstructType = ReconstructionFunctionType::BallPivoting, 
            const VertexSelectionType& vertSelType = VertexSelectionType::UniformRandom,
            const MeshRenderFunction& renderCallback = [](const Geometry::BaseMeshGeometryData&) {});

        /// ==================================================================
        /// \brief Samples vertices from the mesh using m_Dispatcher->m_VertexSamplingStrategy.
        /// \param[in] seed       the seed value for the vertex sampler.
        /// \param[in] nThreads   the preference for the amount of threads used. 0: means one thread will be used, >= nAvailableThreads: means nAvailableThreads will be used.
        /// ================================================================== 
        void DispatchAndSyncWorkers(const std::optional<unsigned int>& seed = std::nullopt, const unsigned int& nThreads = 0);

    private:
        /// \brief default constructor.
        IncrementalMeshBuilder() = default;

        /// \brief Reconstructs the mesh using m_MeshingStrategy.
        void UpdateMesh(const std::vector<pmp::Point>& newVertices);

        /// \brief Invoke destructor of all owned objects including the file mapping.
        void Terminate();

        //
        // ================================================================
        //

        std::atomic<bool> m_IsWorking{ false }; //>! a safety status flag to avoid re-dispatching threads before they're finished.

        bool m_IsInitialized{ false }; //>! initialization flag. Need to call Init if false.

        Geometry::BaseMeshGeometryData m_MeshData{}; //>! mesh data structure.
        std::mutex m_MeshDataMutex{};                //>! Mutex for protecting m_MeshData
        std::unique_ptr<PointCloudMeshingStrategy> m_MeshingStrategy{ nullptr }; //>! a strategy to convert point cloud to mesh.

        std::unique_ptr<Utils::IFileMappingWrapper> m_FileMapping{ nullptr }; //>! file mapping wrapper.
        std::unique_ptr<IncrementalMeshBuilderDispatcher> m_Dispatcher{nullptr}; //>! mesh builder dispatcher.
        MeshRenderFunction m_RenderCallback{}; //>! mesh render callback.
    };
} // namespace IMB