#pragma once

#include "pmp/SurfaceMesh.h"

#include "utils/IFileMappingWrapper.h"

#include "MeshUpdateHandler.h"
#include "IncrementalProgressUtils.h"

namespace IMB
{
	/// \brief enumerator for mesh reconstruction function type.
	enum class [[nodiscard]] ReconstructionFunctionType
	{
		BallPivoting             = 0, //>! reconstructs a mesh using the ball-pivoting algorithm.
		Poisson                  = 1, //>! reconstructs a mesh using the Poisson surface reconstruction algorithm (requires normals).
		MarchingCubes            = 2, //>! reconstructs a mesh using the marching cubes algorithm.
		LagrangianShrinkWrapping = 3, //>! reconstructs a mesh using the Lagrangian shrink-wrapping algorithm.
	};

	/// \brief enumerator for mesh simplification function type.
	enum class [[nodiscard]] VertexSelectionType
	{
		Sequential              = 0, //>! selects vertices sequentially.
		UniformRandom           = 1, //>! selects vertices uniformly at random.
		NormalRandom            = 2, //>! selects vertices with a normal distribution.
		SoftMaxFeatureDetecting = 3, //>! selects vertices using a softmax function with feature detection.
	};

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

        void Init(const std::string& fileName, const unsigned int& completionFrequency, 
            const ReconstructionFunctionType& reconstructType = ReconstructionFunctionType::BallPivoting, 
            const VertexSelectionType& vertSelType = VertexSelectionType::UniformRandom);

        //void UpdateMesh();

        void SetRenderCallback(const MeshRenderFunction& renderCallback) 
        {
            m_RenderCallback = renderCallback; 
        }

    private:
        /// \brief default constructor.
        IncrementalMeshBuilder() = default;

        /// \brief Samples vertices from the mesh using m_SelectVertices.
        void ProcessVertices(const std::optional<unsigned int>& seed = std::nullopt, const unsigned int& nThreads = 0);

        void UpdateMesh(const std::vector<pmp::Point>& vertices, const std::vector<std::vector<unsigned int>>& polyIndices);

        //
        // ================================================================
        //

        Geometry::BaseMeshGeometryData m_MeshData; //>! mesh data structure.
        std::mutex m_MeshDataMutex;                 // Mutex for protecting m_MeshData

        std::unique_ptr<Utils::IFileMappingWrapper> m_FileMapping; //>! file mapping wrapper.
        std::unique_ptr<IncrementalMeshBuilderDispatcher> m_Dispatcher{nullptr}; //>! mesh builder dispatcher.
        MeshRenderFunction m_RenderCallback; //>! mesh render callback.
    };
} // namespace IMB