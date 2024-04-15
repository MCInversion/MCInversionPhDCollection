#pragma once

#include "pmp/SurfaceMesh.h"

#include "utils/IFileMappingWrapper.h"

#include "MeshUpdateHandler.h"

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

    /// \brief A functor for mesh reconstruction.
    using ReconstructionFunction = std::function<void(std::vector<unsigned int>&, const std::vector<pmp::Point>&, const std::vector<pmp::Normal>&)>;

    /// \brief A functor for vertex selection from a file mapping range.
    using VertexSelectionFunction = std::function<void(const char* /* start */, const char* /* end */, std::vector<pmp::Point>& /* result */, const std::optional<unsigned int>& /* seed */)>;

    /// \brief A functor for vertex and normal selection from a file mapping range.
    //using VertexAndNormalSelectionFunction = std::function<void(const char*, const char*, std::vector<pmp::Point>&, std::vector<pmp::Normal>&)>;

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

        /// \brief triangulates the mesh using m_ReconstructPtCloud.
        void Triangulate();

        /// \brief Samples vertices from the mesh using m_SelectVertices.
        void SampleVertices();

        /// \brief Samples vertices and normals from the mesh using m_SelectVerticesAndNormals.
        //void SampleVerticesAndNormals();

        // TODO: New data structure for mesh data?
        /// \brief Refines the mesh using m_MeshData.
        //void RefineMesh();

        Geometry::BaseMeshGeometryData m_MeshData; //>! mesh data structure.

        ReconstructionFunction m_ReconstructPtCloud;                 //>! mesh reconstruction function.
        VertexSelectionFunction m_SelectVertices;                    //>! vertex selection function.
        //VertexAndNormalSelectionFunction m_SelectVerticesAndNormals; //>! vertex and normal selection function.

        std::unique_ptr<Utils::IFileMappingWrapper> m_FileMapping; //>! file mapping wrapper.

        MeshRenderFunction m_RenderCallback; //>! mesh render callback.
    };
} // namespace IMB