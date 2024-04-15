#include "MeshUpdateHandler.h"

#include "IncrementalMeshBuilder.h"

using namespace IMB;

MeshUpdateHandler::MeshUpdateHandler(const MeshRenderFunction& renderFunc)
{
    auto& meshBuilder = IncrementalMeshBuilder::GetInstance();
    meshBuilder.SetRenderCallback(renderFunc);
}
