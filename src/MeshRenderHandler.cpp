#include "MeshRenderHandler.h"

#include "IncrementalMeshBuilder.h"

using namespace IMB;

MeshRenderHandler::MeshRenderHandler(const MeshRenderFunction& renderFunc)
{
    auto& meshBuilder = IncrementalMeshBuilder::GetInstance();
    meshBuilder.SetRenderCallback(renderFunc);
}
