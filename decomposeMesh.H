// const word decompTypeName
// (
//     meshDict.found("decompositionType")
//   ? meshDict.get<word>("decompositionType")
//   : "faceCentre"
// );

tetDecomposer::decompositionType decompType
(
//     tetDecomposer::decompositionTypeNames.found(decompTypeName)
//   ? tetDecomposer::decompositionTypeNames.get(decompTypeName)
//   : tetDecomposer::decompositionType::FACE_CENTRE_TRIS
    tetDecomposer::decompositionType::FACE_CENTRE_TRIS
);

Info<< nl << "Decomposing hexahedrons using the decomposition type of "
    << tetDecomposer::decompositionTypeNames.get(decompType)
    << endl;

{
    tetDecomposer tetMesh(mesh);

    polyTopoChange meshMod(mesh);
    bitSet decomposeCell(identity(mesh.nCells()));

    tetMesh.setRefinement(decompType, decomposeCell, meshMod);

    meshMod.changeMesh(mesh, false);

}

#include "printMeshSummary.H"