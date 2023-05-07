Info<< nl << "Converting the dual of the tet mesh" << endl;

//const word oldInstance(mesh.pointsInstance());

// Mark boundary edges and points.
// (Note: in 1.4.2 we can use the built-in mesh point ordering
//  facility instead)
bitSet isBoundaryEdge(mesh.nEdges());
for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
{
    const labelList& fEdges = mesh.faceEdges()[facei];

    forAll(fEdges, i)
    {
        isBoundaryEdge.set(fEdges[i]);
    }
}

const scalar featureAngle(args.get<scalar>(1));
const scalar minCos(Foam::cos(degToRad(featureAngle)));

Info<< "Feature:" << featureAngle << endl
    << "minCos :" << minCos
    << endl;

const bool splitAllFaces(args.found("splitAllFaces"));
if (splitAllFaces)
{
    Info<< "Splitting all internal faces to create multiple faces"
        << " between two cells." << nl
        << endl;
}

const bool doNotPreserveFaceZones(args.found("doNotPreserveFaceZones"));
const bool concaveMultiCells(args.found("concaveMultiCells"));
if (concaveMultiCells)
{
    Info<< "Generating multiple cells for points on concave feature edges."
        << endl;
}

// Face(centre)s that need inclusion in the dual mesh
labelList featureFaces;
// Edge(centre)s  ,,
labelList featureEdges;
// Points (that become a single cell) that need inclusion in the dual mesh
labelList singleCellFeaturePoints;
// Points (that become a multiple cells)
labelList multiCellFeaturePoints;

// Sample implementation of feature detection
simpleMarkFeatures
(
    mesh,
    isBoundaryEdge,
    featureAngle,
    concaveMultiCells,
    doNotPreserveFaceZones,

    featureFaces,
    featureEdges,
    singleCellFeaturePoints,
    multiCellFeaturePoints
);

// If we want to split all polyMesh faces into one dualface per cell
// we are passing through we also need a point
// at the polyMesh facecentre and edgemid of the faces we want to
// split.
if (splitAllFaces)
{
    featureEdges = identity(mesh.nEdges());
    featureFaces = identity(mesh.nFaces());
}

// Topo change container
polyTopoChange dualMeshMod(mesh.boundaryMesh().size());

// Mesh dualiser engine
meshDualiser dualMaker(mesh);

// Insert all commands into polyTopoChange to create dual of mesh. This does
// all the hard work.
dualMaker.setRefinement
(
    splitAllFaces,
    featureFaces,
    featureEdges,
    singleCellFeaturePoints,
    multiCellFeaturePoints,
    dualMeshMod
);

// Create mesh, return map from old to new mesh
//autoPtr<mapPolyMesh> map(dualMeshMod.changeMesh(mesh, false));
dualMeshMod.changeMesh(mesh, false);

// Update fields
//mesh.updateMesh(map());

#include "printMeshSummary.H"