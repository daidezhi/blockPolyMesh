Info<< nl << "Creating a block mesh";

const blockMesh::mergeStrategy strategy
(
    args.found("merge-points")
  ? blockMesh::MERGE_POINTS
  : blockMesh::DEFAULT_MERGE
);

blockMesh blocks(meshDict, regionName, strategy, args.verbose());

if (!blocks.valid())
{
    // Could/should be Fatal?

    WarningIn(args.executable())
        << "Did not generate any blocks. Stopping." << nl << endl;

    return 1;
}

blocks.verbose(false);

autoPtr<polyMesh> meshPtr
(
    blocks.mesh(IOobject(regionName, meshInstance, runTime))
);

polyMesh& mesh(*meshPtr);

// Merge patch pairs (dictionary entry "mergePatchPairs")
#include "mergePatchPairs.H"

// Handle cyclic patches
#include "handleCyclicPatches.H"

// Set the precision of the points data to 10
IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

// mesh.removeFiles();
// if (!mesh.write())
// {
//     FatalErrorInFunction
//         << "Failed writing polyMesh."
//         << exit(FatalError);
// }

#include "printMeshSummary.H"