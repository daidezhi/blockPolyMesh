/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockMesh

Group
    grpMeshGenerationUtilities

Description
    A multi-block mesh generator.

    Uses the block mesh description found in
      - \c system/blockMeshDict
      - \c system/\<region\>/blockMeshDict
      - \c constant/polyMesh/blockMeshDict
      - \c constant/\<region\>/polyMesh/blockMeshDict

Usage
    \b blockMesh [OPTION]

    Options:
      - \par -write-obj
        Write topology as a set of edges in OBJ format and exit.

      - \par -write-vtk
        Write topology as VTK file (xml, ascii) and exit.

      - \par -merge-points
        Merge points instead of default topological merge

      - \par -region \<name\>
        Specify alternative mesh region.

      - \par -dict \<filename\>
        Alternative dictionary for the block mesh description.

      - \par -sets
        Write cellZones as cellSets too (for processing purposes)

      - \par -no-clean
        Do not remove polyMesh/ directory or files

      - \par -time
        Write resulting mesh to a time directory (instead of constant)

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"
#include "unitConversion.H"
#include "blockMesh.H"
#include "foamVtkInternalMeshWriter.H"
#include "foamVtkSurfaceWriter.H"
#include "attachPolyTopoChanger.H"
#include "polyTopoChange.H"
#include "tetDecomposer.H"
#include "cyclicPolyPatch.H"
#include "cellSet.H"

#include "meshDualiser.H"
#include "meshTools.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "wordPair.H"
#include "slidingInterface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Naive feature detection. All boundary edges with angle > featureAngle become
// feature edges. All points on feature edges become feature points. All
// boundary faces become feature faces.
void simpleMarkFeatures
(
    const polyMesh& mesh,
    const bitSet& isBoundaryEdge,
    const scalar featureAngle,
    const bool concaveMultiCells,
    const bool doNotPreserveFaceZones,

    labelList& featureFaces,
    labelList& featureEdges,
    labelList& singleCellFeaturePoints,
    labelList& multiCellFeaturePoints
)
{
    scalar minCos = Foam::cos(degToRad(featureAngle));

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Working sets
    labelHashSet featureEdgeSet;
    labelHashSet singleCellFeaturePointSet;
    labelHashSet multiCellFeaturePointSet;


    // 1. Mark all edges between patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        const labelList& meshEdges = pp.meshEdges();

        // All patch corner edges. These need to be feature points & edges!
        for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
        {
            label meshEdgeI = meshEdges[edgeI];
            featureEdgeSet.insert(meshEdgeI);
            singleCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][0]);
            singleCellFeaturePointSet.insert(mesh.edges()[meshEdgeI][1]);
        }
    }



    // 2. Mark all geometric feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Make distinction between convex features where the boundary point becomes
    // a single cell and concave features where the boundary point becomes
    // multiple 'half' cells.

    // Addressing for all outside faces
    primitivePatch allBoundary
    (
        SubList<face>
        (
            mesh.faces(),
            mesh.nBoundaryFaces(),
            mesh.nInternalFaces()
        ),
        mesh.points()
    );

    // Check for non-manifold points (surface pinched at point)
    allBoundary.checkPointManifold(false, &singleCellFeaturePointSet);

    // Check for non-manifold edges (surface pinched at edge)
    const labelListList& edgeFaces = allBoundary.edgeFaces();
    const labelList& meshPoints = allBoundary.meshPoints();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            const edge& e = allBoundary.edges()[edgeI];

            //Info<< "Detected non-manifold boundary edge:" << edgeI
            //    << " coords:"
            //    << allBoundary.points()[meshPoints[e[0]]]
            //    << allBoundary.points()[meshPoints[e[1]]] << endl;

            singleCellFeaturePointSet.insert(meshPoints[e[0]]);
            singleCellFeaturePointSet.insert(meshPoints[e[1]]);
        }
    }

    // Check for features.
    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() == 2)
        {
            label f0 = eFaces[0];
            label f1 = eFaces[1];

            // check angle
            const vector& n0 = allBoundary.faceNormals()[f0];
            const vector& n1 = allBoundary.faceNormals()[f1];

            if ((n0 & n1) < minCos)
            {
                const edge& e = allBoundary.edges()[edgeI];
                label v0 = meshPoints[e[0]];
                label v1 = meshPoints[e[1]];

                label meshEdgeI = meshTools::findEdge(mesh, v0, v1);
                featureEdgeSet.insert(meshEdgeI);

                // Check if convex or concave by looking at angle
                // between face centres and normal
                vector c1c0
                (
                    allBoundary[f1].centre(allBoundary.points())
                  - allBoundary[f0].centre(allBoundary.points())
                );

                if (concaveMultiCells && (c1c0 & n0) > SMALL)
                {
                    // Found concave edge. Make into multiCell features
                    Info<< "Detected concave feature edge:" << edgeI
                        << " cos:" << (c1c0 & n0)
                        << " coords:"
                        << allBoundary.points()[v0]
                        << allBoundary.points()[v1]
                        << endl;

                    singleCellFeaturePointSet.erase(v0);
                    multiCellFeaturePointSet.insert(v0);
                    singleCellFeaturePointSet.erase(v1);
                    multiCellFeaturePointSet.insert(v1);
                }
                else
                {
                    // Convex. singleCell feature.
                    if (!multiCellFeaturePointSet.found(v0))
                    {
                        singleCellFeaturePointSet.insert(v0);
                    }
                    if (!multiCellFeaturePointSet.found(v1))
                    {
                        singleCellFeaturePointSet.insert(v1);
                    }
                }
            }
        }
    }


    // 3. Mark all feature faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Face centres that need inclusion in the dual mesh
    labelHashSet featureFaceSet(mesh.nBoundaryFaces());
    // A. boundary faces.
    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        featureFaceSet.insert(facei);
    }

    // B. face zones.
    const faceZoneMesh& faceZones = mesh.faceZones();

    if (doNotPreserveFaceZones)
    {
        if (faceZones.size() > 0)
        {
            WarningInFunction
                << "Detected " << faceZones.size()
                << " faceZones. These will not be preserved."
                << endl;
        }
    }
    else
    {
        if (faceZones.size() > 0)
        {
            Info<< "Detected " << faceZones.size()
                << " faceZones. Preserving these by marking their"
                << " points, edges and faces as features." << endl;
        }

        forAll(faceZones, zoneI)
        {
            const faceZone& fz = faceZones[zoneI];

            Info<< "Inserting all faces in faceZone " << fz.name()
                << " as features." << endl;

            forAll(fz, i)
            {
                label facei = fz[i];
                const face& f = mesh.faces()[facei];
                const labelList& fEdges = mesh.faceEdges()[facei];

                featureFaceSet.insert(facei);
                forAll(f, fp)
                {
                    // Mark point as multi cell point (since both sides of
                    // face should have different cells)
                    singleCellFeaturePointSet.erase(f[fp]);
                    multiCellFeaturePointSet.insert(f[fp]);

                    // Make sure there are points on the edges.
                    featureEdgeSet.insert(fEdges[fp]);
                }
            }
        }
    }

    // Transfer to arguments
    featureFaces = featureFaceSet.toc();
    featureEdges = featureEdgeSet.toc();
    singleCellFeaturePoints = singleCellFeaturePointSet.toc();
    multiCellFeaturePoints = multiCellFeaturePointSet.toc();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Polyhedral mesh generator based on a block mesh.\n"
    );

    argList::noParallel();
    argList::noFunctionObjects();

    argList::addVerboseOption
    (
        "Force verbose output. (Can be used multiple times)"
    );

    argList::addOption("dict", "file", "Alternative blockMeshDict");

    // argList::addOption
    // (
    //     "decompType",
    //     "name",
    //     "Specify tet decomposition type"
    // );

    argList::addArgument
    (
        "featureAngle",
        "in degrees [0-180]"
    );

    argList::addBoolOption
    (
        "splitAllFaces",
        "Have multiple faces in between cells"
    );

    argList::addBoolOption
    (
        "concaveMultiCells",
        "Split cells on concave boundary edges into multiple cells"
    );

    argList::addBoolOption
    (
        "doNotPreserveFaceZones",
        "Disable the default behaviour of preserving faceZones by having"
        " multiple faces in between cells"
    );

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const word regionName
    (
        args.found("region")
      ? args.get<word>("region")
      : polyMesh::defaultRegion
    );
    Info << "Generating mesh for region: " << regionName << endl;

    word meshInstance(runTime.constant());

    #include "findBlockPolyMeshDict.H"

    #include "createBlockMesh.H"
    #include "decomposeMesh.H"
    #include "convertDualMesh.H"

    Info << nl << "Writing mesh to constant" << endl;

    mesh.setInstance(meshInstance);
    mesh.removeFiles();
    mesh.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
