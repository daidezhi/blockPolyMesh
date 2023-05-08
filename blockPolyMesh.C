/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
    Copyright (C) 2023      Dezhi Dai
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
    blockPolyMesh

Description
    A structured polyhedral mesh generator.

    Uses the block mesh description found in
      - \c system/blockMeshDict
      - \c system/\<region\>/blockMeshDict
      - \c constant/polyMesh/blockMeshDict
      - \c constant/\<region\>/polyMesh/blockMeshDict

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
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "combineFaces.H"
#include "removePoints.H"
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

#include "functions.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Structured polyhedral mesh generator, "
        "which is based on `blockMesh`. The $n$ cells of the structured "
        "hexahedral mesh are decomposed into $24n$ tetrahedrons and then "
        "converted to polyhedra following the algorithm in `polyDualMesh`.\n"
    );

    argList::noParallel();
    argList::noFunctionObjects();

    #include "addOptions.H"

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
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //