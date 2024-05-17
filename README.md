# blockPolyMesh

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11212395.svg)](https://doi.org/10.5281/zenodo.11212395)

## Overview

`blockPolyMesh` is a structured polyhedral mesh generator, which is based on `blockMesh`. The $n$ cells of the structured hexahedral mesh are decomposed into $24n$ tetrahedrons and then converted to polyhedra following the algorithm in `polyDualMesh`.

*Key features:*

* structured polyhedral mesh
* built using blocks
* supports cell size grading
* supports curved block edges

Same as the `blockMesh` mesher, `blockPolyMesh` is well suited to simple geometries that can be described by a few blocks, but challenging to apply to cases with a large number of blocks due to book-keeping requirements, *i.e.*, the need to manage point connectivity and ordering.

## Compatibility

The source code of `blockPolyMesh` is developed and maintained for OpenFOAM-v2212.


## Usage

`blockPolyMesh [OPTIONS] <featureAngle>`

Arguments:

* *featureAngle*

    in degrees [0-180]


Options:

* `-case` *dir*

    Case directory (instead of current directory)

* `-concaveAngle` *degrees*

    Specify concave angle [0..180] (default: 30 degrees)

* `-concaveMultiCells`

    Split cells on concave boundary edges into multiple cells

* `-dict` *file*

    Alternative blockMeshDict

* `-doNotPreserveFaceZones`

    Disable the default behaviour of preserving faceZones by having multiple faces in between cells

* `-region` *name*

    Specify alternative mesh region

* `-splitAllFaces`

    Have multiple faces in between cells

* `-write-vtk`

    Write topology as VTU file and exit

* `-doc`

    Display documentation in browser

* `-help`

    Display short help and exit

* `-help-compat`

    Display compatibility options and exit

* `-help-full`

    Display full help and exit


Advanced options:

* `-debug-switch` *name=val*

    Set named DebugSwitch (default value: 1). (Can be used multiple times)

* `-info-switch` *name=val*

    Set named InfoSwitch (default value: 1). (Can be used multiple times)

* `-lib` *name*

    Additional library or library list to load. (Can be used multiple times)

* `-no-libs`

    Disable use of the controlDict 'libs' entry

* `-doc-source`

    Display source code in browser

* `-help-man`

    Display full help (manpage format) and exit

* `-help-notes`

    Display help notes (description) and exit


## Tutorials

Three tutorial cases are available at the `tutorials` directory.

|Case|Mesh|
|:---:|:---:|
|`box`|<img src="./tutorials/box/mesh.png" alt="mesh" width="50%">|
|`cylinder`|<img src="./tutorials/cylinder/mesh.png" alt="mesh" width="50%">|
|`sphere`|<img src="./tutorials/sphere/mesh.png" alt="mesh" width="50%">|

## Contributor(s)

* Dezhi Dai, Argonne National Laboratory, daid@anl.gov (Developer)
