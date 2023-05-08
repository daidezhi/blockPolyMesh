// argList::addVerboseOption
// (
//     "Force verbose output. (Can be used multiple times)"
// );

argList::addOption("dict", "file", "Alternative blockMeshDict");

// argList::addOption
// (
//     "decompType",
//     "name",
//     "Specify tet decomposition type"
// );

argList::addBoolOption
(
    "write-vtk",
    "Write topology as VTU file and exit"
);

argList::addArgument
(
    "featureAngle",
    "in degrees [0-180]"
);

argList::addOption
(
    "concaveAngle",
    "degrees",
    "Specify concave angle [0..180] (default: 30 degrees)"
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