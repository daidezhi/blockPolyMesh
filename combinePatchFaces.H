{
    // Sin of angle between two consecutive edges on a face.
    // If sin(angle) larger than this the face will be considered concave.
    const scalar concaveAngle = args.getOrDefault<scalar>("concaveAngle", 30);

    const scalar concaveSin = Foam::sin(degToRad(concaveAngle));

    Info<< "Merging all faces of a cell" << nl
        << "    - which are on the same patch" << nl
        << "    - which make an angle < " << featureAngle << " degrees"
        << nl
        << "      (cos:" << minCos << ')' << nl
        << "    - even when resulting face becomes concave by more than "
        << concaveAngle << " degrees" << nl
        << "      (sin:" << concaveSin << ')' << nl
        << endl;

    // Merge faces on same patch
    label nChanged = mergePatchFaces
    (
        minCos,
        concaveSin,
        runTime,
        mesh
    );

    nChanged += mergeEdges(minCos, mesh);
}