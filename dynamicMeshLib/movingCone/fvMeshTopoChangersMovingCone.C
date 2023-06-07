/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "fvMeshTopoChangersMovingCone.H"
#include "Time.H"
#include "polyTopoChangeMap.H"
#include "layerAdditionRemoval.H"
#include "meshTools.H"
#include "OFstream.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(movingCone, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, movingCone, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::fvMeshTopoChangers::movingCone::vertexMarkup
(
    const pointField& p,
    const scalar curLeft,
    const scalar curRight
) const
{
    Info<< "Updating vertex markup.  curLeft: "
        << curLeft << " curRight: " << curRight << endl;

    tmp<scalarField> tvertexMarkup(new scalarField(p.size(), 0));
    scalarField& vertexMarkup = tvertexMarkup.ref();
    
    word mc("movingCells");
    cellZoneID movingCellsID(mc, mesh().cellZones());
    
    if (movingCellsID.active())
    {
        // Get cell-point addressing
        const labelListList& cp = mesh().cellPoints();

        // Get labels of all moving cells
        const labelList& movingCells = mesh().cellZones()[movingCellsID.index()];

        forAll (movingCells, cellI)
        {
            const labelList& curCp = cp[movingCells[cellI]];

            forAll (curCp, pointI)
            {
                vertexMarkup[curCp[pointI]] = 1;
            }
        }
    }
    
    word rf("rightExtrusionFaces");
    faceZoneID frontFacesID(rf, mesh().faceZones());

    if (frontFacesID.active())
    {
        const faceZone& frontFaces = mesh().faceZones()[frontFacesID.index()];

        const labelList& mp = frontFaces().meshPoints();

        forAll (mp, mpI)
        {
            vertexMarkup[mp[mpI]] = 1;
        }
    }
    
    word lf("leftExtrusionFaces");
    faceZoneID backFacesID(lf, mesh().faceZones());

    if (backFacesID.active())
    {
        const faceZone& backFaces = mesh().faceZones()[backFacesID.index()];

        const labelList& mp = backFaces().meshPoints();

        forAll (mp, mpI)
        {
            vertexMarkup[mp[mpI]] = 1;
        }
    }

    return tvertexMarkup;
}


void Foam::fvMeshTopoChangers::movingCone::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        mesh().pointZones().size()
     || mesh().faceZones().size()
     || mesh().cellZones().size()
     || topoChanger_.size()
    )
    {
        InfoInFunction
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    Info<< "Time = " << mesh().time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    const vectorField& fc = mesh().faceCentres();
    const vectorField& fa = mesh().faceAreas();
    const vectorField& cc = mesh().cellCentres();

    labelList zone1(fc.size());
    boolList flipZone1(fc.size(), false);
    label nZoneFaces1 = 0;

    labelList zone2(fc.size());
    boolList flipZone2(fc.size(), false);
    label nZoneFaces2 = 0;
    
    labelList zone3(cc.size());
    label nZoneCells = 0;    

    forAll(fc, facei)
    {
        if
        (
            fc[facei].x() > curRight_ - delta_/2 + small
         && fc[facei].x() < curRight_ + delta_/2 - small
         && fc[facei].y() < curTop_ - small
        )
        {
            if ((fa[facei] & vector(1, 0, 0)) < 0)
            {
                flipZone1[nZoneFaces1] = true;
            }

            zone1[nZoneFaces1] = facei;
            Info<< "face " << facei << " for zone 1.  Flip: "
                << flipZone1[nZoneFaces1] << endl;
            nZoneFaces1++;
        }
        else if
        (
            fc[facei].x() > curLeft_ - delta_/2 + small
         && fc[facei].x() < curLeft_ + delta_/2 - small
         && fc[facei].y() < curTop_ - small
        )
        {
            zone2[nZoneFaces2] = facei;

            if ((fa[facei] & vector(1, 0, 0)) > 0)
            {
                flipZone2[nZoneFaces2] = true;
            }

            Info<< "face " << facei << " for zone 2.  Flip: "
                << flipZone2[nZoneFaces2] << endl;
            nZoneFaces2++;
        }
    }
    
    forAll(cc, celli)
    {
        if 
        (
            cc[celli].x() > curLeft_  + small
         && cc[celli].x() < curRight_ - small
         && cc[celli].y() < curTop_   - small
        )
        {
            zone3[nZoneCells] = celli;
            nZoneCells++;
        }
    }

    zone1.setSize(nZoneFaces1);
    flipZone1.setSize(nZoneFaces1);

    zone2.setSize(nZoneFaces2);
    flipZone2.setSize(nZoneFaces2);
    
    zone3.setSize(nZoneCells);

    List<pointZone*> pz(0);
    List<faceZone*> fz(2);
    List<cellZone*> cz(1);

    label nFz = 0;
    label nCz = 0;

    fz[nFz] =
        new faceZone
        (
            "rightExtrusionFaces",
            zone1,
            flipZone1,
            nFz,
            mesh().faceZones()
        );
    nFz++;

    fz[nFz] =
        new faceZone
        (
            "leftExtrusionFaces",
            zone2,
            flipZone2,
            nFz,
            mesh().faceZones()
        );
    nFz++;
    
    cz[nCz] = 
        new cellZone
        (
            "movingCells",
            zone3,
            0,
            mesh().cellZones()
        );
    nCz++;

    fz.setSize(nFz);
    cz.setSize(nCz);

    Info<< "Adding mesh zones." << endl;
    mesh().addZones(pz, fz, cz);

    // Add layer addition/removal interfaces

    List<polyMeshModifier*> tm(2);
    label nMods = 0;

    tm[nMods] =
        new layerAdditionRemoval
        (
            "right",
            nMods,
            topoChanger_,
            "rightExtrusionFaces",
            motionDict_.subDict("right").lookup<scalar>("minThickness"),
            motionDict_.subDict("right").lookup<scalar>("maxThickness")
        );
    nMods++;

    tm[nMods] = new layerAdditionRemoval
    (
        "left",
        nMods,
        topoChanger_,
        "leftExtrusionFaces",
        motionDict_.subDict("left").lookup<scalar>("minThickness"),
        motionDict_.subDict("left").lookup<scalar>("maxThickness")
    );
    nMods++;
    tm.setSize(nMods);

    Info<< "Adding " << nMods << " mesh modifiers" << endl;
    topoChanger_.addTopologyModifiers(tm);

    write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::movingCone::movingCone
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    fvMeshTopoChanger(mesh),
    topoChanger_(mesh),
    motionDict_(dict),
    motionVelAmplitude_(motionDict_.lookup("motionVelAmplitude")),
    motionVelPeriod_(motionDict_.lookup<scalar>("motionVelPeriod")),
    curMotionVel_
    (
        motionVelAmplitude_*sin(mesh.time().value()*pi/motionVelPeriod_)
    ),
    curLeft_(motionDict_.lookup<scalar>("leftObstacleEdge")),
    curRight_(motionDict_.lookup<scalar>("rightObstacleEdge")),
    curTop_(motionDict_.lookup<scalar>("upperObstacleEdge")),
    delta_(motionDict_.lookup<scalar>("minMeshWidth"))
{
    addZonesAndModifiers();

    curLeft_ = average
    (
        mesh.faceZones()
        [
            mesh.faceZones().findZoneID("leftExtrusionFaces")
        ]().localPoints()
    ).x() - small;

    curRight_ = average
    (
        mesh.faceZones()
        [
            mesh.faceZones().findZoneID("rightExtrusionFaces")
        ]().localPoints()
    ).x() + small;

    motionMask_ = vertexMarkup
    (
        mesh.points(),
        curLeft_,
        curRight_
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::movingCone::~movingCone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::movingCone::update()
{
    // Do mesh changes (use inflation - put new points in topoChangeMap)
    autoPtr<polyTopoChangeMap> topoChangeMap = topoChanger_.changeMesh(true);

    // Calculate the new point positions depending on whether the
    // topological change has happened or not
    pointField newPoints;

    vector curMotionVel_ =
        motionVelAmplitude_*sin(mesh().time().value()*pi/motionVelPeriod_);

    Pout<< "time:" << mesh().time().value()
        << " curMotionVel_:" << curMotionVel_
        << " curLeft:" << curLeft_
        << " curRight:" << curRight_
        << endl;

    if (topoChangeMap.valid())
    {
        Info<< "Topology change. Calculating motion points" << endl;

        if (topoChangeMap().hasMotionPoints())
        {
            Info<< "Topology change. Has premotion points" << endl;

            motionMask_ =
                vertexMarkup
                (
                    topoChangeMap().preMotionPoints(),
                    curLeft_,
                    curRight_
                );

            // Move points inside the motionMask
            newPoints =
                topoChangeMap().preMotionPoints()
              + motionMask_*curMotionVel_*mesh().time().deltaT().value();
        }
        else
        {
            Info<< "Topology change. Already set mesh points" << endl;

            motionMask_ =
                vertexMarkup
                (
                    mesh().points(),
                    curLeft_,
                    curRight_
                );

            // Move points inside the motionMask
            newPoints =
                mesh().points()
              + motionMask_*curMotionVel_*mesh().time().deltaT().value();
        }
    }
    else
    {
        Info<< "No topology change" << endl;
        // Set the mesh motion
        newPoints =
            mesh().points()
          + motionMask_*curMotionVel_*mesh().time().deltaT().value();
    }

    // The mesh now contains the cells with zero volume
    Info << "Executing mesh motion" << endl;
    mesh().movePoints(newPoints);

    //  The mesh now has got non-zero volume cells

    curLeft_ = average
    (
        mesh().faceZones()
        [
            mesh().faceZones().findZoneID("leftExtrusionFaces")
        ]().localPoints()
    ).x() - small;

    curRight_ = average
    (
        mesh().faceZones()
        [
            mesh().faceZones().findZoneID("rightExtrusionFaces")
        ]().localPoints()
    ).x() + small;

    return true;
}


void Foam::fvMeshTopoChangers::movingCone::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fvMeshTopoChangers::movingCone::mapMesh
(
    const polyMeshMap& map
)
{}


void Foam::fvMeshTopoChangers::movingCone::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
