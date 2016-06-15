/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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

Application
    calcDropletForces

Description
    Calculates the position and the force acting on a single droplet
    in VoF simulations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sampledIsoSurface.H"
#include "interpolationCellPoint.H"
#include "vtkSurfaceWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

//#   include "createMesh.H"
//
// createMesh.H
// ~~~~~~~~~~~~

    Foam::Info
        << "Create mesh for time = "
        << runTime.timeName() << Foam::nl << Foam::endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    IOdictionary transportProperties
    (
        IOobject
        (   
            "transportProperties", 
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    word alphaName("alpha1");

    if (transportProperties.found("phases"))
    {
        word phase1Name = wordList(transportProperties.lookup("phases"))[0];
    
        alphaName = "alpha."+phase1Name;
    }

    
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject alpha1header
        (
            alphaName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check data exists
        if ( alpha1header.headerOk())
        {

            mesh.readUpdate(); //update the mesh

            Info<< "    Reading " << alphaName << endl;
            volScalarField alpha1(alpha1header, mesh);

            Info << "    Calculating droplet position (mean volume)" << endl; 

            vectorField aCV = (alpha1.internalField() * mesh.V() * mesh.C());
            vector centre = gSum(((alpha1.internalField()) * mesh.V() * mesh.C())) / gSum(((alpha1.internalField()) * mesh.V()));

            Info << "Droplet centre = " << centre << endl;
        }
        else
        {
            Info<< "    No data" << endl;
        }

        Info<< endl;
    }

    return 0;
}


// ************************************************************************* //
