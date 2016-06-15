/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    uprime

Description
    Calculates and writes the scalar field of uprime (sqrt(2/3 k)).

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "fvc.H"
//#include "interfaceProperties.H"
//#include "twoPhaseMixture.H"
#include "typeInfo.H"
#include "oilParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

//    #include "addRegionOption.H"

//    argList::addBoolOption
//    (
//        "compressible",
//        "calculate compressible y+"
//    );

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
//    #include "createNamedMesh.H"
    #include "createMesh.H"

//    const bool compressible = args.optionFound("compressible");

//    bool writeResults = !args.optionFound("noWrite");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        fvMesh::readUpdateState state = mesh.readUpdate();

        IOobject alpha1header
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (alpha1header.headerOk())
        {
            Info<< "    Reading alpha1" << endl;
            volScalarField alpha1(alpha1header, mesh);

    Info<< "\nCreating particle cloud" << endl;
    oilParticleCloud parts(mesh); //,"defaultCloud",false);

            parts.readFields();

            label nParts = 0;

            nParts = parts.size();

            reduce(nParts,sumOp<label>());

            Info << nParts << " particles in cloud" << endl;
         
            scalar partVol = 0.;

            forAllIter(Cloud<oilParticle>, parts, pIter)
            {
                oilParticle& p = pIter();
                partVol += (M_PI/6.) * p.d()  * p.d()  * p.d();
            }

            reduce(partVol,sumOp<scalar>());

//        Info<< "    Calculating uprime" << endl;
//        volScalarField uprime
//        (
//            IOobject
//            (
//                "uprime",
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ
//            ),
//            sqrt((2.0/3.0)*k)
//        );

            Info<< "  Mesh Volume (actual)     : " << gSum(mesh.V()) << endl;
            Info<< "  Fluid Volume (liquid)    : " << gSum(mesh.V() * alpha1.internalField()) << endl;
            Info<< "  Fluid Volume (particles) : " << partVol << endl;
            Info<< "  Saturation               : " << (gSum(mesh.V() * alpha1.internalField()) + partVol )/ gSum(mesh.V()) << endl;

//        if (writeResults)
//        {
//            uprime.write();
//        }
        }
        else
        {
            Info<< "    No alpha1" << endl;
        }
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
