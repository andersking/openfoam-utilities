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
    initTBL

Description
    Apply a simplified boundary-layer model to the velocity and
    turbulence fields based on the 1/7th power-law.

    The uniform boundary-layer thickness is either provided via the -ybl option
    or calculated as the average of the distance to the wall scaled with
    the thickness coefficient supplied via the option -Cbl.  If both options
    are provided -ybl is used.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "singlePhaseTransportModel.H"
//#include "RASModel.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// turbulence constants - file-scope
//static const scalar Cmu(0.09);
//static const scalar kappa(0.41);


int main(int argc, char *argv[])
{

    argList::validArgs.append("Uprime");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Random rand(1);

    Info<< "Time = " << runTime.timeName() << nl << endl;

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );



    scalar maguprime(args.argRead<scalar>(1));
    
    forAll(U.internalField(),cI)
    {
        vector  uprime((2*rand.vector01() - vector::one)*maguprime);

        U.internalField()[cI] += uprime;
        
    }


    Info<< "Writing U\n" << endl;
    U.write();

    // Update/re-write phi
//    phi = fvc::interpolate(U) & mesh.Sf();
//    phi.write();

    // Calculate nut
//    tmp<volScalarField> tnut = turbulence->nut();
//    volScalarField& nut = tnut();
//    volScalarField S(mag(dev(symm(fvc::grad(U)))));
//    nut = sqr(kappa*min(y, ybl))*::sqrt(2)*S;
//
//    if (args.optionFound("writenut"))
 //   {
  //      Info<< "Writing nut" << endl;
//        nut.write();
//    }

    // Create G field - used by RAS wall functions
//    volScalarField G("RASModel::G", nut*2*sqr(S));


    //--- Read and modify turbulence fields

    // Turbulence k
//    tmp<volScalarField> tk = turbulence->k();
//    volScalarField& k = tk();
//    scalar ck0 = pow025(Cmu)*kappa;
//    k = sqr(nut/(ck0*min(y, ybl)));
//    k.correctBoundaryConditions();

//    Info<< "Writing k\n" << endl;
//    k.write();


    // Turbulence epsilon
//    tmp<volScalarField> tepsilon = turbulence->epsilon();
//    volScalarField& epsilon = tepsilon();
//    scalar ce0 = ::pow(Cmu, 0.75)/kappa;
//    epsilon = ce0*k*sqrt(k)/min(y, ybl);
//    epsilon.correctBoundaryConditions();

//    Info<< "Writing epsilon\n" << endl;
//    epsilon.write();


    // Turbulence omega
//    IOobject omegaHeader
 //   (
  //      "omega",
   //     runTime.timeName(),
    //    mesh,
     //   IOobject::MUST_READ,
      //  IOobject::NO_WRITE,
      //  false
 //   );

 //   if (omegaHeader.headerOk())
 //   {
//        volScalarField omega(omegaHeader, mesh);
 //       omega = epsilon/(Cmu*k);
  //      omega.correctBoundaryConditions();
//
 //       Info<< "Writing omega\n" << endl;
  //      omega.write();
  //  }


    // Turbulence nuTilda
//    IOobject nuTildaHeader
//    (
//        "nuTilda",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::NO_WRITE,
//        false
//    );

//    if (nuTildaHeader.headerOk())
//    {
//        volScalarField nuTilda(nuTildaHeader, mesh);
//        nuTilda = nut;
 //       nuTilda.correctBoundaryConditions();

//        Info<< "Writing nuTilda\n" << endl;
//        nuTilda.write();
//    }


    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
