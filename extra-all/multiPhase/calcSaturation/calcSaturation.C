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
    calcSaturation

Description
    Calculates the layered saturation for combined particle/VoF filtration

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sampledIsoSurface.H"
#include "interpolationCellPoint.H"
#include "vtkSurfaceWriter.H"
#include "passiveParticle.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    bool gnuplotFormat = true; // plot each time as unbroken group, double spaced at end. First value is time

    argList::validArgs.append("nSections"); 

// TODO items
//    argList::addBoolOption
//    (
//        "tagDroplets",
//        "tag contiguous droplets - not implemented (yet!)"
//    );
//   argList::validArgs.append("Udir");// define flow direction - not implemented 
//   argList::validArgs.append("deadLength");// length of mesh in flow-dir not included as filter - not implemented 



#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

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

    Info << "Using " << alphaName << " for alphaName." << endl;

    const boundBox &bounds(mesh.bounds());
    
    label nSections = args.argRead<label>(1);
    
    //total volume
    scalar emptyVol = bounds.volume();
    scalar splitLength = bounds.span().z()/nSections;
    scalar splitVol=bounds.span().x()*bounds.span().y()*splitLength;
    
    autoPtr<OFstream> outfileAlphaGP;
    autoPtr<OFstream> outfileParticlesGP;
    autoPtr<OFstream> outfileTotalGP;
    
    if (gnuplotFormat)
    {
        outfileAlphaGP.reset(new OFstream("saturation-alpha-gnuplot.gxy"));
        outfileParticlesGP.reset(new OFstream("saturation-particles-gnuplot.gxy"));
        outfileTotalGP.reset(new OFstream("saturation-total-gnuplot.gxy"));
    }
    
    OFstream outfileAlpha("saturation-alpha.gxy");
    OFstream outfileParticles("saturation-particles.gxy");
    OFstream outfileTotal("saturation-total.gxy");
    
    bool firstWrite = true;

    // Time averages
    scalar lastTime = -1;
    List<scalar> avgAlpha(nSections,0.);
    List<scalar> avgParticles(nSections,0.);
    List<scalar> avgTotal(nSections,0.);
    
    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;
        if (firstWrite)
        {   
            firstWrite = false;
        } else {
            if (gnuplotFormat)
            {   
                outfileAlphaGP() << nl << nl << endl;
                outfileParticlesGP() << nl << nl << endl;
                outfileTotalGP() << nl << nl << endl;
            }    
            outfileAlpha << endl;
            outfileParticles << endl;
            outfileTotal << endl;
        }
        
        if (gnuplotFormat)
        {
            outfileAlphaGP() << runTime.timeName();
            outfileParticlesGP() << runTime.timeName();
            outfileTotalGP() << runTime.timeName();
        }
        
        outfileAlpha << runTime.timeName();
        outfileParticles << runTime.timeName();
        outfileTotal << runTime.timeName();

        IOobject alpha1header
        (
            alphaName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        List<scalar> splitAlpha(nSections,0.);
        List<scalar> splitParticles(nSections,0.);
        List<scalar> splitTotal(nSections,0.);

        // Check data exists
        if ( alpha1header.headerOk())
        {

            mesh.readUpdate(); //update the mesh

            Info<< "    Reading " << alphaName << endl;
            volScalarField alpha1(alpha1header, mesh);
            label split = 0;

            forAll(mesh.C(),cI)
            {
                split = (mesh.C()[cI].z()-bounds.min().z())/splitLength; //which partition
                splitAlpha[split] += alpha1.internalField()[cI]*mesh.V()[cI]; //volume of fluid
            }

        }
        else
        {
            Info<< "    No data" << endl;
        }

        // particles, read without reading data
        Cloud<passiveParticle> particles(mesh,"defaultCloud",false);
        
        Pout << particles.size() << " particles." <<endl;

        label split = 0;
        // header for diameter field
        
        IOobject dHeader
        (
            "d",
            mesh.time().timeName(),
            "lagrangian/defaultCloud",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );
        
        if (dHeader.headerOk())
        {

            IOField<scalar> d(dHeader);

            label pI = -1;
         
            forAllConstIter(Cloud<passiveParticle>, particles, p)
            {
                pI++;
                split = (p().position().z()-bounds.min().z())/splitLength;
                splitParticles[split] += M_PI*d[pI]*d[pI]*d[pI]/6.; 
            }
        }
    
        scalar splitI = 0;
    
        forAll(splitAlpha,sI)
        {
            splitI = splitAlpha[sI];
            Foam::reduce(splitI, sumOp<scalar>());
            splitAlpha[sI] = splitI;
        }

        forAll(splitParticles,sI)
        {
            splitI = splitParticles[sI];
            Foam::reduce(splitI, sumOp<scalar>());
            splitParticles[sI] = splitI;
        }

        forAll(splitTotal,sI)
        {
            splitTotal[sI] = splitAlpha[sI] + splitParticles[sI];
        }
        // Calc Averages
        if (lastTime < 0.)
        {
            forAll(avgAlpha,sI)
            {
                avgAlpha[sI] = splitAlpha[sI];
                avgParticles[sI] = splitParticles[sI];
                avgTotal[sI] = splitTotal[sI];
            }
            lastTime = runTime.value();
        } else {
            scalar currTime = runTime.value();
            forAll(avgAlpha,sI)
            {
                avgAlpha[sI] += (splitAlpha[sI] - avgAlpha[sI])*(1. - (lastTime/currTime));
                avgParticles[sI] += (splitParticles[sI] - avgParticles[sI])*(1. - (lastTime/currTime));
                avgTotal[sI] += (splitTotal[sI] - avgTotal[sI])*(1. - (lastTime/currTime));
            }
            lastTime = currTime;
            
        }
            
        
        if (gnuplotFormat)
        {
            forAll(splitAlpha,sI)
            {
                outfileAlphaGP() << nl << (splitAlpha[sI]/splitVol);
                outfileParticlesGP() << nl << (splitParticles[sI]/splitVol);
                outfileTotalGP() << nl << (splitTotal[sI]/splitVol);
            }
        } 

        forAll(splitAlpha,sI)
        {
            outfileAlpha << " " << (splitAlpha[sI]/splitVol);
            outfileParticles << " " << (splitParticles[sI]/splitVol);
            outfileTotal << " " << (splitTotal[sI]/splitVol);
        }
    

        forAll(avgAlpha,sI)
        {
            outfileAlpha << " " << (avgAlpha[sI]/splitVol);
            outfileParticles << " " << (avgParticles[sI]/splitVol);
            outfileTotal << " " << (avgTotal[sI]/splitVol);
        }
    
        
        Info<< endl;
    }

    return 0;
}


// ************************************************************************* //
