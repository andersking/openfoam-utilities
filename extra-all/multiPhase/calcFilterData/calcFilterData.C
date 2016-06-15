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
    calcAlphaBoxVOF

Description
    Calculates the position and the force acting on a single droplet
    in VoF simulations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "sampledIsoSurface.H"
#include "interpolationCellPoint.H"
#include "vtkSurfaceWriter.H"
#include "boundBox.H"
    #include "argList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    
    timeSelector::addOptions();

    argList::addOption
    (
      "bb",
      "boundBox",
      "provide boxMin and boxMax coordinates '(0 0 0)' '(1 1 1)'"
    );

    argList::addOption
    (
      "box",
      "vector",
      "provide boxMin and boxMax coordinates '(0 0 0)' '(1 1 1)'"
    );

    argList::addOption
    (
      "xcoordinate",
      "scalar",
      "provide xMin and xMax coordinates '0 1 "
    );

    argList::addOption
    (
      "ycoordinate",
      "scalar",
      "provide yMin and yMax coordinates '0 1 "
    );

    argList::addOption
    (
      "zcoordinate",
      "scalar",
      "provide zMin and zMax coordinates '0 1 "
    );
    
    argList::addBoolOption("zeroMissing","use zeros for missing data");
    argList::addBoolOption("noOutputFile","Don't write a summary file");
    argList::addOption("outputFilename","<filename>","file to write data to");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    if (args.options().empty())
    {
       FatalErrorIn(args.executable())
	<< "No options supplied, please use the following format"
	   "-box (xmin ymin zmin)' '(xmax ymax zmax) or "
	   "-xcoordinate xmin xmax or "
	   "-ycoordinate ymin ymax or "
	   "-zcoordinate zmin zmax"
	<< exit(FatalError);
    }

    // read a bounding box
    // format ((xmin ymin zmin) (xmax ymax zmax))
    boundBox bb(vector(0,0,0),vector(1,1,1));
    if (args.optionFound("bb"))
    {
        bb = args.optionRead<boundBox>("bb");
    }
    // read a list of two vectors
    // (same as above, but not a bounding box)
    // "((xmin ymin zmin) (xmax ymax zmax))"
    FixedList<vector,2> box;
    if (args.optionFound("box"))
    {
        box = args.optionRead<FixedList<vector,2> >("box");
        Info << box << endl;
        bb.min() = box[0];
        bb.max() = box[1];
    }
    
    // read a list of 2 scalars, can also be a "tuple"
    // (xmin xmax)
    FixedList<scalar,2> minmax;
    if (args.optionFound("xcoordinate"))
    {
        minmax = args.optionRead<FixedList<scalar,2> >("xcoordinate");
        Info << minmax << endl;
        bb.min().x() = minmax[0];
        bb.max().x() = minmax[1];
    }

    if (args.optionFound("ycoordinate"))
    {
        minmax = args.optionRead<FixedList<scalar,2> >("ycoordinate");
        Info << minmax << endl;
        bb.min().y() = minmax[0];
        bb.max().y() = minmax[1];
    }

    if (args.optionFound("zcoordinate"))
    {
        minmax = args.optionRead<FixedList<scalar,2> >("zcoordinate");
        Info << minmax << endl;
        bb.min().z() = minmax[0];
        bb.max().z() = minmax[1];
    }

    Info << bb << endl;

    fileName outfileName("filterData.xy");

    bool zeroMissing = args.optionFound("zeroMissing");

    args.optionReadIfPresent("outputFilename",outfileName);
    bool noWriteFile = args.optionFound("noOutputFile");

    autoPtr<OFstream> outfile;

    if (!noWriteFile)
    {
        Info << "Writing to " << outfileName << "." << endl;
        outfile.set(new OFstream(outfileName));
        outfile() 
            << "#time " 
            << "vofFull " 
            << "avgVoFFull " 
            << "centreFull.x " 
            << "centreFull.y " 
            << "centreFull.z " 
            << "vofBox " 
            << "avgVoFBox " 
            << "centreBox.x " 
            << "centreBox.y " 
            << "centreBox.z " 
            << endl;
    }




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
            mesh.readUpdate();
            volScalarField alpha1(alpha1header, mesh);

            scalar sumAlpha1 = 0;
            scalar sumAlpha1Box = 0;
            scalar Vboxsum = 0;
            scalar Vsum = 0;
            vector centre = vector::zero;
            vector boxcentre = vector::zero;
            
            forAll(mesh.V(),celli)
            {
                scalar vcell = mesh.V()[celli];
                scalar acell = alpha1.internalField()[celli];
                vector ccell = mesh.C()[celli];
                sumAlpha1 += vcell*acell;
                Vsum += vcell;
                centre += vcell*acell*ccell;

                if (bb.contains(mesh.C()[celli])) 
                {
                    Vboxsum += vcell;
                    sumAlpha1Box += vcell*acell;
                    boxcentre += vcell*acell*ccell;
                }
            }
            reduce(sumAlpha1,sumOp<scalar>());
            reduce(Vsum,sumOp<scalar>());
            reduce(centre,sumOp<vector>());

            reduce(sumAlpha1Box,sumOp<scalar>());
            reduce(Vboxsum,sumOp<scalar>());
            reduce(boxcentre,sumOp<vector>());
            
            Info << "time = " << runTime.timeName() << " alpha1boxVOF = " << (sumAlpha1Box/Vboxsum) << endl;
            
            if (!noWriteFile)
            {
                outfile() 
                    << runTime.timeName() << " " 
                    << sumAlpha1 << " " 
                    << (sumAlpha1/Vboxsum) << " " 
                    << centre.x() << " "
                    << centre.y() << " "
                    << centre.z() << " "
                    // box components
                    << sumAlpha1Box << " " 
                    << (sumAlpha1Box/Vboxsum) << " "
                    << boxcentre.x() << " "
                    << boxcentre.y() << " "
                    << boxcentre.z() 
                << endl;
            }
        }
        else
        {
            if (zeroMissing && !noWriteFile)
            {
                outfile() 
                    << runTime.timeName() << " " 
                    << 0. << " " 
                    << 0. << " " 
                    << 0. 
                << endl;
            }
            Info<< "    No data" << endl;
        }

        Info<< endl;
    }

    return 0;
}


// ************************************************************************* //
