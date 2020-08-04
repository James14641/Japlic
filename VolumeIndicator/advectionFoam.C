/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    advectionFoam

Description
    Solves a transport equation for a passive scalar using an explicit and/or
    implicit time-stepping method.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "velocityField.H"
#include "CourantNoFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("VOL","New Volume method");
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()

    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<label>("nCorr", label(2));
    const scalar CoLimit = readScalar(mesh.schemesDict().lookup("CoLimit"));

    #include "createFields.H"

    Info<< "Maximum total Courant number: " << max(Co).value() << endl;

    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << endl;
        if (args.options().found("VOL"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                
                + INDI*fvm::div(phi, T, "upwind")
                + INDE*fvc::div(phi, T, "explicit")
                );
                TEqn.solve();
            }
        }
        else
        {
            FatalErrorIn("JimpExpEulerFoam")
                 << " no valid advection scheme given as an argument"
                 << exit(FatalError);
        }
        
        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;
        runTime.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
