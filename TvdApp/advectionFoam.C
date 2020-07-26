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
    Foam::argList::addBoolOption("timeVaryingWind", "read the wind field (U/Uf/phi) at every timestep");
    Foam::argList::addBoolOption("FEBE","1st order explicit euler, or 1st order implicit euler where needed");
    Foam::argList::addBoolOption("FEBEDC","Deffered correction");
    Foam::argList::addBoolOption("FEBEHO","High space order implicit");
    Foam::argList::addBoolOption("SSP2BE","SSP2, or 1st order implicit euler where needed");
    Foam::argList::addBoolOption("CNDC","CNDCjt");
    Foam::argList::addBoolOption("SSP2CN","SSP2, or CrankNicholson where needed");
    Foam::argList::addBoolOption("SSP104BE","SSP104BE");
    Foam::argList::addBoolOption("SSP104BEDC","SSP104BEDC");
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()

    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<label>("nCorr", label(2));
    const scalar CoLimit = readScalar(mesh.schemesDict().lookup("CoLimit"));

    
    // Output file to write error measures each time step
    OFstream os("errorMeasures.dat");
    os << "#Time minT maxT sumT TV" << endl;
    ///////////////////////////////////////
    
    
    #include "createFields.H"

    Info<< "Maximum total Courant number: " << max(Co).value()
        << "\nMaximum explicit Courant number: " << max(CoExp).value()
        << "\nMaximum implicit Courant number: " << max(CoImp).value() << endl;
    Info<< "\nCalculating advection\n" << endl;
    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << endl;
        
        if (args.options().found("FEBE"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                
                + fvm::div(phiBig, T, "upwind")
                + fvc::div(phiSmall, T, "explicit")
                );
                TEqn.solve();
            }
        }
        else if (args.options().found("FEBEHO"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                  + fvm::div(phiBig, T, "implicit")
                  + fvc::div(phiSmall, T, "explicit")
                );
                TEqn.solve();
            }
        }
        else if (args.options().found("FEBEDC"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                + fvm::div(phiBig, T, "upwind")
                //      Correction term      //
                - fvc::div(phiBig, T, "upwind")
                + fvc::div(phiBig, T, "implicit")
               
                + fvc::div(phiSmall, T, "explicit")
                );
                TEqn.solve();
            }
        }
        
        else if (args.options().found("SSP2BE"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                KK = T + dt*fvc::div(phiSmall, T, "explicit");
            
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                + fvm::div(phiBig, T, "implicit")
                + 0.5*fvc::div(phiSmall, T, "explicit")
                + 0.5*fvc::div(phiSmall, KK, "explicit")
                );
                TEqn.solve();
            }
        }
        else if (args.options().found("SSP2CN"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                KK = T + dt*fvc::div(phiSmall, T, "explicit");
            
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                + 0.5*fvm::div(phiBig, T, "implicit")
                + 0.5*fvc::div(phiBig, T, "implicit")
                + 0.5*fvc::div(phiSmall, T, "explicit")
                + 0.5*fvc::div(phiSmall, KK, "explicit")
                );
                TEqn.solve();
            }
        }
        else if (args.options().found("CNDC"))
        {
            for(label corr = 0; corr < nCorr; corr++)
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                + 0.5*fvm::div(phiBig, T, "upwind")
                - 0.5*fvc::div(phiBig, T, "upwind")
                + 0.5*fvc::div(phiBig, T, "implicit")
                + 0.5*fvc::div(phiBig, T, "implicit")

                + 0.5*fvm::div(phiSmall, T, "upwind")
                - 0.5*fvc::div(phiSmall, T, "upwind")
                + 0.5*fvc::div(phiSmall, T, "implicit")
                + 0.5*fvc::div(phiSmall, T, "implicit")
                );
                TEqn.solve();
            }
        }
        
        else if (args.options().found("SSP104BE"))
        {
            // Explicit SSP104 update
            for (int i = 0; i < 4; i++)
            {
                T -= 1./6.*dt*fvc::div(phiSmall, T, "explicit");
            };
            KK = 3./5.*T.oldTime()
               + 2./5.*T
               - 1./15.*dt*fvc::div(phiSmall, T, "explicit");
            for (int i = 0; i < 4; i++)
            {
                KK -= 1./6.*dt*fvc::div(phiSmall, KK, "explicit");
            };
            
            T = 1./25.*T.oldTime()
              + 9./25.*T
              + 3./5.*KK 
              - 3./50.*dt*fvc::div(phiSmall, T, "explicit")
              - 1./10.*dt*fvc::div(phiSmall, KK, "explicit");
              
            // Implicit BE part (has to be used with 1st order upwind)
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phiBig, T, "upwind")
              - (T - T.oldTime())/dt
            );
            TEqn.solve();
        }
        
        
        else if (args.options().found("SSP104BEDC"))
        {
            // Explicit SSP104 update
            for (int i = 0; i < 4; i++)
            {
                T -= 1./6.*dt*fvc::div(phiSmall, T, "explicit");
            };
            KK = 3./5.*T.oldTime()
               + 2./5.*T
               - 1./15.*dt*fvc::div(phiSmall, T, "explicit");
            for (int i = 0; i < 4; i++)
            {
                KK -= 1./6.*dt*fvc::div(phiSmall, KK, "explicit");
            };
            
            T = 1./25.*T.oldTime()
              + 9./25.*T
              + 3./5.*KK 
              - 3./50.*dt*fvc::div(phiSmall, T, "explicit")
              - 1./10.*dt*fvc::div(phiSmall, KK, "explicit");
              
            // Lets try doing the defered correction. 
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phiBig, T, "upwind")
              - fvc::div(phiBig, T, "upwind")
              + fvc::div(phiBig, T, "implicit")
              - (T - T.oldTime())/dt
            );
            TEqn.solve();
        }
        else
        {
            FatalErrorIn("JimpExpEulerFoam")
                 << " no valid advection scheme given as an argument"
                 << exit(FatalError);
        }
        
        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;
        ///////////////////////
        scalarField TV = 0.25*fvc::surfaceSum
        (
            mag(fvc::snGrad(T))/mesh.deltaCoeffs()
        )().primitiveField();
        const scalarField& V = mesh.V().field();

        os << runTime.timeName() << " " << min(T.internalField()).value() << " "
           << max(T.internalField()).value() << " "
           << gSum(T.primitiveField()*V)/gSum(V) << " "
           << gSum(TV) << endl;
        ///////////////////////
        runTime.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
