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
#include "CourandNoFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::argList::addBoolOption("timeVaryingWind", "read the wind field (U/Uf/phi) at every timestep");
    Foam::argList::addBoolOption("FEBE","1st order explicit euler, or 1st order implicit euler where needed");
    Foam::argList::addBoolOption("FEBEDC","Deffered correction");
    Foam::argList::addBoolOption("FEBEHO","High space order implicit");
    Foam::argList::addBoolOption("SSP2BE","SSP2, or 1st order implicit euler where needed");
    Foam::argList::addBoolOption("SSP2CN","SSP2, or CrankNicholson where needed");
    Foam::argList::addBoolOption("SSP104BE","SSP104, or backward euler where needed");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #define dt runTime.deltaT()

    // Read the number of iterations each time-step
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nCorr = itsDict.lookupOrDefault<label>("nCorr", label(2));
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));
    const scalar CoLimit = readScalar(mesh.schemesDict().lookup("CoLimit"));

    #include "createFields.H"

    Info<< "\nCalculating advection\n" << endl;

    IOdictionary dict
    (
        IOobject
        (
            "advectionDict",
            mesh.time().system(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    bool timeVaryingWind = dict.lookupOrDefault<bool>("timeVaryingWind", false);
    const dictionary& velocityDict = dict.subOrEmptyDict("velocity");
    autoPtr<velocityField> v;
    if (velocityDict.size() > 0)
    {
        v = velocityField::New(dict.subOrEmptyDict("velocity"));
    }
    
    while (runTime.loop())
    {
        Info<< "\nTime = " << runTime.timeName() << endl;
        if (timeVaryingWind)
        {
            v->applyTo(phi);

            forAll(phi, faceI)
            {
                if (mag(phi[faceI]) > phiLimit[faceI])
                {
                    phiSmall[faceI] = 0;
                    phiBig[faceI] = phi[faceI];
                }
                else
                {
                    phiSmall[faceI] = phi[faceI];
                    phiBig[faceI] = 0;
                }
            }

            U = fvc::reconstruct(phi);
            Uf = linearInterpolate(U);
            Uf += (phi - (Uf & mesh.Sf()))*mesh.Sf()/sqr(mesh.magSf());
        }

        #include "CourantNo.H"
        
        
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
        if (args.options().found("FEBEHO"))
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
        if (args.options().found("FEBEDC"))
        {
            for(label corr = 0; corr < nCorr; corr++)// this scheme is very sensitive to the number of corrective passes. 
            {
                fvScalarMatrix TEqn
                (
                    fvm::ddt(T)
                + fvm::div(phiBig, T, "upwind")
                //   Correction term   //
                - fvc::div(phiBig, T, "upwind")
                + fvc::div(phiBig, T, "implicit")
               
                + fvc::div(phiSmall, T, "upwind")
                );
                TEqn.solve();
            }
        }
        
        if (args.options().found("SSP2BE"))
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
        if (args.options().found("SSP2CN"))
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
        
        if (args.options().found("SSP104BE"))// not functional
        {
            for(label corr = 0; corr < nCorr; corr++)
            {   
                
                cout << "Currently not functional" << endl;
                KK = T;
                for (int i = 1; i < 5; i++)
                {
                KK = KK + 1.0/6.0*dt*fvc::div(phiSmall, KK, "explicit");
                };
                KK2 = 3.0/5.0*T +2.0/5.0*KK + 1.0/15.0*dt*fvc::div(phiSmall, KK, "explicit");
                for (int i = 1; i < 5; i++)
                {
                KK2 = KK2 +1.0/6.0*dt*fvc::div(phiSmall, KK, "explicit");
                };
                
                fvScalarMatrix TEqn
                (
                  fvm::ddt(T)  
                + fvm::div(phiBig, T, "upwind")
                
                + -24.0/25.0*T*1.0/dt
                + 9.0/25.0*KK*1.0/dt
                + 3.0/5.0*KK2*1.0/dt
                + 3.0/50.0*fvc::div(phiSmall, KK, "explicit")
                + 1.0/10.0*fvc::div(phiSmall, KK2, "explicit")
                );
                TEqn.solve();
            }
        }
        
        Info << " T goes from " << min(T.internalField()).value() << " to "
             << max(T.internalField()).value() << endl;
        runTime.write();
        
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
