/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    lamvisc

Description
    Calculates the laminar viscosity (mu) from the temperature field T

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "fluidThermo.H"
#include "singlePhaseTransportModel.H"
//#include "viscosityModel.H"
//#include "viscosityModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions(false, true);   // no -constant, with -zeroTime
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

/*
    autoPtr<fluidThermo> thermo
    (
        fluidThermo::New(mesh)
    );
*/

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();


        IOobject Theader
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
/*
        IOobject rhoHeader
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        );
*/
/*
        if (!Theader.headerOk())
        {
            Info<< "No T .. skipping" << endl;
            continue;
        }

        // no density, likely no solution
        if (!rhoHeader.headerOk())
        {
            Info<< "No rho ... skipping" << endl;
            continue;
        }
*/
        volScalarField T(Theader, mesh);


	Info<< "Reading field U\n" << endl;
	volVectorField U
	(
	    IOobject
	    (
		"U",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
	);

    Info<< "    Adding to phiFluid\n" << endl;
       surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(U)
        );
        

	volTensorField gradUx = fvc::grad(U);
//Info<< "gradUx before" << gradUx<< endl;	
	forAll(gradUx, i)
	{
//	Info<< "i" << i<< endl;	
		forAll(gradUx[i],j)
		{
//			Info<< "j" << gradUx[i][j]<< endl;
			if ((j == 3) or (j == 6))
			{continue;}
			else
			{gradUx[i][j] = 0;}
		}
	}
//Info<< "gradUx after" << gradUx<< endl;	

        volScalarField shearRate
        (
            IOobject
            (
                "StrainRate",
                runTime.timeName(),
                mesh
            ),
//            strainRate()
            Foam::sqrt(2.0)*mag(symm(gradUx))
 //           sqrt(2.0)
        );

        Info<< "Writing: " << " shearRate" << endl;
        shearRate.write();
    
    IOdictionary transportProperties
	(
	IOobject
	(
	    "transportProperties",
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE
	)
	);

	dimensionedScalar k0
	(
	"k0",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("k0")
	);
	    
	dimensionedScalar k1
	(
	"k1",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("k1")
	);
	        
    
 	dimensionedScalar k2
	(
	"k2",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("k2")
	);
	       
 	dimensionedScalar k_rho
	(
	"k_rho",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("k_rho")
	);
	
	dimensionedScalar n0
	(
	"n0",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("n0")
	);
	    
	dimensionedScalar n1
	(
	"n1",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("n1")
	);
	        
    
 	dimensionedScalar n2
	(
	"n2",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("n2")
	);
	
	
	 	dimensionedScalar nuMax
	(
	"nuMax",
	transportProperties.optionalSubDict("tempdeppowerLawPolynomialCoeffs").lookup("nuMax")
	);
	        
	       
            
  

            volScalarField viscosity_sym
        (
            IOobject
            (
                "nu_sym",
                runTime.timeName(),
                mesh
            ),
                        
            min(
            	nuMax,
            
               (k0+k1*(T-dimensionedScalar(dimTemperature, 273.15))+k2*pow((T-dimensionedScalar(dimTemperature, 273.15)),2))*k_rho*
                pow(
                max
                (
//                    dimensionedScalar(dimTime, 1.0)*strainRate(),
//			strainRate_new,
			dimensionedScalar(dimTime, 1.0)*shearRate,			

                    dimensionedScalar(dimless, small)
                ),	
                (n0+n1*T+n2*pow(T,2.0) - 1.0)

            )
            )
        );
//	                   Info<< "n0: " << n0 << endl;   
//	                   Info<< "n1: " << n1 << endl;   
//	                   Info<< "n2: " << n2 << endl;   	                   	                   
    viscosity_sym.write();
    
    
    
    
    
    
        singlePhaseTransportModel fluid(U, phi);
        volScalarField viscosity
        (
            IOobject
            (
                "nu",
                runTime.timeName(),
                mesh
            ),
            fluid.nu()
        );

        Info<< "Writing: " << " nu" << endl;
        viscosity.write();


//        Info<< "gradUx " << Foam::fvc::grad(U.component(vector::X))<< endl;
//         Info<< "gradUy " << Foam::fvc::grad(U.component(vector::Y))<< endl;       
//        Info<< "gradU " << Foam::fvc::grad(U)<< endl;        
//        Info<< "point " << U.mesh().C()<< endl;


    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
