/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "PrandltkEqn.H"
#include "bound.H"
#include "wallDist.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PrandltkEqn, 0);
addToRunTimeSelectionTable(RASModel, PrandltkEqn, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void PrandltkEqn::correctNut()
{

	volVectorField C = mesh().C();
	volScalarField Y = C & vector(0,1,0);

	Info << "Correct nut" << nl;

//    nut_ =Ck_*sqrt(k_)*MixingLength_*Y;
// 	   nut_.correctBoundaryConditions();

	Info << "End Correct nut" << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PrandltkEqn::PrandltkEqn
(
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName,
            const word& type 
)
:    eddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ck",
            coeffDict_,
			dimless,  
            1.44
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
			dimless,  
            1.92
        )
    ),
    Ustar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ustar",
            coeffDict_,
			dimVelocity,     	       	
			1.0			
        )
    ),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    MixingLength_
    (
        IOobject
        (
            "MixingLength",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,        
        dimensionedScalar("MixingLength",  dimensionSet(0,1,0,0,0,0,0), 0.0)
    ),
    Cdelta_0
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdelta0",
            coeffDict_,
			dimLength,  
            0
        )
    ),
    Cdelta_1
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdelta1",
            coeffDict_,
			dimLength,  
            0
        )
    ),
    Cdelta_2
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdelta2",
            coeffDict_,
			dimLength,  
            0
        )
    ),
    Cdelta_3
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdelta3",
            coeffDict_,
			dimLength,  
            0
        )
    ),
    Cdelta_4
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdelta4",
            coeffDict_,
			dimLength,  
            0
        )
    ),
    Cdelta_5
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdelta5",
            coeffDict_,
			dimLength,  
            0
        )
    ),
    Cdelta_6
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cdelta6",
            coeffDict_,
			dimLength,  
            0
        )
    ),
	y0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "y0",
            coeffDict_,
			dimLength,  
            0
        )
    )
{
    bound(k_, kMin_);

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool PrandltkEqn::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        Ck_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
        Ustar_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

 

void PrandltkEqn::correct()
{
    if (!turbulence_)
    {
        return;
    }


    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
	volVectorField C = mesh().C();  

 
	volScalarField Y = wallDist::New(mesh_).y()+ y0_;


	scalar c_y;

	forAll(Y,i){
		if(Y[i] > 2.0)
			Y[i] = 2.0-Y[i] ;
	}
/*
	forAll(MixingLength_,i){
		c_y = Y[i];
		MixingLength_[i] = Cdelta_0 ;//+  Cdelta_1*c_y +  Cdelta_2*c_y*c_y +  Cdelta_3*c_y*c_y*c_y+  Cdelta_4*c_y*c_y*c_y*c_y +  Cdelta_5*c_y*c_y*c_y*c_y*c_y   + Cdelta_6*c_y*c_y*c_y*c_y*c_y*c_y;
	}
*/
  	dimensionedScalar c1("c1", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1);
 	dimensionedScalar c2("c2", dimensionSet(0, -2, 0, 0, 0, 0, 0), 1);
 	dimensionedScalar c3("c3", dimensionSet(0, -3, 0, 0, 0, 0, 0), 1);
 	dimensionedScalar c4("c4", dimensionSet(0, -4, 0, 0, 0, 0, 0), 1);
 	dimensionedScalar c5("c5", dimensionSet(0, -5, 0, 0, 0, 0, 0), 1);
 	dimensionedScalar c6("c6", dimensionSet(0, -6, 0, 0, 0, 0, 0), 1);

	MixingLength_ = Cdelta_0  +  Cdelta_1*Y*c1 +  Cdelta_2*Y*Y*c2 +  Cdelta_3*Y*Y*Y*c3 +  Cdelta_4*Y*Y*Y*Y*c4 +  Cdelta_5*Y*Y*Y*Y*Y*c5   + Cdelta_6*Y*Y*Y*Y*Y*Y*c6;
	MixingLength_.correctBoundaryConditions();


 	dimensionedScalar c001("c2", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0);

	volScalarField mutk = Cmu_*MixingLength_*sqrt(k_)+ nu();


	nut_ = Ck_*sqrt(k_)*MixingLength_;
	nut_.correctBoundaryConditions();





	volScalarField G(GName(), 2.0*nut*magSqr(symm(fvc::grad(U))));

    tmp<fvScalarMatrix> kEqn
    (
       fvm::ddt(k_)
     + fvm::div(phi_, k_)
     - fvm::laplacian(mutk, k_)
    ==
       G
     - fvm::Sp(sqrt(k_)/(MixingLength_+ c001), k_)
    );

    kEqn.ref().relax();
   	solve(kEqn);

    bound(k_, kMin_);


//    PrandltkEqn.ref().relax();

 //   solve(PrandltkEqn);

    bound(k_, this->kMin_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
