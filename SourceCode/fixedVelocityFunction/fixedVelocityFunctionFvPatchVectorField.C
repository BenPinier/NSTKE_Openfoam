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

#include "fixedVelocityFunctionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedVelocityFunctionFvPatchVectorField::fixedVelocityFunctionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    tau0_(Zero),
	alpha_(0.0001)
{}


Foam::fixedVelocityFunctionFvPatchVectorField::fixedVelocityFunctionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    tau0_(dict.lookupOrDefault<vector>("tau", Zero)),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.0001)),
    y0_(dict.lookupOrDefault<scalar>("y0", 0.0001))
{
    fvPatchField<vector>::operator=(patchInternalField());
}


Foam::fixedVelocityFunctionFvPatchVectorField::fixedVelocityFunctionFvPatchVectorField
(
    const fixedVelocityFunctionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    tau0_(ptf.tau0_),
	alpha_(ptf.alpha_),
	y0_(ptf.y0_)
{}


Foam::fixedVelocityFunctionFvPatchVectorField::fixedVelocityFunctionFvPatchVectorField
(
    const fixedVelocityFunctionFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    tau0_(ptf.tau0_),
	alpha_(ptf.alpha_),
	y0_(ptf.y0_)
{}


Foam::fixedVelocityFunctionFvPatchVectorField::fixedVelocityFunctionFvPatchVectorField
(
    const fixedVelocityFunctionFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    tau0_(ptf.tau0_),
	alpha_(ptf.alpha_),
	y0_(ptf.y0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedVelocityFunctionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

	const word wordi = patch().name();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    scalarField nuEff(turbModel.nuEff(patch().index()));


    vectorField Uc(patchInternalField());

    //vector tauHat = tau0_/(mag(tau0_) + ROOTVSMALL);

    vector tauHat(1,0,0);
    vector e1(1,0,0);
    vector e3(0,0,1);
    const scalarField& ry = patch().deltaCoeffs();


	vectorField U_0(patchInternalField());
	vectorField U_C(patchInternalField());
	const scalar kappa = 0.41;
	vectorField U_g(patchInternalField());

    const tmp<scalarField> tnuw = turbModel.nu(patchi);

    const scalarField& nuw = tnuw();
    tmp<scalarField> tnutw(new scalarField(y));
    scalarField& nutw = tnutw.ref();

	scalarField nut(turbModel.nut(patch().index()));

	forAll(U_0,faceI){

	    U_0[faceI][0] = tau0_[0]*(1./kappa*log(y0_*tau0_[0]/nuw[faceI])+5.5);
		U_0[faceI][1] = 0.;
		U_0[faceI][2] = 0.;	

 
		
	}


	Info << wordi << " " << average(tauHat*(tauHat & (-tau0_*(1./(ry*(nut+nuw))) + Uc))  )<< nl;
 
 operator==(tauHat*(tauHat & (alpha_*mag(U_0-U_C)*(e1 * (e1 & (U_0-U_C)))*(1.0/(ry*(nut+nuw))) + Uc)) ) ;
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::fixedVelocityFunctionFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
    os.writeKeyword("tau") << tau0_ << token::END_STATEMENT << nl;
    os.writeKeyword("alpha") << alpha_ << token::END_STATEMENT << nl;
    os.writeKeyword("y0") << y0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        fixedVelocityFunctionFvPatchVectorField
    );
}

// ************************************************************************* //
