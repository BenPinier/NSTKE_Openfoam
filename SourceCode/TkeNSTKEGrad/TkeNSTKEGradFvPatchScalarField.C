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

#include "TkeNSTKEGradFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "turbulenceModel.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TkeNSTKEGradFvPatchScalarField::
TkeNSTKEGradFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),    
	Ustar_(p.size(),0.0),
	Cmu_(p.size(),0.0)
{}


Foam::TkeNSTKEGradFvPatchScalarField::
TkeNSTKEGradFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    Ustar_("Ustar",dict,p.size()),
    Cmu_("Cmu",dict,p.size())
{}


Foam::TkeNSTKEGradFvPatchScalarField::
TkeNSTKEGradFvPatchScalarField
(
    const TkeNSTKEGradFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
	Ustar_(ptf.Ustar_),
    Cmu_(ptf.Cmu_)
{}


Foam::TkeNSTKEGradFvPatchScalarField::
TkeNSTKEGradFvPatchScalarField
(
    const TkeNSTKEGradFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
	Ustar_(wbppsf.Ustar_),
    Cmu_(wbppsf.Cmu_)
{}


Foam::TkeNSTKEGradFvPatchScalarField::
TkeNSTKEGradFvPatchScalarField
(
    const TkeNSTKEGradFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
	Ustar_(wbppsf.Ustar_),
    Cmu_(wbppsf.Cmu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TkeNSTKEGradFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

	const scalarField& y = turbModel.y()[patchi];

    const scalarField k_internal(patchInternalField());

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

	scalarField nut(turbModel.nut(patch().index()));

	Info << nuw << nl;

	forAll(nuw,faceI){
			if (y[faceI] > 1)
				this->gradient()[faceI] = -0.00025*k_internal[faceI]*sqrt(k_internal[faceI])/(nuw[faceI]+nut[faceI]);
			else
				this->gradient()[faceI] = -0.00025*k_internal[faceI]*sqrt(k_internal[faceI])/(nuw[faceI]+nut[faceI]);


	}

	   	Info << average(-0.5*k_internal*sqrt(k_internal)/(nuw+nut)) << " " << average(k_internal) << nl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::TkeNSTKEGradFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    Ustar_.writeEntry("Ustar", os);
    this->writeEntry("value", os);
    Cmu_.writeEntry("Cmu", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        TkeNSTKEGradFvPatchScalarField
    );
}

// ************************************************************************* //
