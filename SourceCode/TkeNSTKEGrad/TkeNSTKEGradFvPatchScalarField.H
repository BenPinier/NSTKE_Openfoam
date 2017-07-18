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

Class
    Foam::TkeNSTKEGradFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a free-stream condition for pressure.
    It is a zero-gradient condition that constrains the flux across the patch
    based on the free-stream velocity.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | none
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            TkeNSTKEGrad;
    }
    \endverbatim

Note
    This condition is designed to operate with a freestream velocity condition

See also
    Foam::zeroGradientFvPatchField
    Foam::freestreamFvPatchField

SourceFiles
    TkeNSTKEGradFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef TkeNSTKEGradFvPatchScalarFields_H
#define TkeNSTKEGradFvPatchScalarFields_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class TkeNSTKEGradFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class TkeNSTKEGradFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        scalarField	Ustar_;

		scalarField Cmu_;


public:

    //- Runtime type information
    TypeName("TkeNSTKEGrad");


    // Constructors

        //- Construct from patch and internal field
        TkeNSTKEGradFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        TkeNSTKEGradFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given TkeNSTKEGradFvPatchScalarField onto
        //  a new patch
        TkeNSTKEGradFvPatchScalarField
        (
            const TkeNSTKEGradFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        TkeNSTKEGradFvPatchScalarField
        (
            const TkeNSTKEGradFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new TkeNSTKEGradFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        TkeNSTKEGradFvPatchScalarField
        (
            const TkeNSTKEGradFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new TkeNSTKEGradFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
