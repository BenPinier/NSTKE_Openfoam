/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::fixedVelocityFunctionFvPatchVectorField

Group
    grpTurbulenceBoundaryConditions

Description
    Set a constant shear stress as tau0 = -nuEff dU/dn.

SourceFiles
    fixedVelocityFunctionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedVelocityFunctionFvPatchVectorField_H
#define fixedVelocityFunctionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class fixedVelocityFunctionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedVelocityFunctionFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Constant shear stress
        const vector tau0_;

		const scalar alpha_;



		 scalar y0_;
public:

    //- Runtime type information
    TypeName("fixedVelocityFunction");


    // Constructors

        //- Construct from patch and internal field
        fixedVelocityFunctionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedVelocityFunctionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        fixedVelocityFunctionFvPatchVectorField
        (
            const fixedVelocityFunctionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedVelocityFunctionFvPatchVectorField
        (
            const fixedVelocityFunctionFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedVelocityFunctionFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedVelocityFunctionFvPatchVectorField
        (
            const fixedVelocityFunctionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new fixedVelocityFunctionFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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
