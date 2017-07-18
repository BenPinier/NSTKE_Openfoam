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
    Foam::incompressible::RASModels::PrandltkEqn

Group
    grpIcoRASTurbulence

Description
    Low Reynolds-number k-kl-omega turbulence model for
    incompressible flows.

    This turbulence model is described in:
    \verbatim
        Walters, D. K., & Cokljat, D. (2008).
        A three-equation eddy-viscosity model for Reynolds-averaged
        Navier–Stokes simulations of transitional flow.
        Journal of Fluids Engineering, 130(12), 121401.
    \endverbatim

    however the paper contains several errors which must be corrected for the
    model to operation correctly as explained in

    \verbatim
        Furst, J. (2013).
        Numerical simulation of transitional flows with laminar kinetic energy.
        Engineering MECHANICS, 20(5), 379-388.
    \endverbatim

    All these corrections and updates are included in this implementation.

    The default model coefficients are
    \verbatim
        PrandltkEqnCoeffs
        {
            A0             4.04
            As             2.12
            Av             6.75
            Abp            0.6
            Anat           200
            Ats            200
            CbpCrit        1.2
            Cnc            0.1
            CnatCrit       1250
            Cint           0.75
            CtsCrit        1000
            CrNat          0.02
            C11            3.4e-6
            C12            1.0e-10
            CR             0.12
            CalphaTheta    0.035
            Css            1.5
            CtauL          4360
            Cw1            0.44
            Cw2            0.92
            Cw3            0.3
            CwR            1.5
            Clambda        2.495
            CmuStd         0.09
            Prtheta        0.85
            Sigmak         1
            Sigmaw         1.17
        }
    \endverbatim

SourceFiles
    PrandltkEqn.C

\*---------------------------------------------------------------------------*/

#ifndef PrandltkEqn_H
#define PrandltkEqn_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class PrandltkEqn Declaration
\*---------------------------------------------------------------------------*/

class PrandltkEqn
:
    public eddyViscosity<incompressible::RASModel>
{
    // Private member functions




protected:

 
  // Protected data

        // Model coefficients

            dimensionedScalar Ck_;

            dimensionedScalar Cmu_;
            dimensionedScalar Ustar_;


			dimensionedScalar Cdelta_0;
			dimensionedScalar Cdelta_1;
			dimensionedScalar Cdelta_2;
			dimensionedScalar Cdelta_3;
			dimensionedScalar Cdelta_4;
			dimensionedScalar Cdelta_5;
			dimensionedScalar Cdelta_6;
			dimensionedScalar y0_;
        // Fields

            volScalarField k_;
            volScalarField epsilon_;



			volScalarField MixingLength_;

    // Protected Member Functions



        virtual void correctNut();




public:

    //- Runtime type information
    TypeName("PrandltkEqn");


    // Constructors

        //- Construct from components
        PrandltkEqn
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~PrandltkEqn()
    {}


        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_+ nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
