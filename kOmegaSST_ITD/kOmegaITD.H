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
    Foam::RASModels::kOmegaITD

Group
    grpRASTurbulence

Description
    Standard high Reynolds-number k-omega turbulence model for
    incompressible and compressible flows.

    *** Implemented also the turbulence damping formulation (Egorov 2004) for
    *** turbulence reduction on the interface.

    References:
    \verbatim
        Wilcox, D. C. (1998).
        Turbulence modeling for CFD
        (Vol. 2, pp. 103-217). La Canada, CA: DCW industries.
    \endverbatim

    The default model coefficients correspond to the following:
    \verbatim
        kOmegaITDCoeffs
        {
            Cmu         0.09;  // Equivalent to betaStar
            gamma       0.55556;
            beta        0.075;
            alphak      0.5;
            alphaOmega  0.5;

            // specific parameters for the turbulence damping Egorov model (mandatory, no predefined value - except for B)
            B             B   [0 0  0 0 0 0 0]  100;     // Damping factor, suggested value: 100
            dn            dn  [0 1  0 0 0 0 0]  0.05;    // cell height normal to interface (average) - may be computed but usually it is better to generate a mesh rather constant
            nu1           nu1 [0 2 -1 0 0 0 0]  1.2e-06; // viscosity of phase1 (the same in transportProperties dict)
            nu2           nu2 [0 2 -1 0 0 0 0]  5.0e-04; // viscosity of phase2 (the same in transportProperties dict)

            kAmb         kAmb [0 2 -2 0 0 0 0]  1e-05;   // kinetic energy ambient level (usually 1e-6*Uinf^2)
            omegaAmb omegaAmb [0 0 -1 0 0 0 0]  10;      // turbulent vorticity ambient level (usually 5 Uinf/L)

            vortProdForm  no;                            // TKE production based on vorticity formulation
        }
    \endverbatim

    *****UPDATE 31/08/2012

    revision of coefficients and TKE production limiter (as suggested in CFX
    solver for this turb. model, also yet implemented in the openfoam
    kw-SST turb. model).
    
    See:
    "Wilcox (1988) k-omega Two-Equation Model (Wilcox1988)" at the address:
    http://turbmodels.larc.nasa.gov/wilcox.html
    
    *****END UPDATE 31/08/2012
SourceFiles
    kOmegaITD.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaITD_H
#define kOmegaITD_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class kOmegaITD Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaITD
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar beta_;
            dimensionedScalar gamma_;
            dimensionedScalar alphaK_;
            dimensionedScalar alphaOmega_;

            dimensionedScalar B_;
            dimensionedScalar dn_;
            dimensionedScalar nu1_;
            dimensionedScalar nu2_;

            dimensionedScalar kAmb_;
            dimensionedScalar omegaAmb_;

            Switch vortProdForm_;


        // Fields

            volScalarField k_;
            volScalarField omega_;


    // Protected Member Functions

        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaITD");


    // Constructors

        //- Construct from components
        kOmegaITD
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaITD()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    alphaK_*this->nut_ + this->nu()
                )
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DomegaEff",
                    alphaOmega_*this->nut_ + this->nu()
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epsilon",
                        this->mesh_.time().timeName(),
                        this->mesh_
                    ),
                    Cmu_*k_*omega_,
                    omega_.boundaryField().types()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaITD.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
