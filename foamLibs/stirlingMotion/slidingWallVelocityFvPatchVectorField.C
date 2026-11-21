/*
Copyright (c) 2026 Zhi Sheng Chen

Licensed under the PolyForm Noncommercial License 1.0.0 (the "License").
You may not use this file except in compliance with the License. Any
commercial use requires a separate commercial license from the copyright
holder.

You may obtain a copy of the License at:
    https://polyformproject.org/licenses/noncommercial/1.0.0

Or see the LICENSE file in the root of this repository.

This software is provided "as is", without warranty of any kind, express or
implied. See the License for the specific language governing permissions
and limitations.
*/

#include "slidingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slidingWallVelocityFvPatchVectorField::slidingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::slidingWallVelocityFvPatchVectorField::slidingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false)
{
    // value is optional (zero = wall at rest)
    if (dict.found("value"))
    {
        fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(vector::zero);
    }
}


Foam::slidingWallVelocityFvPatchVectorField::slidingWallVelocityFvPatchVectorField
(
    const slidingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::slidingWallVelocityFvPatchVectorField::slidingWallVelocityFvPatchVectorField
(
    const slidingWallVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf)
{}


Foam::slidingWallVelocityFvPatchVectorField::slidingWallVelocityFvPatchVectorField
(
    const slidingWallVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::slidingWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvMesh& mesh = internalField().mesh();

    if (mesh.moving())
    {
        // volume swept by the face this step / face area = wall-normal velocity of the face
        const fvPatch& p = patch();
        const scalarField& phip = p.patchField<surfaceScalarField, scalar>(mesh.phi());
        const vectorField n(p.nf());
        vectorField::operator=(n*phip/(p.magSf() + VSMALL));
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        slidingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
