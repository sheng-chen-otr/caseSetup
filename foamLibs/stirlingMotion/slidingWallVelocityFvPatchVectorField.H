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

// static wall in a moving mesh (ESI OpenFOAM v2412): its points slide along the wall with the mesh, so its faces
// sweep a little volume (snapped faces are not exactly parallel to the motion). noSlip (U = 0) then lets gas through
// the wall, movingWallVelocity closes it but drags the gas with the moving points. Here the velocity is normal to
// the face only and matches the mesh flux: U = n meshPhi/|Sf|, no tangential velocity (zero when the mesh is static)

#ifndef slidingWallVelocityFvPatchVectorField_H
#define slidingWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class slidingWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    TypeName("slidingWallVelocity");

    slidingWallVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&
    );

    slidingWallVelocityFvPatchVectorField
    (
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&
    );

    slidingWallVelocityFvPatchVectorField
    (
        const slidingWallVelocityFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&
    );

    slidingWallVelocityFvPatchVectorField
    (
        const slidingWallVelocityFvPatchVectorField&
    );

    slidingWallVelocityFvPatchVectorField
    (
        const slidingWallVelocityFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&
    );

    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>
        (
            new slidingWallVelocityFvPatchVectorField(*this)
        );
    }

    virtual tmp<fvPatchVectorField> clone
    (
        const DimensionedField<vector, volMesh>& iF
    ) const
    {
        return tmp<fvPatchVectorField>
        (
            new slidingWallVelocityFvPatchVectorField(*this, iF)
        );
    }

    virtual void updateCoeffs();
};

} // End namespace Foam

#endif
