from symplyphysics import (
    symbols, Eq, pretty, solve, Quantity, units,
    validate_input, validate_output, expr_to_quantity
)

# Description
## Field circulation along closed curve is curvilinear integral of field by this curve.
## C = CurveIntegral(F * dl, Curve), where
## C is circulation
## F is field vector function
## l if radius-vector of the curve
## In other words, circulation of the field along the closed curve is flow of the rotor of this field 
## through any area surrounded by this curve.
## CurveIntegral(F * dl, Curve) = Integral(Integral(Rot(F) * dS)), where
## S is area surrounded by Curve.
## Potential field is the field with zero rotor. Also potential field is called irrotational field.
## Work to move the object along the closed curve in the potential field is zero.