from sympy import Eq
from symplyphysics import (
    units,
    angle_type,
    Quantity,
    SymbolNew,
    validate_input,
    validate_output,
    symbols,
)
from symplyphysics.core.symbols.quantities import scale_factor

# Description
## When an object is under a shearing stress, the shear stress applied to the object is proportional
## to the shear modulus of the object and the angle formed by the face of the object perpendicular to
## the shearing force before (AB) and after the deformation (AB'), see Note for reference.

# Law: tau = G * gamma
## tau - shear stress
## G - shear modulus of the object
## gamma - shear strain (angle of deformation, see Note for reference)

# Conditions
## - The deformation is elastic

# Note
## For a visual representation of shear stress visit [this](https://vuzdoc.org/htm/img/3/6176/223.png)

# Links: Wikipedia <https://en.wikipedia.org/wiki/Shear_stress#Pure>
# TODO: update documentation

shear_stress = symbols.stress
shear_modulus = SymbolNew("G", units.pressure)
# TODO: add to `symbols`

shear_strain = SymbolNew("gamma", angle_type, display_latex="\\gamma")
# TODO: add to `symbols`

law = Eq(shear_stress, shear_modulus * shear_strain)


@validate_input(
    shear_modulus_=shear_modulus,
    shear_strain_=shear_strain,
)
@validate_output(shear_stress)
def calculate_shear_stress(
    shear_modulus_: Quantity,
    shear_strain_: Quantity | float,
) -> Quantity:
    result = law.rhs.subs({
        shear_modulus: shear_modulus_,
        shear_strain: scale_factor(shear_strain_),
    })
    return Quantity(result)
