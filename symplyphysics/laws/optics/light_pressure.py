from sympy import (Eq, solve)
from sympy.physics.units import speed_of_light
from symplyphysics import (units, Quantity, Symbol, validate_input,
    validate_output, dimensionless)

# Description
## Light pressure is the pressure exerted by light radiation incident on the surface of a body.
## Under normal light incidence on the surface, the light pressure depends only on the reflection
## coefficient and the intensity of the radiation.

## Law is: p = I * (1 + R) / c, where
## p - light pressure,
## I - intensity of incident radiation,
## R - surface reflection coefficient,
## c - speed of light in vacuum.

# Links: BYJU's, "Radiation pressure formula" <https://byjus.com/physics/radiation-pressure/>

pressure = Symbol("pressure", units.pressure)

intensity = Symbol("intensity", units.power / units.area)
reflection_coefficient = Symbol("reflection_coefficient", dimensionless)

law = Eq(pressure, intensity * (1 + reflection_coefficient) / speed_of_light)


@validate_input(intensity_=intensity, reflection_coefficient_=reflection_coefficient)
@validate_output(pressure)
def calculate_pressure(intensity_: Quantity, reflection_coefficient_: float) -> Quantity:
    result_expr = solve(law, pressure, dict=True)[0][pressure]
    result_expr = result_expr.subs({
        intensity: intensity_,
        reflection_coefficient: reflection_coefficient_,
    })
    return Quantity(result_expr)
