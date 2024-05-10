from typing import Optional
from sympy import Eq
from symplyphysics import (
    units,
    Symbol,
    scale_vector,
)
from symplyphysics.core.fields.scalar_field import ScalarField
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.core.fields.operators import divergence_operator, gradient_operator

# Description
## Heat equation governs heat diffusion, as well as other diffusive processes. It describes
## the evolution of heat transferred from hotter to colder environments in time and space.
## TODO: write about 3d operators

# Law: rho * c_v * dT/dt = div(kappa * grad(T)) + q
## rho - density of medium
## c_v - mass-specific isochoric heat capacity of medium
## kappa - [thermal conductivity](https://en.wikipedia.org/wiki/Thermal_conductivity_and_resistivity#Definition) of medium
## q - density of heat sources
## T - temperature
## t - time
## x - position (spatial variable)
## d/dt - time derivative
## div - divergence operator
## grad - gradient operator

medium_density = Symbol("medium_density", units.mass / units.volume)
medium_specific_heat_capacity = Symbol(
    "medium_specific_heat_capacity",
    units.energy / (units.temperature * units.mass),
)
time = Symbol("time", units.time)


def heat_equation(
    temperature_field: ScalarField,
    thermal_conductivity_field: ScalarField,
    heat_source_density_field: Optional[ScalarField] = None,
) -> Eq:
    temp_diff_time = temperature_field.diff(time).apply_to_basis()

    lhs = medium_density * medium_specific_heat_capacity * temp_diff_time

    grad_temp = gradient_operator(temperature_field)
    kappa = thermal_conductivity_field.apply_to_basis()

    # FIXME: kappa disappears completely, possible mistake in `scale_vector`
    kappa_times_grad_temp = scale_vector(kappa, grad_temp)

    kappa_times_grad_temp_field = VectorField(
        kappa_times_grad_temp.components,
        temperature_field.coordinate_system,
    )

    div_product = divergence_operator(kappa_times_grad_temp_field)

    if heat_source_density_field:
        q = heat_source_density_field.apply_to_basis()
    else:
        q = 0

    rhs = div_product + q

    return Eq(lhs, rhs)
