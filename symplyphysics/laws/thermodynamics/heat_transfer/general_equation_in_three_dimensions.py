from typing import Sequence
from sympy import Eq, Derivative, S, Expr
from symplyphysics import (
    units,
    Symbol,
    Function,
)
from symplyphysics.core.dimensions import ScalarValue
from symplyphysics.core.fields.scalar_field import ScalarField, T as PointType
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
    units.energy / (units.temperature * units.mass)
)
thermal_conductivity = Function(
    "thermal_conductivity",
    units.power / (units.length * units.temperature),
)
time = Symbol("time", units.time)

def heat_equation(temperature_field: ScalarField, heat_source_density: ScalarField) -> Eq:
    temperature_function = temperature_field.field_function
    _callable = callable(temperature_function)

    def _temperature_diff_time(point: PointType) -> ScalarValue:
        temperature_value = temperature_function(point) if _callable else temperature_function
        return Derivative(temperature_value, time)

    temperature_diff_time = ScalarField(
        _temperature_diff_time,
        temperature_field.coordinate_system,
    )

    grad_temperature = gradient_operator(temperature_field).components

    base_scalars = temperature_field.coordinate_system.coord_system.base_scalars()

    def _kappa_times_grad(point: PointType) -> Sequence[ScalarValue]:
        kappa_applied = thermal_conductivity(*point.coordinates)
        return tuple(
            kappa_applied * grad_component.subs(base_scalar, coordinate)
            for grad_component, base_scalar, coordinate in zip(
                grad_temperature, base_scalars, point.coordinates
            )
        )
    
    kappa_times_grad = VectorField(_kappa_times_grad, temperature_field.coordinate_system)

    div_product = divergence_operator(kappa_times_grad)

    # FIXME: finish this function
