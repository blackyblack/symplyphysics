from sympy import (
    Eq,
    Derivative,
    symbols as sym_symbols,
    Function as SymFunction,
    solve,
)
from symplyphysics import (
    symbols,
    units,
    Quantity,
    Symbol,
    Function,
    validate_input,
    validate_output,
)
from symplyphysics.core.expr_comparisons import expr_equals
from symplyphysics.laws.thermodynamics import gibbs_energy_differential

# Description
## Maxwell relations are a set of equations that unite the most common thermodynamic quantities between
## one another. They are derived from the fundamental themodynamic relations featuring differentials
## of thermodynamic potentials, and this method of derivation is called the method of thermodynamic
## potentials.

# Law: (dS/dp)_T = -(dV/dT)_p
## S - entropy
## p - pressure
## V - volume
## T - absolute temperature

# Conditions
## - Changes in particle count are not taken into account and it is assumed to be constant.

entropy = Function("entropy", units.energy / units.temperature)
pressure = Symbol("pressure", units.pressure)
volume = Function("volume", units.volume)
temperature = symbols.thermodynamics.temperature

law = Eq(
    Derivative(entropy(temperature, pressure), pressure),
    -1 * Derivative(volume(temperature, pressure), temperature)
)

# Derive from expression of Gibbs energy differential

_gibbs_energy_expr = sym_symbols("gibbs_energy", cls=SymFunction)(temperature, pressure)

_gibbs_energy_diff_temperature = gibbs_energy_differential.law.rhs.subs({
    gibbs_energy_differential.temperature_change: 1,
    gibbs_energy_differential.pressure_change: 0,
    gibbs_energy_differential.particle_count_change: 0,
    gibbs_energy_differential.entropy: entropy(temperature, pressure)
})

_gibbs_energy_diff_temperature_diff_pressure = _gibbs_energy_diff_temperature.diff(pressure)

_eqn_temperature_pressure = Eq(
    _gibbs_energy_expr.diff(temperature, pressure),
    _gibbs_energy_diff_temperature_diff_pressure
)

_gibbs_energy_diff_pressure = gibbs_energy_differential.law.rhs.subs({
    gibbs_energy_differential.temperature_change: 0,
    gibbs_energy_differential.pressure_change: 1,
    gibbs_energy_differential.particle_count_change: 0,
    gibbs_energy_differential.volume: volume(temperature, pressure)
})

_gibbs_energy_diff_pressure_diff_temperature = _gibbs_energy_diff_pressure.diff(temperature)

_eqn_pressure_temperature = Eq(
    _gibbs_energy_expr.diff(pressure, temperature),
    _gibbs_energy_diff_pressure_diff_temperature,
)

_entropy_diff_pressure = solve(
    (_eqn_temperature_pressure, _eqn_pressure_temperature),
    (entropy(temperature, pressure).diff(pressure), _gibbs_energy_expr.diff(temperature, pressure)),
    dict=True,
)[0][entropy(temperature, pressure).diff(pressure)]

assert expr_equals(_entropy_diff_pressure, law.rhs)


@validate_input(
    volume_change_=volume,
    temperature_change_=temperature,
)
@validate_output(units.energy / (units.temperature * units.pressure))
def calculate_entropy_differential(
    volume_change_: Quantity,
    temperature_change_: Quantity,
) -> Quantity:
    volume_ = (volume_change_ / temperature_change_) * temperature
    result = law.rhs.subs(volume(temperature, pressure), volume_).doit()
    return Quantity(result)
