from collections import namedtuple
from pytest import approx, fixture, raises

from sympy import pi
from sympy.physics.units import planck as planck_constant
from symplyphysics import (
    errors,
    units,
    convert_to,
    dimensionless,
    Quantity,
    SI,
)
from symplyphysics.laws.сondensed_matter import effective_mass_of_the_electron_from_the_energy as ef_mass_el

# Description
## Consider the case of an arbitrary dependence of the electron energy on the component of the wave vector.


@fixture(name="test_args")
def test_args_fixture():
    propagation_vec_axis = Quantity(2e8 * pi * (1/units.meter))
    period_length = Quantity(2e-8 * (units.meter))
    mass_1 = Quantity(2.1218e-11 * (units.kilogram))
    mass_2 = Quantity(6.18e-10 * (units.kilogram))
    energy_function = ((planck_constant/2/pi)**2) * (((propagation_vec_axis**2) / 2 / mass_1) + (period_length/mass_2) * propagation_vec_axis**3)
    Args = namedtuple("Args", ["energy_function", "propagation_vec_axis"])
    return Args(energy_function=energy_function,
                propagation_vec_axis=propagation_vec_axis)


def test_basic_mass(test_args):
    result = ef_mass_el.calculate_mass(test_args.energy_function, test_args.propagation_vec_axis)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass)
    result = convert_to(result, units.mass).evalf(6)
    assert result == approx(5.91249e-9, 0.00001)


def test_bad_propagation_vec(test_args):
    propagation_vec_axis = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ef_mass_el.calculate_mass(test_args.energy_function, propagation_vec_axis)
    with raises(TypeError):
        ef_mass_el.calculate_mass(test_args.energy_function, 100)
