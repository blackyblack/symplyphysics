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
from symplyphysics.laws.сondensed_matter_physics import effective_mass_of_the_electron_from_the_energy as ef_mass_el

# Description
## Consider the case of the quadratic isotropic law of dispersion of charge carriers.


@fixture(name="test_args")
def test_args_fixture():
    propogation_vec_axis = Quantity(2e8 * pi * (dimensionless/units.length))
    energy_func = ((planck_constant/2/pi)**2) * (propogation_vec_axis**2) / 2 / Quantity(1 * units.mass) / 2.1218e-11
    Args = namedtuple("Args", ["energy_func", "propogation_vec_axis"])
    return Args(energy_func=energy_func,
                propogation_vec_axis=propogation_vec_axis)


def test_basic_mass(test_args):
    result = ef_mass_el.calculate_mass(test_args.energy_func, test_args.propogation_vec_axis)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass)
    result_current = convert_to(result, units.mass).evalf(6)
    assert result_current == approx(2.1218e-11, 0.00001)


def test_bad_propogation_vec(test_args):
    propogation_vec_axis = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        ef_mass_el.calculate_mass(test_args.energy_func, propogation_vec_axis)
    with raises(TypeError):
        ef_mass_el.calculate_mass(test_args.energy_func, 100)
