from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics import average_square_of_velocity as velocity_law


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(4 * units.atomic_mass_unit)
    T = to_kelvin_quantity(Celsius(25))
    Args = namedtuple("Args", ["m", "T"])
    return Args(m=m, T=T)


def test_basic_energy(test_args):
    result = velocity_law.calculate_rsm_velocity(test_args.m, test_args.T, test_args.M)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity ** 2)
    result_pressure = convert_to(result, units.velocity ** 2).evalf(5)
    assert result_pressure == approx(3 * 1.38 * 1E24 * (25+273) / 4 * 1.66 * 1E-27, 0.01)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_rsm_velocity(test_args.T, mb)
    with raises(TypeError):
        velocity_law.calculate_rsm_velocity(test_args.T, 100)


def test_bad_temperature(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_rsm_velocity(tb, test_args.m)
    with raises(TypeError):
        velocity_law.calculate_rsm_velocity(100, test_args.m)