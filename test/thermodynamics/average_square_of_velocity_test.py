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

# Description
## The ideal gas has a temperature of 298 K. The mass of an atom of a given gas is equal to 4 atomic units of mass. Then the average square of velocity should be 1.858 * 10^6 (m/s)^2


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(4 * units.atomic_mass_unit)
    temperature = Celsius(25)
    t = to_kelvin_quantity(temperature)
    Args = namedtuple("Args", ["m", "t"])
    return Args(m=m, t=t)


def test_basic_average_square_velocity(test_args):
    result = velocity_law.calculate_rsm_velocity(test_args.t, test_args.m)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity ** 2)
    result_pressure = convert_to(result, units.velocity ** 2).evalf(5)
    assert result_pressure == approx(1.858 * 1E6, 0.001)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_rsm_velocity(test_args.t, mb)
    with raises(TypeError):
        velocity_law.calculate_rsm_velocity(test_args.t, 100)


def test_bad_temperature(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_rsm_velocity(tb, test_args.m)
    with raises(TypeError):
        velocity_law.calculate_rsm_velocity(100, test_args.m)
