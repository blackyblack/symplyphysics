from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.thermodynamics import average_square_of_velocity as velocity_law

# Description
## The ideal gas has a temperature of 300 K. The mass of gas is equal to 32 gram. Then the average square of velocity should be 3.881 * 1E-21 (m/s)^2


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(32 * units.gram)
    t = Quantity(300 * units.kelvin)
    Args = namedtuple("Args", ["m", "t"])
    return Args(m=m, t=t)


def test_basic_average_square_velocity(test_args):
    result = velocity_law.calculate_average_square_velocity(test_args.t,
                                                            test_args.m)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.velocity ** 2)
    result_average_square_velocity = convert_to(result, (units.meter / units.second) ** 2).evalf(5)
    assert result_average_square_velocity == approx(3.881 * 1E-21, 0.001)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_average_square_velocity(test_args.t, mb)
    with raises(TypeError):
        velocity_law.calculate_average_square_velocity(test_args.t, 100)


def test_bad_temperature(test_args):
    tb = Quantity(100 * units.coulomb)
    with raises(errors.UnitsError):
        velocity_law.calculate_average_square_velocity(tb, test_args.m)
    with raises(TypeError):
        velocity_law.calculate_average_square_velocity(100, test_args.m)
