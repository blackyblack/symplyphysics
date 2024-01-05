from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.kinematic import impuls as impuls_law

# Description
## Test example from http://kornev-school.ru/f9_momentum.html
## the momentum of a body weighing 6 kg moving at a speed of 4 m/s should be equal to 24 kg * m/s


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(6 * units.kilogram)
    v = Quantity(4 * units.meter / units.second)
    Args = namedtuple("Args", ["m", "v"])
    return Args(m=m, v=v)


def test_basic_impuls(test_args):
    result = impuls_law.calculate_impuls(test_args.m, test_args.v)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass * units.velocity)
    result_acceleration = convert_to(result, units.kilogram * units.meter / units.second).evalf(4)
    assert result_acceleration == approx(24, 0.01)


def test_bad_velocity(test_args):
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impuls_law.calculate_impuls(test_args.m, vb)
    with raises(TypeError):
        impuls_law.calculate_impuls(test_args.m, 100)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        impuls_law.calculate_impuls(mb, test_args.v)
    with raises(TypeError):
        impuls_law.calculate_impuls(100, test_args.v)
