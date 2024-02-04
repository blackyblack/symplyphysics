from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.definitions import momentum_is_mass_times_velocity as momentum_def


@fixture(name="test_args")
def test_args_fixture():
    m = Quantity(1 * units.kilogram)
    v = Quantity(5 * units.meter / units.second)
    Args = namedtuple("Args", ["m", "v"])
    return Args(m=m, v=v)


def test_basic_momentum(test_args):
    result = momentum_def.calculate_momentum(test_args.m, test_args.v)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.momentum)
    result_momentum = convert_to(result, momentum_def.definition_units_SI).evalf(2)
    assert_approx(result_momentum, 5)


def test_bad_mass(test_args):
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(mb, test_args.v)
    with raises(TypeError):
        momentum_def.calculate_momentum(100, test_args.v)


def test_bad_velocity(test_args):
    vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        momentum_def.calculate_momentum(test_args.m, vb)
    with raises(TypeError):
        momentum_def.calculate_momentum(test_args.m, 100)
