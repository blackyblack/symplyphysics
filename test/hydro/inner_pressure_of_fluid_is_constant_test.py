from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
import symplyphysics.laws.hydro.inner_pressure_of_fluid_is_constant as bernoullis_equation


@fixture(name="test_args")
def test_args_fixture():
    inner_pressure_before = Quantity(1 * units.pascal)
    Args = namedtuple("Args", "inner_pressure_before")
    return Args(inner_pressure_before=inner_pressure_before)


def test_bernoullis_equation(test_args):
    result = bernoullis_equation.calculate_inner_pressure(test_args.inner_pressure_before)
    assert_equal(result, 1 * units.pascal)


def test_bad_inner_pressure_before():
    bad_inner_pressure_before = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_inner_pressure(bad_inner_pressure_before)
    with raises(TypeError):
        bernoullis_equation.calculate_inner_pressure(100)
