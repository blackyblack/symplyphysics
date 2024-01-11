from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
import symplyphysics.laws.hydro.bernoullis_equation as bernoullis_equation


@fixture(name="test_args")
def test_args_fixture():
    inner_pressure_before = Quantity(1 * units.pascal)
    Args = namedtuple("Args", "inner_pressure_before")
    return Args(inner_pressure_before=inner_pressure_before)


def test_bernoullis_equation(test_args):
    result = bernoullis_equation.calculate_inner_pressure(test_args.inner_pressure_before)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).evalf(3)
    assert result_pressure == approx(1, 1e-3)


def test_bad_inner_pressure_before(test_args):
    bad_inner_pressure_before = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        bernoullis_equation.calculate_inner_pressure(bad_inner_pressure_before)
    with raises(TypeError):
        bernoullis_equation.calculate_inner_pressure(100)
