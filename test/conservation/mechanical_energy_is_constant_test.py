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
from symplyphysics.laws.conservation import mechanical_energy_is_constant as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    Es = Quantity(5 * units.joule)
    Args = namedtuple("Args", ["Es"])
    return Args(Es=Es)


def test_basic_conservation(test_args):
    result_expr = conservation_law.calculate_energy_after(test_args.Es)
    assert SI.get_dimension_system().equivalent_dims(result_expr.dimension, units.energy)
    result = convert_to(result_expr, units.joule).evalf(2)
    assert_approx(result, 5)


def test_bad_energy():
    Eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_energy_after(Eb)
    with raises(TypeError):
        conservation_law.calculate_energy_after(100)
