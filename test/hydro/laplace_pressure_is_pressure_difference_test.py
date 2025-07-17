from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, errors, units, Quantity
from symplyphysics.laws.hydro import laplace_pressure_is_pressure_difference as law

Args = namedtuple("Args", "p_o p_i")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p_o = Quantity(100 * units.pascal)
    p_i = Quantity(99 * units.pascal)

    return Args(p_o=p_o, p_i=p_i)


def test_law(test_args: Args) -> None:
    result = law.calculate_laplace_pressure(test_args.p_o, test_args.p_i)
    assert_equal(result, 1 * units.pascal)


def test_bad_pressure(test_args: Args) -> None:
    bad_scalar = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_laplace_pressure(bad_scalar, test_args.p_i)
    with raises(TypeError):
        law.calculate_laplace_pressure(100, test_args.p_i)
    with raises(errors.UnitsError):
        law.calculate_laplace_pressure(test_args.p_o, bad_scalar)
    with raises(TypeError):
        law.calculate_laplace_pressure(test_args.p_o, 100)
