from collections import namedtuple
from pytest import fixture, raises
from sympy import I
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.electricity import admittance_is_conductance_and_susceptance as law

Args = namedtuple("Args", "g b")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    g = Quantity(10 * units.siemens)
    b = Quantity(3 * units.siemens)
    return Args(g=g, b=b)


def test_law(test_args: Args) -> None:
    result = law.calculate_admittance(test_args.g, test_args.b)
    assert_equal(result, (10 + 3 * I) * units.siemens)


def test_bad_conductance(test_args: Args) -> None:
    gb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_admittance(gb, test_args.b)
    with raises(TypeError):
        law.calculate_admittance(100, test_args.b)
    with raises(errors.UnitsError):
        law.calculate_admittance(test_args.g, gb)
    with raises(TypeError):
        law.calculate_admittance(test_args.g, 100)
