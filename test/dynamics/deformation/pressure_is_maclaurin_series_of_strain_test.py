from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.dynamics.deformation import (
    pressure_is_maclaurin_series_of_strain as law,)

Args = namedtuple("Args", "e a b s")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(10 * units.pascal)
    a = Quantity(-3 * units.pascal)
    b = Quantity(1 * units.pascal)
    s = 0.05
    return Args(e=e, a=a, b=b, s=s)


def test_law(test_args: Args) -> None:
    result = law.calculate_pressure(test_args.e, test_args.a, test_args.b, test_args.s)
    assert_equal(result, 0.493 * units.pascal)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_pressure(pb, test_args.a, test_args.b, test_args.s)
    with raises(TypeError):
        law.calculate_pressure(100, test_args.a, test_args.b, test_args.s)
    with raises(errors.UnitsError):
        law.calculate_pressure(test_args.e, pb, test_args.b, test_args.s)
    with raises(TypeError):
        law.calculate_pressure(test_args.e, 100, test_args.b, test_args.s)
    with raises(errors.UnitsError):
        law.calculate_pressure(test_args.e, test_args.a, pb, test_args.s)
    with raises(TypeError):
        law.calculate_pressure(test_args.e, test_args.a, 100, test_args.s)


def test_bad_strain(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_pressure(test_args.e, test_args.a, test_args.b, eb)
