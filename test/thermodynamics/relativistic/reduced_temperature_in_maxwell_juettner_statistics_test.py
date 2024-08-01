from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics.relativistic import (
    reduced_temperature_in_maxwell_juettner_statistics as law,
)

Args = namedtuple("Args", "t m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(1e9 * units.kelvin)
    m = units.electron_rest_mass
    return Args(t=t, m=m)


def test_law(test_args: Args) -> None:
    result = law.calculate_reduced_temperature(test_args.t, test_args.m)
    assert_equal(result, 0.17, tolerance=9e-3)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_reduced_temperature(tb, test_args.m)
    with raises(TypeError):
        law.calculate_reduced_temperature(100, test_args.m)


def test_bad_mass(test_args: Args) -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_reduced_temperature(test_args.t, mb)
    with raises(TypeError):
        law.calculate_reduced_temperature(test_args.t, 100)
