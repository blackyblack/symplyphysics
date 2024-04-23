from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import temperature_is_derivative_of_internal_energy as temperature_law

# Description
## The internal energy of the system increases by 10 J when the entropy of the system increases by 0.01 J/K.
## The temperature of the system is T = 1000 K.

Args = namedtuple("Args", "du ds")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    du = Quantity(10 * units.joule)
    ds = Quantity(0.01 * units.joule / units.kelvin)
    return Args(du=du, ds=ds)


def test_law(test_args: Args) -> None:
    result = temperature_law.calculate_temperature(test_args.du, test_args.ds)
    assert_equal(result, 1000 * units.kelvin)


def test_bad_energy(test_args: Args) -> None:
    ub = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        temperature_law.calculate_temperature(ub, test_args.ds)
    with raises(TypeError):
        temperature_law.calculate_temperature(100, test_args.ds)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        temperature_law.calculate_temperature(test_args.du, sb)
    with raises(TypeError):
        temperature_law.calculate_temperature(test_args.du, 100)
