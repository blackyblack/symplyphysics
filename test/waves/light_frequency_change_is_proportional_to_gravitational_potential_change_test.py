from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import (
    light_frequency_change_is_proportional_to_gravitational_potential_change as law,)

Args = namedtuple("Args", "nu dphi")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    nu = Quantity(6e14 * units.hertz)
    dphi = Quantity(1e-9 * (units.meter / units.second)**2)
    return Args(nu=nu, dphi=dphi)


def test_law(test_args: Args) -> None:
    result = law.calculate_frequency_change(test_args.nu, test_args.dphi)
    assert_equal(result, -6.7e-12 * units.hertz, tolerance=4e-3)


def test_bad_frequency(test_args: Args) -> None:
    nub = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_frequency_change(nub, test_args.dphi)
    with raises(TypeError):
        law.calculate_frequency_change(100, test_args.dphi)


def test_bad_potential(test_args: Args) -> None:
    phib = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_frequency_change(test_args.nu, phib)
    with raises(TypeError):
        law.calculate_frequency_change(test_args.nu, 100)
