from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    Quantity,
    errors,
    units,
)
from symplyphysics.laws.gravity import keplers_constant_via_attracting_body_mass as law
from symplyphysics.quantities import solar_mass

Args = namedtuple("Args", "m")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    return Args(m=solar_mass)


def test_law(test_args: Args) -> None:
    result = law.calculate_keplers_constant(test_args.m)
    assert_equal(result, 7.5e-6 * units.astronomical_unit**3 / units.day**2)


def test_bad_mass() -> None:
    mb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_keplers_constant(mb)
    with raises(TypeError):
        law.calculate_keplers_constant(100)
