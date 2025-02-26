from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.gravity.radial_motion import semimajor_axis_via_kepler_constant_and_total_energy as law

Args = namedtuple("Args", "k e")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(7.5e-6 * units.astronomical_unit**3 / units.day**2)
    e = Quantity(-4.4e8 * units.joule / units.kilogram)
    return Args(k=k, e=e)


def test_law(test_args: Args) -> None:
    result = law.calculate_semimajor_axis(test_args.k, test_args.e)
    assert_equal(result, 1.0 * units.astronomical_unit, relative_tolerance=9e-3)


def test_bad_kepler_constant(test_args: Args) -> None:
    kb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_semimajor_axis(kb, test_args.e)
    with raises(TypeError):
        law.calculate_semimajor_axis(100, test_args.e)


def test_bad_energy_per_mass(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_semimajor_axis(test_args.k, eb)
    with raises(TypeError):
        law.calculate_semimajor_axis(test_args.k, 100)
