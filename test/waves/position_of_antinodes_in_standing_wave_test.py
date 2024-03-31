from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import position_of_antinodes_in_standing_wave as antinode_law

# Description
## For a standing wave of wavelength 8 nm the 4th antinode is located at x = 18 nm.

Args = namedtuple("Args", "m l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = 4
    l = Quantity(8.0 * units.nanometer)
    return Args(m=m, l=l)


def test_law(test_args: Args) -> None:
    result = antinode_law.calculate_antinode_position(test_args.m, test_args.l)
    assert_equal(result, 18.0 * units.nanometer)


def test_bad_integer_factor(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        antinode_law.calculate_antinode_position(mb, test_args.l)


def test_bad_wavelength(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        antinode_law.calculate_antinode_position(test_args.m, lb)
    with raises(TypeError):
        antinode_law.calculate_antinode_position(test_args.m, 100)
