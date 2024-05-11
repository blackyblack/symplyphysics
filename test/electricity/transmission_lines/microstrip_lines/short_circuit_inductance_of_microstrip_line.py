from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
    prefixes
)

from symplyphysics.laws.electricity.transmission_lines.microstrip_lines import short_circuit_inductance_of_microstrip_line as inductance_law

## The thickness of the substrate is 600 millimeter. The radius of the hole in the substrate is 2.5 millimeter.
## Then the inductance of the metallized hole is equal to 562 nanohenry.

Args = namedtuple("Args", ["thickness_of_substrate", "radius_of_hole",])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    thickness_of_substrate = Quantity(600 * units.millimeter)
    radius_of_hole = Quantity(2.5 * units.millimeter)
    return Args(thickness_of_substrate=thickness_of_substrate,
        radius_of_hole=radius_of_hole,)


def test_basic_inductance(test_args: Args) -> None:
    result = inductance_law.calculate_inductance(test_args.thickness_of_substrate, test_args.radius_of_hole)
    assert_equal(result, 562 * prefixes.nano * units.henry)


def test_bad_thickness_of_substrate(test_args: Args) -> None:
    bad_thickness_of_substrate = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(bad_thickness_of_substrate, test_args.radius_of_hole)
    with raises(TypeError):
        inductance_law.calculate_inductance(100, test_args.radius_of_hole)


def test_bad_radius_of_hole(test_args: Args) -> None:
    bad_radius_of_hole = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        inductance_law.calculate_inductance(test_args.thickness_of_substrate, bad_radius_of_hole)
    with raises(TypeError):
        inductance_law.calculate_inductance(test_args.thickness_of_substrate, 100)
