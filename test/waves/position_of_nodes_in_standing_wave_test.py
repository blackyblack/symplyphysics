from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.waves import position_of_nodes_in_standing_wave as node_law

# Description
## For a standing wave of wavelength 5 nm the 1000th node is located at x = 2.5 µm.

Args = namedtuple("Args", "m l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m = 1000
    l = Quantity(5.0 * units.nanometer)
    return Args(m=m, l=l)


def test_law(test_args: Args) -> None:
    result = node_law.calculate_node_position(test_args.m, test_args.l)
    assert_equal(result, 2.5 * units.micrometer)


def test_bad_integer_factor(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        node_law.calculate_node_position(mb, test_args.l)


def test_bad_wavelength(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        node_law.calculate_node_position(test_args.m, lb)
    with raises(TypeError):
        node_law.calculate_node_position(test_args.m, 100)
