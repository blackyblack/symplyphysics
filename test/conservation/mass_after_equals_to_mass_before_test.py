from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.conservation import mass_after_equals_to_mass_before as conservation_law

Args = namedtuple("Args", ["ms"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    ms = Quantity(5 * units.kilograms)
    return Args(ms=ms)


def test_basic_conservation(test_args: Args) -> None:
    result = conservation_law.calculate_mass_after(test_args.ms)
    assert_equal(result, 5 * units.kilograms)


def test_bad_mass() -> None:
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_mass_after(mb)
    with raises(TypeError):
        conservation_law.calculate_mass_after(100)
