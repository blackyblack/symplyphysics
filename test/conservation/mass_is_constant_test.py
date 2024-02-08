from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.conservation import mass_is_constant as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    ms = Quantity(5 * units.kilograms)
    Args = namedtuple("Args", ["ms"])
    return Args(ms=ms)


def test_basic_conservation(test_args):
    result = conservation_law.calculate_mass_after(test_args.ms)
    assert_equal(result, 5 * units.kilograms)


def test_bad_mass():
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_mass_after(mb)
    with raises(TypeError):
        conservation_law.calculate_mass_after(100)
