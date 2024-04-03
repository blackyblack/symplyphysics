from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import heat_capacity_ratio

# Description
## For a gas with isobaric heat capacity of 5 J/K and isochoric heat capacity of 3 J/K
## the heat capacity ratio is approximately 1.67.

Args = namedtuple("Args", "cp cv")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    cp = Quantity(5 * units.joule / units.kelvin)
    cv = Quantity(3 * units.joule / units.kelvin)
    return Args(cp=cp, cv=cv)


def test_law(test_args: Args) -> None:
    result = heat_capacity_ratio.calculate_heat_capacity_ratio(test_args.cp, test_args.cv)
    assert_equal(result, 1.67, tolerance=2e-3)


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        heat_capacity_ratio.calculate_heat_capacity_ratio(cb, test_args.cv)
    with raises(TypeError):
        heat_capacity_ratio.calculate_heat_capacity_ratio(100, test_args.cv)
    with raises(errors.UnitsError):
        heat_capacity_ratio.calculate_heat_capacity_ratio(test_args.cp, cb)
    with raises(TypeError):
        heat_capacity_ratio.calculate_heat_capacity_ratio(test_args.cp, 100)
