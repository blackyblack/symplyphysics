from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.hydro import weight_in_fluid_via_ratio_of_densities as weight

Args = namedtuple("Args", ["weight_vacuum", "liquid_density", "body_density"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    weight_vacuum = Quantity(2.6 * units.newton)
    liquid_density = Quantity(1000 * units.kilogram / units.meter**3)
    body_density = Quantity(2600 * units.kilogram / units.meter**3)
    return Args(weight_vacuum=weight_vacuum, liquid_density=liquid_density, body_density=body_density)


def test_basic_weight(test_args: Args) -> None:
    result = weight.calculate_weight(test_args.weight_vacuum, test_args.liquid_density,
        test_args.body_density)
    assert_equal(result, 1.6 * units.newton)


def test_bad_density(test_args: Args) -> None:
    at = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        weight.calculate_weight(test_args.weight_vacuum, at, test_args.body_density)
    with raises(errors.UnitsError):
        weight.calculate_weight(test_args.weight_vacuum, test_args.liquid_density, at)
    with raises(TypeError):
        weight.calculate_weight(test_args.weight_vacuum, 100, test_args.body_density)
    with raises(TypeError):
        weight.calculate_weight(test_args.weight_vacuum, test_args.liquid_density, 100)


def test_bad_weight(test_args: Args) -> None:
    bad_weight = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        weight.calculate_weight(bad_weight, test_args.liquid_density, test_args.body_density)
    with raises(TypeError):
        weight.calculate_weight(100, test_args.liquid_density, test_args.body_density)
