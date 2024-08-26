from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.hydro import submerged_volume_of_floating_body_via_density_ratio as submerged

Args = namedtuple("Args", ["body_volume", "body_density", "liquid_density"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    body_volume = Quantity(12 * units.meter**3)
    body_density = Quantity(500 * units.kilogram / units.meter**3)
    liquid_density = Quantity(1000 * units.kilogram / units.meter**3)
    return Args(body_volume=body_volume, body_density=body_density, liquid_density=liquid_density)


def test_basic_volume(test_args: Args) -> None:
    result = submerged.calculate_submerged_volume(test_args.body_volume, test_args.body_density,
        test_args.liquid_density)
    assert_equal(result, 6 * units.meter**3)


def test_bad_density(test_args: Args) -> None:
    bad_density = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        submerged.calculate_submerged_volume(test_args.body_volume, bad_density,
            test_args.liquid_density)
    with raises(errors.UnitsError):
        submerged.calculate_submerged_volume(test_args.body_volume, test_args.body_density,
            bad_density)
    with raises(TypeError):
        submerged.calculate_submerged_volume(test_args.body_volume, 100, test_args.liquid_density)
    with raises(TypeError):
        submerged.calculate_submerged_volume(test_args.body_volume, test_args.body_density, 100)


def test_bad_volume(test_args: Args) -> None:
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        submerged.calculate_submerged_volume(bad_volume, test_args.body_density,
            test_args.liquid_density)
    with raises(TypeError):
        submerged.calculate_submerged_volume(100, test_args.body_density, test_args.liquid_density)
