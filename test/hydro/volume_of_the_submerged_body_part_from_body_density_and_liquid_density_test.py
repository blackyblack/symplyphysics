from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.hydro import volume_of_the_submerged_body_part_from_body_density_and_liquid_density as submerged

## If the density of the body is greater than the density of the liquid, then the volume of the submerged part of the body will be greater than the part of the body on the surface, since gravity is greater than the Archimedes force.
## then the volume of the submerged part of the body will be greater than the part of the body on the surface, since gravity is greater than the Archimedes force.

@fixture(name="test_args")
def test_args_fixture():
    body_volume = Quantity(0.04 * units.meter**3)
    body_density = Quantity(2600 * units.kilogram / units.meter**3)
    liquid_density = Quantity(1000 * units.kilogram / units.meter**3)
    Args = namedtuple("Args", ["body_volume", "body_density", "liquid_density"])
    return Args(
        body_volume=body_volume,
        body_density=body_density,
        liquid_density=liquid_density
               )


def test_basic_volume(test_args):
    result = submerged.calculate_submerged_volume(test_args.body_volume, test_args.body_density, test_args.liquid_density)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume)
    result_volume = convert_to(result, units.meter**3).evalf(5)
    assert result_volume == approx(0.104, 0.001)


def test_bad_density(test_args):
    bad_density = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        submerged.calculate_submerged_volume(test_args.body_volume, bad_density, test_args.liquid_density)
    with raises(errors.UnitsError):
        submerged.calculate_submerged_volume(test_args.body_volume, test_args.body_density, bad_density)
    with raises(TypeError):
        submerged.calculate_submerged_volume(test_args.body_volume, 100, test_args.liquid_density)
    with raises(TypeError):
        submerged.calculate_submerged_volume(test_args.body_volume, test_args.body_density, 100)


def test_bad_volume(test_args):
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        submerged.calculate_submerged_volume(bad_volume, test_args.body_density, test_args.liquid_density)
    with raises(TypeError):
        submerged.calculate_submerged_volume(100, test_args.body_density, test_args.liquid_density)
