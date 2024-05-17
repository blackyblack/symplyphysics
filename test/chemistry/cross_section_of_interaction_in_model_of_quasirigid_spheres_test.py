from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.chemistry import cross_section_of_interaction_in_model_of_quasirigid_spheres as cross_section_law

# Description
## The distance of the closest approach of two particles is 1.63e-10 meter.
## Then the cross-sectional area of the interaction is 8.35e-20 meter^2.

Args = namedtuple("Args", ["distance_of_convergence_of_two_particles"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    distance_of_convergence_of_two_particles = Quantity(1.63e-10 * units.meter)

    return Args(distance_of_convergence_of_two_particles=distance_of_convergence_of_two_particles)


def test_basic_cross_sectional_area_of_interaction(test_args: Args) -> None:
    result = cross_section_law.calculate_cross_sectional_area_of_interaction(
        test_args.distance_of_convergence_of_two_particles)
    assert_equal(result, 8.35e-20 * units.meter**2)


def test_bad_distance_of_convergence_of_two_particles() -> None:
    distance_of_convergence_of_two_particles = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(
            distance_of_convergence_of_two_particles)
    with raises(TypeError):
        cross_section_law.calculate_cross_sectional_area_of_interaction(100)
