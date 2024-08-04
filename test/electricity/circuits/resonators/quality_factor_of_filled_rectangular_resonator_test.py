from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.circuits.resonators import quality_factor_of_filled_rectangular_resonator as factor_law

# Description
## The quality factor of an empty resonator is 1000. The tangent of the dielectric loss angle of the material
## filling the resonator is 1e-4. Then the quality factor of the filled resonator will be equal to 909.

Args = namedtuple("Args", ["empty_resonator_quality_factor", "tangent_dielectric_loss_angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    empty_resonator_quality_factor = 1000
    tangent_dielectric_loss_angle = 1e-4

    return Args(empty_resonator_quality_factor=empty_resonator_quality_factor,
        tangent_dielectric_loss_angle=tangent_dielectric_loss_angle)


def test_basic_quality_factor(test_args: Args) -> None:
    result = factor_law.calculate_quality_factor(test_args.empty_resonator_quality_factor,
        test_args.tangent_dielectric_loss_angle)
    assert_equal(result, 909)


def test_bad_empty_resonator_quality_factor(test_args: Args) -> None:
    empty_resonator_quality_factor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(empty_resonator_quality_factor,
            test_args.tangent_dielectric_loss_angle)


def test_bad_tangent_dielectric_loss_angle(test_args: Args) -> None:
    tangent_dielectric_loss_angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.empty_resonator_quality_factor,
            tangent_dielectric_loss_angle)
