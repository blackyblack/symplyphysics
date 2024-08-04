from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (
    errors,
    units,
    Quantity,
    assert_equal,
)

from symplyphysics.laws.electricity.circuits.waveguides import attenuation_coefficient_in_dielectric as coefficient_law

## Parameters of the coaxial waveguide: the relative permittivity of the dielectric is 2.2, the relative permeability of the dielectric is 1,
## the tangent of the dielectric loss angle is 1e-4. The angular frequency of signal is 2 * pi * 100e6 radians per second.
## The attenuation coefficient will be 155.4e-6 [1 / meter].
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", [
    "relative_permittivity", "relative_permeability", "angular_frequency",
    "tangent_dielectric_loss_angle"
])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permittivity = 2.2
    relative_permeability = 1
    angular_frequency = Quantity(2 * pi * 100e6 * (units.radian / units.second))
    tangent_dielectric_loss_angle = 1e-4
    return Args(relative_permittivity=relative_permittivity,
        relative_permeability=relative_permeability,
        angular_frequency=angular_frequency,
        tangent_dielectric_loss_angle=tangent_dielectric_loss_angle)


def test_basic_attenuation_coefficient(test_args: Args) -> None:
    result = coefficient_law.calculate_attenuation_coefficient(
        test_args.relative_permittivity, test_args.relative_permeability,
        test_args.angular_frequency, test_args.tangent_dielectric_loss_angle)
    assert_equal(result, 155.4e-6 * (1 / units.meter))


def test_bad_relative_permittivity(test_args: Args) -> None:
    bad_relative_permittivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(bad_relative_permittivity,
            test_args.relative_permeability, test_args.angular_frequency,
            test_args.tangent_dielectric_loss_angle)


def test_bad_relative_permeability(test_args: Args) -> None:
    bad_relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(bad_relative_permeability,
            test_args.relative_permeability, test_args.angular_frequency,
            test_args.tangent_dielectric_loss_angle)


def test_bad_angular_frequency(test_args: Args) -> None:
    bad_angular_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
            test_args.relative_permeability, bad_angular_frequency,
            test_args.tangent_dielectric_loss_angle)
    with raises(TypeError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
            test_args.relative_permeability, 100, test_args.tangent_dielectric_loss_angle)


def test_bad_tangent_dielectric_loss_angle(test_args: Args) -> None:
    tangent_dielectric_loss_angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_attenuation_coefficient(test_args.relative_permittivity,
            test_args.relative_permeability, test_args.angular_frequency,
            tangent_dielectric_loss_angle)
