from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes)

from symplyphysics.laws.electricity.circuits import specific_conductivity_of_coaxial_waveguide as conductivity_law

## Parameters of the coaxial waveguide: the angent of the dielectric loss angle is 1e-4, the specific capacitance is 64.5 picofarad per meter.
## The frequency of signal is 2 * pi * 100e6 radians per second. The specific conductivity will be 4.053 microsiemens per meter.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf
## https://test-energy.ru/tangens-ugla-dielektricheskih-poter/

Args = namedtuple("Args", ["frequency", "specific_capacity", "tangent_dielectric_loss_angle"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    frequency = Quantity(2 * pi * 100e6 * (1 / units.second))
    specific_capacity = Quantity(64.5 * prefixes.pico * units.farad / units.meter)
    tangent_dielectric_loss_angle = 1e-4
    return Args(frequency=frequency, specific_capacity=specific_capacity, tangent_dielectric_loss_angle=tangent_dielectric_loss_angle)


def test_basic_specific_conductivity(test_args: Args) -> None:
    result = conductivity_law.calculate_specific_conductivity(test_args.frequency, test_args.specific_capacity,
        test_args.tangent_dielectric_loss_angle)
    assert_equal(result, 4.053 * prefixes.micro * units.siemens / units.meter)


def test_bad_frequency(test_args: Args) -> None:
    bad_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        conductivity_law.calculate_specific_conductivity(bad_frequency, test_args.specific_capacity, test_args.tangent_dielectric_loss_angle)


def test_bad_tangent_dielectric_loss_angle(test_args: Args) -> None:
    tangent_dielectric_loss_angle = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        conductivity_law.calculate_specific_conductivity(test_args.frequency, test_args.specific_capacity, tangent_dielectric_loss_angle)