from collections import namedtuple
from pytest import fixture, raises
from sympy import pi
from symplyphysics import (errors, units, Quantity, assert_equal, prefixes, quantities)

from symplyphysics.laws.electricity.circuits.waveguides import surface_resistance_of_metal as resistance_law

## Parameters of the coaxial waveguide: the relative permeability of the dielectric is 1,
## the specific conductivity of the conductor is 59.5e6 siemens per meter.
## The angular frequency of signal is 2 * pi * 100e6 radians per second.
## The surface_resistance will be 2.576 milliohm.
## https://old.study.urfu.ru/view/aid/67/1/resonators.pdf

Args = namedtuple("Args", ["absolute_permeability", "angular_frequency", "specific_conductivity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permeability_ = 1
    absolute_permeability = Quantity(relative_permeability_ * quantities.vacuum_permeability)
    angular_frequency = Quantity(2 * pi * 100e6 * (units.radian / units.second))
    specific_conductivity = Quantity(59.5e6 * (units.siemens / units.meter))
    return Args(absolute_permeability=absolute_permeability,
        angular_frequency=angular_frequency,
        specific_conductivity=specific_conductivity)


def test_basic_surface_resistance(test_args: Args) -> None:
    result = resistance_law.calculate_surface_resistance(test_args.absolute_permeability,
        test_args.angular_frequency, test_args.specific_conductivity)
    assert_equal(result, 2.576 * prefixes.milli * units.ohm)


def test_bad_absolute_permeability(test_args: Args) -> None:
    bad_absolute_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_surface_resistance(bad_absolute_permeability,
            test_args.angular_frequency, test_args.specific_conductivity)
    with raises(TypeError):
        resistance_law.calculate_surface_resistance(100,
            test_args.angular_frequency, test_args.specific_conductivity)


def test_bad_angular_frequency(test_args: Args) -> None:
    bad_angular_frequency = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_surface_resistance(test_args.absolute_permeability,
            bad_angular_frequency, test_args.specific_conductivity)
    with raises(TypeError):
        resistance_law.calculate_surface_resistance(test_args.absolute_permeability, 100,
            test_args.specific_conductivity)


def test_bad_specific_conductivity(test_args: Args) -> None:
    specific_conductivity = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        resistance_law.calculate_surface_resistance(test_args.absolute_permeability,
            test_args.angular_frequency, specific_conductivity)
    with raises(TypeError):
        resistance_law.calculate_surface_resistance(test_args.absolute_permeability,
            test_args.angular_frequency, 100)
