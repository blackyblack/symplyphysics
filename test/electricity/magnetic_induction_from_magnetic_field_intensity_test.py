from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes, quantities)
from symplyphysics.laws.electricity import magnetic_induction_from_magnetic_field_intensity as induction_law

# Description
## With a relative magnetic permeability of 1592 and a magnetic field intensity of 1 ampere per meter,
## the magnetic induction is 2 millitesla.
## https://www.calculatoratoz.com/en/magnetic-permeability-calculator/Calc-2144

Args = namedtuple("Args", ["absolute_permeability", "intensity"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permeability_ = 1592
    absolute_permeability = Quantity(relative_permeability_ * quantities.vacuum_permeability)
    intensity = Quantity(1 * (units.ampere / units.meter))
    return Args(absolute_permeability=absolute_permeability, intensity=intensity)


def test_basic_induction(test_args: Args) -> None:
    result = induction_law.calculate_induction(test_args.absolute_permeability, test_args.intensity)
    assert_equal(result, 2 * prefixes.milli * units.tesla)


def test_bad_absolute_permeability(test_args: Args) -> None:
    absolute_permeability = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(absolute_permeability, test_args.intensity)
    with raises(TypeError):
        induction_law.calculate_induction(100, test_args.intensity)


def test_bad_intensity(test_args: Args) -> None:
    intensity = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.absolute_permeability, intensity)
    with raises(TypeError):
        induction_law.calculate_induction(test_args.absolute_permeability, 100)
