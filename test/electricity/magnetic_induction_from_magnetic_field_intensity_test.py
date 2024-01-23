from collections import namedtuple
from pytest import approx, fixture, raises
from sympy.physics.units import prefixes
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import magnetic_induction_from_magnetic_field_intensity as induction_law

# Description
## With a relative magnetic permeability of 1592 and a magnetic field intensity of 1 ampere per meter,
## the magnetic induction is 2 millitesla.
## https://www.calculatoratoz.com/en/magnetic-permeability-calculator/Calc-2144


@fixture(name="test_args")
def test_args_fixture():
    relative_permeability = 1592
    intensity = Quantity(1 * (units.ampere / units.meter))

    Args = namedtuple("Args", ["relative_permeability", "intensity"])
    return Args(relative_permeability=relative_permeability, intensity=intensity)


def test_basic_induction(test_args):
    result = induction_law.calculate_induction(test_args.relative_permeability, test_args.intensity)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.magnetic_density)
    result = convert_to(result, prefixes.milli * units.tesla).evalf(5)
    assert result == approx(2, rel=0.01)


def test_bad_relative_permeability(test_args):
    relative_permeability = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(relative_permeability, test_args.intensity)
    with raises(TypeError):
        induction_law.calculate_induction(True, test_args.intensity)


def test_bad_intensity(test_args):
    intensity = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.relative_permeability, intensity)
    with raises(TypeError):
        induction_law.calculate_induction(test_args.relative_permeability, 100)
