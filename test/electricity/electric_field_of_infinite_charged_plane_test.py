from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import electric_field_of_infinite_charged_plane as intensity_law

# Description
## It is known that with a surface charge density of 2.5 [coulomb / meter^2], the electric field strength is 1.4e11 [volt / meter].
## https://www.calculatoratoz.com/en/electric-field-due-to-infinite-sheet-calculator/Calc-672


@fixture(name="test_args")
def test_args_fixture():
    surface_charge_density = Quantity(2.5 * (units.coulomb / units.meter**2))

    Args = namedtuple("Args", ["surface_charge_density"])
    return Args(surface_charge_density=surface_charge_density)


def test_basic_electric_intensity(test_args):
    result = intensity_law.calculate_electric_intensity(test_args.surface_charge_density)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.voltage / units.length)
    result = convert_to(result, units.volt / units.meter).evalf(5)
    assert result == approx(1.4e11, rel=0.01)


def test_bad_surface_charge_density():
    surface_charge_density = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        intensity_law.calculate_electric_intensity(surface_charge_density)
    with raises(TypeError):
        intensity_law.calculate_electric_intensity(100)
