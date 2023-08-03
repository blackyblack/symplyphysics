from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.thermodynamics import irradiation_of_black_body_from_temperature as stefan_boltzmann_law


@fixture(name="test_args")
def test_args_fixture():
    t = Quantity(20 * units.kelvin)
    Args = namedtuple("Args", ["t"])
    return Args(t=t)

def test_basic_radiant_heat_energy(test_args):
    result = stefan_boltzmann_law.calculate_irradiance(test_args.t)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power/units.length**2)
    result_radiant_heat_energy_joule = convert_to(result, units.joule/units.meter**2/units.second).evalf(6)
    assert result_radiant_heat_energy_joule == approx(0.009073, 0.000001)

def test_bad_temperature(test_args):
    tb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        stefan_boltzmann_law.calculate_irradiance(tb)
    with raises(TypeError):
        stefan_boltzmann_law.calculate_irradiance(100)