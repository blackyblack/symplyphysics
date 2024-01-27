from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.core.symbols.celsius import Celsius, to_kelvin_quantity
from symplyphysics.laws.thermodynamics import heat_transfer_equation as heat_flow

# What amount of heat is transferred through the wall from one side to the other if the material has a thermal conductivity of 0.4 Watt/(Kelvin*meter), the cross section thickness is 1 cm, its area is 100 cm**2, the temperature inside is 60 Celsius, the temperature outside is 20 Celsius

@fixture(name="test_args")
def test_args_fixture():
    k = Quantity(0.4 * units.watt / (units.kelvin * units.meter))
    inside_temperature = Celsius(60)
    T1 = to_kelvin_quantity(inside_temperature)
    outside_temperature = Celsius(20)
    T2 = to_kelvin_quantity(outside_temperature)
    d = Quantity(1 * units.centimeter)
    S = Quantity(100 * units.centimeter**2)
    Args = namedtuple("Args", ["k", "T1", "T2", "d", "S"])
    return Args(k=k, T1=T1, T2=T2, d=d, S=S)

def test_basic_amount(test_args):
    result = heat_flow.calculate_heat_flow(test_args.k, test_args.T1, test_args.T2,
        test_args.d, test_args.S)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.power)
    result_energy = convert_to(result, units.watt).evalf(7)
    assert result_energy == approx(16, 0.000001)

def test_bad_material_thermal_conductivity(test_args):
    kb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_flow.calculate_heat_flow(kb, test_args.T1, test_args.T2, test_args.d, test_args.S)
    with raises(TypeError):
        heat_flow.calculate_heat_flow(100, test_args.T1, test_args.T2, test_args.d, test_args.S)

def test_bad_temperature(test_args):
    Tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_flow.calculate_heat_flow(test_args.k, Tb, test_args.T2, test_args.d, test_args.S)
    with raises(TypeError):
        heat_flow.calculate_heat_flow(test_args.k, 100, test_args.T2, test_args.d, test_args.S)
    with raises(errors.UnitsError):
        heat_flow.calculate_heat_flow(test_args.k, test_args.T1, Tb, test_args.d, test_args.S)
    with raises(TypeError):
        heat_flow.calculate_heat_flow(test_args.k, test_args.T1, 100, test_args.d, test_args.S)

def test_bad_cross_section_thickness(test_args):
    db = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_flow.calculate_heat_flow(test_args.k, test_args.T1, test_args.T2, db, test_args.S)
    with raises(TypeError):
        heat_flow.calculate_heat_flow(test_args.k, test_args.T1, test_args.T2, 100, test_args.S)

def test_bad_cross_sectional_area(test_args):
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heat_flow.calculate_heat_flow(test_args.k, test_args.T1, test_args.T2, test_args.d, Sb)
    with raises(TypeError):
        heat_flow.calculate_heat_flow(test_args.k, test_args.T1, test_args.T2, test_args.d, 100)