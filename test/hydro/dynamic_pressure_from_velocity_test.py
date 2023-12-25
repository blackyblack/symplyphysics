from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.hydro import dynamic_pressure_from_velocity as dynamic_pressure_formula

# Description
## Water with density 1000kg/m3 flows with 0.2m/s. This flow should cause dynamic pressure of 20 Pa. 


@fixture(name="test_args")
def test_args_fixture():
    ro =Quantity(1000 * units.kilogram / units.meter**3) 
    v = Quantity(0.2 * units.meter / units.second)
    Args = namedtuple("Args", ["ro", "v"])
    return Args(ro=ro, v=v)


def test_basic_pressure(test_args):
    result = dynamic_pressure_formula.calculate_pressure(test_args.ro, test_args.v)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_velocity = convert_to(result, units.pascal).evalf(5)
    assert result_velocity == approx(20, 0.0001)

def test_bad_density(test_args):
    bd = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dynamic_pressure_formula.calculate_pressure(bd, test_args.v)
    with raises(TypeError):
        dynamic_pressure_formula.calculate_pressure(100, test_args.ro)

def test_bad_velocity(test_args):
    bv = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        dynamic_pressure_formula.calculate_pressure(test_args.ro, bv)
    with raises(TypeError):
        dynamic_pressure_formula.calculate_pressure(test_args.ro, 100)