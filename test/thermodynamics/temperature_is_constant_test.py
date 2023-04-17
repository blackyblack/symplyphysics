from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.thermodynamics import temperature_is_constant as boyles_law

@fixture
def test_args():
    P0 = Quantity(units.pressure, 1 * units.pascal)
    P1 = Quantity(units.pressure, 2 * units.pascal)
    V0 = Quantity(units.volume, 1 * units.liter)
    Args = namedtuple("Args", ["P0", "P1", "V0"])
    return Args(P0=P0, P1=P1, V0=V0)

def test_basic_volume(test_args):
    result = boyles_law.calculate_volume(test_args.P0, test_args.V0, test_args.P1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume)
    result_volume = convert_to(result, units.liter).subs(units.liter, 1).evalf(2)
    assert result_volume == approx(0.5, 0.01)

def test_bad_pressure(test_args):
    Pb = Quantity(units.length)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(Pb, test_args.V0, test_args.P1)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(test_args.P0, test_args.V0, Pb)
    with raises(TypeError):
        boyles_law.calculate_volume(100, test_args.V0, test_args.P1)
    with raises(TypeError):
        boyles_law.calculate_volume(test_args.P0, test_args.V0, 100)

def test_bad_volume(test_args):
    Vb = Quantity(units.length)
    with raises(errors.UnitsError):
        boyles_law.calculate_volume(test_args.P0, Vb, test_args.P1)
    with raises(TypeError):
        boyles_law.calculate_volume(test_args.P0, 100, test_args.P1)
