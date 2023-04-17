from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.thermodynamics import zero_heat_transfer

@fixture
def test_args():
    n = Quantity(units.amount_of_substance, 1 * units.mole)
    t0 = Quantity(units.temperature, 100 * units.kelvin)
    V0 = Quantity(units.volume, 1000 * units.liter)
    V1 = Quantity(units.volume, 2000 * units.liter)
    # Choose specific heats ratio
    y = 1.665
    Args = namedtuple("Args", ["n", "t0", "V0", "V1", "y"])
    return Args(n=n, t0=t0, V0=V0, V1=V1, y=y)

def test_basic_pressure(test_args):
    result = zero_heat_transfer.calculate_pressure(
        test_args.n, test_args.t0, test_args.V0, test_args.V1, test_args.y)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).subs(units.pascal, 1).evalf(8)
    assert result_pressure == approx(262.19, 0.001)

def test_bad_mole_count(test_args):
    nb = Quantity(units.length)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            nb, test_args.t0, test_args.V0, test_args.V1, test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            100, test_args.t0, test_args.V0, test_args.V1, test_args.y)

def test_bad_temperature(test_args):
    tb = Quantity(units.length)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, tb, test_args.V0, test_args.V1, test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, 100, test_args.V0, test_args.V1, test_args.y)

def test_bad_volume(test_args):
    Vb = Quantity(units.length)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, Vb, test_args.V1, test_args.y)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, Vb, test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, 100, test_args.V1, test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, 100, test_args.y)

def test_bad_specific_heats_ratio(test_args):
    yb = Quantity(units.length)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, test_args.V1, yb)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, test_args.V1, 'bad')
