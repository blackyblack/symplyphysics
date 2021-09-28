from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.thermodynamics import zero_heat_transfer

@fixture
def test_args():
    n = units.Quantity('n')
    SI.set_quantity_dimension(n, units.amount_of_substance)
    SI.set_quantity_scale_factor(n, 1 * units.mole)
    t0 = units.Quantity('t0')
    SI.set_quantity_dimension(t0, units.temperature)
    SI.set_quantity_scale_factor(t0, 100 * units.kelvin)
    V0 = units.Quantity('V0')
    SI.set_quantity_dimension(V0, units.volume)
    SI.set_quantity_scale_factor(V0, 1000 * units.liter)
    V1 = units.Quantity('V1')
    SI.set_quantity_dimension(V1, units.volume)
    SI.set_quantity_scale_factor(V1, 2000 * units.liter)
    # Choose specific heats ratio
    y = 1.665

    Args = namedtuple('Args', ['n', 't0', 'V0', 'V1', 'y'])
    return Args(n=n, t0=t0, V0=V0, V1=V1, y=y)

def test_basic_pressure(test_args):
    result = zero_heat_transfer.calculate_pressure(
        test_args.n, test_args.t0, test_args.V0, test_args.V1, test_args.y)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)

    result_pressure = convert_to(
        result, units.pascal).subs(units.pascal, 1).evalf(8)
    assert result_pressure == approx(262.19, 0.001)

def test_bad_mole_count(test_args):
    nb = units.Quantity('nb')
    SI.set_quantity_dimension(nb, units.length)
    SI.set_quantity_scale_factor(nb, 1 * units.meter)

    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            nb, test_args.t0, test_args.V0, test_args.V1, test_args.y)

    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            100, test_args.t0, test_args.V0, test_args.V1, test_args.y)


def test_bad_temperature(test_args):
    t0b = units.Quantity('t0b')
    SI.set_quantity_dimension(t0b, units.length)
    SI.set_quantity_scale_factor(t0b, 1 * units.meter)

    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, t0b, test_args.V0, test_args.V1, test_args.y)

    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, 100, test_args.V0, test_args.V1, test_args.y)


def test_bad_volume(test_args):
    # Make V0 invalid
    V0b = units.Quantity('V0b')
    SI.set_quantity_dimension(V0b, units.length)
    SI.set_quantity_scale_factor(V0b, 1 * units.meter)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, V0b, test_args.V1, test_args.y)

    # Make V1 invalid
    V1b = units.Quantity('V1b')
    SI.set_quantity_dimension(V1b, units.length)
    SI.set_quantity_scale_factor(V1b, 1 * units.meter)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, V1b, test_args.y)

    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, 100, test_args.V1, test_args.y)

    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, 100, test_args.y)

def test_bad_specific_heats_ratio(test_args):
    yb = units.Quantity('yb')
    SI.set_quantity_dimension(yb, units.length)
    SI.set_quantity_scale_factor(yb, 1 * units.meter)

    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, test_args.V1, yb)

    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(
            test_args.n, test_args.t0, test_args.V0, test_args.V1, 'bad')
