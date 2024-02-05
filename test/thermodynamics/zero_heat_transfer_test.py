from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_approx,
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.thermodynamics import zero_heat_transfer


@fixture(name="test_args")
def test_args_fixture():
    n = Quantity(1 * units.mole)
    t0 = Quantity(100 * units.kelvin)
    V0 = Quantity(1000 * units.liter)
    V1 = Quantity(2000 * units.liter)
    # Choose specific heats ratio
    y = 1.665
    Args = namedtuple("Args", ["n", "t0", "V0", "V1", "y"])
    return Args(n=n, t0=t0, V0=V0, V1=V1, y=y)


def test_basic_pressure(test_args):
    result = zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0,
        test_args.V1, test_args.y)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).evalf(8)
    assert_approx(result_pressure, 262.19)


def test_bad_mole_count(test_args):
    nb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(nb, test_args.t0, test_args.V0, test_args.V1,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(100, test_args.t0, test_args.V0, test_args.V1,
            test_args.y)


def test_bad_temperature(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, tb, test_args.V0, test_args.V1,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(test_args.n, 100, test_args.V0, test_args.V1,
            test_args.y)


def test_bad_volume(test_args):
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, Vb, test_args.V1,
            test_args.y)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, Vb,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, 100, test_args.V1,
            test_args.y)
    with raises(TypeError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, 100,
            test_args.y)


def test_bad_specific_heats_ratio(test_args):
    yb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, test_args.V1,
            yb)
    with raises(errors.UnitsError):
        zero_heat_transfer.calculate_pressure(test_args.n, test_args.t0, test_args.V0, test_args.V1,
            'bad')
