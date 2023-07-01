from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)
from symplyphysics.laws.thermodynamics import volume_is_constant as isochoric_law


@fixture(name="test_args")
def test_args_fixture():
    t0 = Quantity(1 * units.kelvin)
    t1 = Quantity(2 * units.kelvin)
    P0 = Quantity(1 * units.pascal)
    Args = namedtuple("Args", ["t0", "t1", "P0"])
    return Args(t0=t0, t1=t1, P0=P0)


def test_basic_pressure(test_args):
    result = isochoric_law.calculate_pressure(test_args.t0, test_args.P0, test_args.t1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.pressure)
    result_pressure = convert_to(result, units.pascal).subs(units.pascal, 1).evalf(2)
    assert result_pressure == approx(2.0, 0.01)


def test_bad_temperature(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(tb, test_args.P0, test_args.t1)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(test_args.t0, test_args.P0, tb)
    with raises(TypeError):
        isochoric_law.calculate_pressure(100, test_args.P0, test_args.t1)
    with raises(TypeError):
        isochoric_law.calculate_pressure(test_args.t0, test_args.P0, 100)


def test_bad_pressure(test_args):
    Pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        isochoric_law.calculate_pressure(test_args.t0, Pb, test_args.t1)
    with raises(TypeError):
        isochoric_law.calculate_pressure(test_args.t0, 100, test_args.t1)
