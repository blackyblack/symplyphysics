from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import pressure_is_constant as gay_lussacs_law


@fixture(name="test_args")
def test_args_fixture():
    t0 = Quantity(1 * units.kelvin)
    t1 = Quantity(2 * units.kelvin)
    V0 = Quantity(1 * units.liter)
    Args = namedtuple("Args", ["t0", "t1", "V0"])
    return Args(t0=t0, t1=t1, V0=V0)


def test_basic_volume(test_args):
    result = gay_lussacs_law.calculate_volume(test_args.t0, test_args.V0, test_args.t1)
    assert_equal(result, 2 * units.liter)


def test_bad_temperature(test_args):
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(tb, test_args.V0, test_args.t1)
    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(test_args.t0, test_args.V0, tb)
    with raises(TypeError):
        gay_lussacs_law.calculate_volume(100, test_args.V0, test_args.t1)
    with raises(TypeError):
        gay_lussacs_law.calculate_volume(test_args.t0, test_args.V0, 100)


def test_bad_volume(test_args):
    Vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gay_lussacs_law.calculate_volume(test_args.t0, Vb, test_args.t1)
    with raises(TypeError):
        gay_lussacs_law.calculate_volume(test_args.t0, 100, test_args.t1)
