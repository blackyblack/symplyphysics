from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.thermodynamics import temperature_derivative_via_volume_derivative as joule_thompson_law

# A certain amount of gas with isobaric heat capacity C_p = 1000 J/K undergoes a differential
## Joule-Thompson process. Initial parameters are V = 1 m**3 at T = 300 K and V = 1.00001 m**3
## at T = 300.001 K. The temperature derivative w.r.t. pressure during this process is 2 mK/Pa.

Args = namedtuple("Args", "v0 v1 t0 t1 cp")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v0 = Quantity(1. * units.meter**3)
    v1 = Quantity(1.00001 * units.meter**3)
    t0 = Quantity(300. * units.kelvin)
    t1 = Quantity(300.001 * units.kelvin)
    cp = Quantity(1000 * units.joule / units.kelvin)
    return Args(v0=v0, v1=v1, t0=t0, t1=t1, cp=cp)


def test_law(test_args: Args) -> None:
    result = joule_thompson_law.calculate_temperature_derivative(test_args.v0, test_args.v1,
        test_args.t0, test_args.t1, test_args.cp)
    assert_equal(result, 2 * prefixes.milli * units.kelvin / units.pascal)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        joule_thompson_law.calculate_temperature_derivative(vb, test_args.v1, test_args.t0,
            test_args.t1, test_args.cp)
    with raises(TypeError):
        joule_thompson_law.calculate_temperature_derivative(100, test_args.v1, test_args.t0,
            test_args.t1, test_args.cp)
    with raises(errors.UnitsError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, vb, test_args.t0,
            test_args.t1, test_args.cp)
    with raises(TypeError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, 100, test_args.t0,
            test_args.t1, test_args.cp)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, test_args.v1, tb,
            test_args.t1, test_args.cp)
    with raises(TypeError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, test_args.v1, 100,
            test_args.t1, test_args.cp)
    with raises(errors.UnitsError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, test_args.v1,
            test_args.t0, tb, test_args.cp)
    with raises(TypeError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, test_args.v1,
            test_args.t0, 100, test_args.cp)


def test_bad_heat_capacity(test_args: Args) -> None:
    cb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, test_args.v1,
            test_args.t0, test_args.t1, cb)
    with raises(TypeError):
        joule_thompson_law.calculate_temperature_derivative(test_args.v0, test_args.v1,
            test_args.t0, test_args.t1, 100)
