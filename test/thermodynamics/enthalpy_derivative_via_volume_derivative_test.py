from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import enthalpy_derivative_via_volume_derivative as enthalpy_law

# Description
## During a short isothermal quasistatic process, the volume changed from V = 1 L at T = 100 K to
## V = 1.001 L at T = 101 K. The isothermal derivative of enthalpy w.r.t. pressure in this process amounts to
## (dH/dp)_T = 92 J/atm.

Args = namedtuple("Args", "v0 v1 t0 t1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v0 = Quantity(1 * units.liter)
    v1 = Quantity(1.001 * units.liter)
    t0 = Quantity(100 * units.kelvin)
    t1 = Quantity(101 * units.kelvin)
    return Args(v0=v0, v1=v1, t0=t0, t1=t1)


def test_law(test_args: Args) -> None:
    result = enthalpy_law.calculate_enthalpy_derivative(test_args.v0, test_args.v1, test_args.t0, test_args.t1)
    assert_equal(result, 92 * units.joule / units.atmosphere, tolerance=1e-2)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_derivative(vb, test_args.v1, test_args.t0, test_args.t1)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_derivative(100, test_args.v1, test_args.t0, test_args.t1)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_derivative(test_args.v0, vb, test_args.t0, test_args.t1)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_derivative(test_args.v0, 100, test_args.t0, test_args.t1)
    

def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_derivative(test_args.v0, test_args.v1, tb, test_args.t1)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_derivative(test_args.v0, test_args.v1, 100, test_args.t1)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy_derivative(test_args.v0, test_args.v1, test_args.t0, tb)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy_derivative(test_args.v0, test_args.v1, test_args.t0, 100)

