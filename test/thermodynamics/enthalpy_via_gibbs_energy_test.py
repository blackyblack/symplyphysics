from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import enthalpy_via_gibbs_energy as gibbs_helmholtz_law

# Description
## A thermodynamic system is at temperature T = 300 K. A thermodynamic process happens around this temperature
## point, in which the value of the Gibbs energy is G = 100 J at the beginning (T = 299 K) and G = 110 J at the
## end (T = 301 K). The value of the enthalpy at T = 300 K is H = -1.4 kJ.

Args = namedtuple("Args", "g0 g1 t0 t1 t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    g0 = Quantity(100 * units.joule)
    g1 = Quantity(110 * units.joule)
    t0 = Quantity(299 * units.kelvin)
    t1 = Quantity(301 * units.kelvin)
    t = Quantity(300 * units.kelvin)
    return Args(g0=g0, g1=g1, t0=t0, t1=t1, t=t)


def test_law(test_args: Args) -> None:
    result = gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, test_args.g1, test_args.t0,
        test_args.t1, test_args.t)
    assert_equal(result, -1.4 * prefixes.kilo * units.joule, relative_tolerance=4e-3)


def test_bad_energy(test_args: Args) -> None:
    gb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_helmholtz_law.calculate_enthalpy(gb, test_args.g1, test_args.t0, test_args.t1,
            test_args.t)
    with raises(TypeError):
        gibbs_helmholtz_law.calculate_enthalpy(100, test_args.g1, test_args.t0, test_args.t1,
            test_args.t)
    with raises(errors.UnitsError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, gb, test_args.t0, test_args.t1,
            test_args.t)
    with raises(TypeError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, 100, test_args.t0, test_args.t1,
            test_args.t)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, test_args.g1, tb, test_args.t1,
            test_args.t)
    with raises(TypeError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, test_args.g1, 100, test_args.t1,
            test_args.t)
    with raises(errors.UnitsError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, test_args.g1, test_args.t0, tb,
            test_args.t)
    with raises(TypeError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, test_args.g1, test_args.t0, 100,
            test_args.t)
    with raises(errors.UnitsError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, test_args.g1, test_args.t0,
            test_args.t1, tb)
    with raises(TypeError):
        gibbs_helmholtz_law.calculate_enthalpy(test_args.g0, test_args.g1, test_args.t0,
            test_args.t1, 100)
