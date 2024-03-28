from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    free_energy_change_via_temperature_change_and_volume_change as free_energy_law,
)

# Description
## A closed system in thermal equlibrium has entropy S = 1 J/k and pressure p = 1 Pa inside.
## The temperature of the system drops by 0.01 K and the volume grows by 0.02 m**3. Then
## the change in Helmholtz free energy of the system is -0.01 J.

Args = namedtuple("Args", "s dt p dv")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    s = Quantity(1 * units.joule / units.kelvin)
    dt = Quantity(-0.01 * units.kelvin)
    p = Quantity(1 * units.pascal)
    dv = Quantity(0.02 * units.meter**3)
    return Args(s=s, dt=dt, p=p, dv=dv)


def test_law(test_args: Args) -> None:
    result = free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p, test_args.dv)
    assert_equal(result, -0.01 * units.joule)


def test_bad_entropy(test_args: Args) -> None:
    sb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(sb, test_args.dt, test_args.p, test_args.dv)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(100, test_args.dt, test_args.p, test_args.dv)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, tb, test_args.p, test_args.dv)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(test_args.s, 100, test_args.p, test_args.dv)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, pb, test_args.dv)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, 100, test_args.dv)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p, vb)
    with raises(TypeError):
        free_energy_law.calculate_free_energy_change(test_args.s, test_args.dt, test_args.p, 100)
