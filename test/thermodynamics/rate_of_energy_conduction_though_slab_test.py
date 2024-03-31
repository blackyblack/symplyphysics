from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import rate_of_energy_conduction_through_slab as conduction_law

# Description
## A copper slab with face area A = 0.5 m**2 and thickness L = 10 cm is conducting
## heat. The temperature difference between the faces of the slab is 30°C. The thermal
## conductivity of the material is 384.1 W/(m*K). Then the rate of energy conductance
## through the slab amounts to 57.6 kW.

Args = namedtuple("Args", "k a l dt")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    k = Quantity(384.1 * units.watt / (units.meter * units.kelvin))
    a = Quantity(0.5 * units.meter**2)
    l = Quantity(10 * units.centimeter)
    dt = Quantity(30 * units.kelvin)
    return Args(k=k, a=a, l=l, dt=dt)


def test_law(test_args: Args) -> None:
    result = conduction_law.calculate_energy_conduction_rate(test_args.k, test_args.a, test_args.l, test_args.dt)
    assert_equal(result, 57.6 * prefixes.kilo * units.watt)


def test_bad_thermal_conductivity(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        conduction_law.calculate_energy_conduction_rate(kb, test_args.a, test_args.l, test_args.dt)
    with raises(TypeError):
        conduction_law.calculate_energy_conduction_rate(100, test_args.a, test_args.l, test_args.dt)


def test_bad_area(test_args: Args) -> None:
    ab = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        conduction_law.calculate_energy_conduction_rate(test_args.k, ab, test_args.l, test_args.dt)
    with raises(TypeError):
        conduction_law.calculate_energy_conduction_rate(test_args.k, 100, test_args.l, test_args.dt)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        conduction_law.calculate_energy_conduction_rate(test_args.k, test_args.a, lb, test_args.dt)
    with raises(TypeError):
        conduction_law.calculate_energy_conduction_rate(test_args.k, test_args.a, 100, test_args.dt)


def test_bad_temperature(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        conduction_law.calculate_energy_conduction_rate(test_args.k, test_args.a, test_args.l, tb)
    with raises(TypeError):
        conduction_law.calculate_energy_conduction_rate(test_args.k, test_args.a, test_args.l, 100)
