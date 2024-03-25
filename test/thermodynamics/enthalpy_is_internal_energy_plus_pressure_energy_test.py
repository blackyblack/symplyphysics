from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
    prefixes,
)
from symplyphysics.laws.thermodynamics import (
    enthalpy_is_internal_energy_plus_pressure_energy as enthalpy_law,
)

# Description
## The internal energy of a system is 1 kJ, its internal pressure is 1 kPa and its volume
## is 0.3 m**3. The enthalpy of the system is 1.3 kJ.

Args = namedtuple("Args", "u p v")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    u = Quantity(1 * prefixes.kilo * units.joule)
    p = Quantity(1 * prefixes.kilo * units.pascal)
    v = Quantity(0.3 * units.meter**3)
    return Args(u=u, p=p, v=v)


def test_law(test_args: Args) -> None:
    result = enthalpy_law.calculate_enthalpy(test_args.u, test_args.p, test_args.v)
    assert_equal(result, 1.3 * prefixes.kilo * units.joule)


def test_bad_energy(test_args: Args) -> None:
    ub = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy(ub, test_args.p, test_args.v)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy(100, test_args.p, test_args.v)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy(test_args.u, pb, test_args.v)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy(test_args.u, 100, test_args.v)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        enthalpy_law.calculate_enthalpy(test_args.u, test_args.p, vb)
    with raises(TypeError):
        enthalpy_law.calculate_enthalpy(test_args.u, test_args.p, 100)
