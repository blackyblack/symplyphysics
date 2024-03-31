from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import thermal_resistance_to_conduction as thermal_resistance_law

# Description
## The thermal resistance to conduction of a fiberglass slab (k = 0.048 W/(m*K)) of thickness
## L = 30 cm is 6.25 (m**2)*K/W.

Args = namedtuple("Args", "l k")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    l = Quantity(30 * units.centimeter)
    k = Quantity(0.048 * units.watt / units.meter / units.kelvin)
    return Args(l=l, k=k)


def test_law(test_args: Args) -> None:
    result = thermal_resistance_law.calculate_thermal_resistance(test_args.l, test_args.k)
    assert_equal(result, 6.25 * units.meter**2 * units.kelvin / units.watt)


def test_bad_length(test_args: Args) -> None:
    lb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        thermal_resistance_law.calculate_thermal_resistance(lb, test_args.k)
    with raises(TypeError):
        thermal_resistance_law.calculate_thermal_resistance(100, test_args.k)


def test_bad_thermal_conductivity(test_args: Args) -> None:
    kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        thermal_resistance_law.calculate_thermal_resistance(test_args.l, kb)
    with raises(TypeError):
        thermal_resistance_law.calculate_thermal_resistance(test_args.l, 100)
