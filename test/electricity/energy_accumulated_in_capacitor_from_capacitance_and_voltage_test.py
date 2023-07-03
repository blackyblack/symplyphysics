from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import energy_accumulated_in_capacitor_from_capacitance_and_voltage as capacitor_law

# Description
## Assert we have 220 uF capacitor charged to 10 volts.
## According to law we should have amount of energy accumulated in this capacitor equals to 220 * 0.000001 * 10**2 / 2 = 0.011 Joules.


@fixture(name="test_args")
def test_args_fixture():
    C = Quantity(220 * units.micro * units.farad)
    V = Quantity(10 * units.volt)
    Args = namedtuple("Args", ["C", "V"])
    return Args(C=C, V=V)


def test_basic_energy(test_args):
    result = capacitor_law.calculate_accumulated_energy(test_args.C, test_args.V)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_power = convert_to(result, units.joule).evalf(5)
    assert result_power == approx(0.011, 0.00001)


def test_bad_capacitance(test_args):
    Cb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitor_law.calculate_accumulated_energy(Cb, test_args.V)
    with raises(TypeError):
        capacitor_law.calculate_accumulated_energy(100, test_args.V)


def test_bad_voltage(test_args):
    Vb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        capacitor_law.calculate_accumulated_energy(test_args.C, Vb)
    with raises(TypeError):
        capacitor_law.calculate_accumulated_energy(test_args.C, 100)
