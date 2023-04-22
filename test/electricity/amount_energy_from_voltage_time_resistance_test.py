from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import amount_energy_from_voltage_time_resistance as joule_lenz_law

#  How much energy can a household electric kettle produce to heat water in one minute
#  at 220 volts and a heater resistance of 36 ohms?


@fixture
def test_args():
    U = Quantity(220 * units.volt)
    t = Quantity(60 * units.second)
    R = Quantity(36 * units.ohm)
    Args = namedtuple("Args", ["U", "t", "R"])
    return Args(U=U, t=t, R=R)


def test_basic_amount(test_args):
    result = joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t,
                                                    test_args.R)
    assert SI.get_dimension_system().equivalent_dims(result.dimension,
                                                     units.energy)
    result_energy = convert_to(result, units.joule).subs(units.joule,
                                                         1).evalf(6)
    assert result_energy == approx(80666.6, 0.000001)


def test_bad_voltage(test_args):
    Ub = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_amount_energy(Ub, test_args.t, test_args.R)
    with raises(TypeError):
        joule_lenz_law.calculate_amount_energy(100, test_args.t, test_args.R)


def test_bad_time(test_args):
    tb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_amount_energy(test_args.U, tb, test_args.R)
    with raises(TypeError):
        joule_lenz_law.calculate_amount_energy(test_args.U, 100, test_args.R)


def test_bad_resistance(test_args):
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, Rb)
    with raises(TypeError):
        joule_lenz_law.calculate_amount_energy(test_args.U, test_args.t, 100)
