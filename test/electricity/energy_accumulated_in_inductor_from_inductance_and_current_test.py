from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.electricity import energy_accumulated_in_inductor_from_inductance_and_current as inductor_law

# Description
## Assert we have 150mH inductor with 0.5A current flowing through it.
## According to law we should have amount of energy accumulated in this inductor equals to 0.150 * 0.5**2 / 2 = 0.01875 Joules.


@fixture
def test_args():
    Inductance = Quantity(150 * units.milli * units.henry)
    Current = Quantity(0.5 * units.ampere)
    Args = namedtuple("Args", ["Inductance", "Current"])
    return Args(Inductance=Inductance, Current=Current)


def test_basic_energy(test_args):
    result = inductor_law.calculate_accumulated_energy(test_args.Inductance, test_args.Current)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_power = convert_to(result, units.joule).subs(units.joule, 1).evalf(5)
    assert result_power == approx(0.01875, 0.00001)


def test_bad_Inductance(test_args):
    Lb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductor_law.calculate_accumulated_energy(Lb, test_args.Current)
    with raises(TypeError):
        inductor_law.calculate_accumulated_energy(100, test_args.Current)


def test_bad_Current(test_args):
    Ib = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        inductor_law.calculate_accumulated_energy(test_args.Inductance, Ib)
    with raises(TypeError):
        inductor_law.calculate_accumulated_energy(test_args.Inductance, 100)
