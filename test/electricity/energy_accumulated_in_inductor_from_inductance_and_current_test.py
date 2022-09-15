# Description
## Assert we have 150mH  inductor with 0.5A current flowing through it.
## According to law we should have amount of energy accumulated in this inductor equals to 0.150 * 0.5**2 / 2 = 0.01875 Joules.

from collections import namedtuple
from pytest import approx, fixture, raises

from symplyphysics import (
    units, convert_to, SI, errors
)
from symplyphysics.laws.electricity import energy_accumulated_in_inductor_from_inductance_and_current as inductor_law

@fixture
def test_args():
    Inductance = units.Quantity('Inductance')
    Current = units.Quantity('Current')
    SI.set_quantity_dimension(Inductance, units.inductance)
    SI.set_quantity_scale_factor(Inductance, 0.15 * units.henry)
    SI.set_quantity_dimension(Current, units.current)
    SI.set_quantity_scale_factor(Current, 0.5 * units.ampere)

    Args = namedtuple('Args', ['Inductance', 'Current'])
    return Args(Inductance = Inductance, Current = Current)

def test_basic_energy(test_args):
    result = inductor_law.calculate_accumulated_energy(test_args.Inductance, test_args.Current)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_power = convert_to(result, units.joule).subs(units.joule, 1).evalf(2)
    assert result_power == approx(0.01875, 0.01)

def test_bad_Inductance(test_args):
    Lb = units.Quantity('Lb')
    SI.set_quantity_dimension(Lb, units.length)
    SI.set_quantity_scale_factor(Lb, 1 * units.meter)

    with raises(errors.UnitsError):
        inductor_law.calculate_accumulated_energy(Lb, test_args.Current)

    with raises(TypeError):
        inductor_law.calculate_accumulated_energy(100, test_args.Current)

def test_bad_Current(test_args):
    Ib = units.Quantity('Ib')
    SI.set_quantity_dimension(Ib, units.length)
    SI.set_quantity_scale_factor(Ib, 1 * units.meter)

    with raises(errors.UnitsError):
        inductor_law.calculate_accumulated_energy(test_args.Inductance, Ib)

    with raises(TypeError):
        inductor_law.calculate_accumulated_energy(test_args.Inductance, 100)
        