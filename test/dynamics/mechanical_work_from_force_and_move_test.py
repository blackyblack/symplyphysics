from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    units, convert_to, SI, errors
)
from sympy.physics.units.definitions.dimension_definitions import angle as angle_type
from symplyphysics.core.symbols.quantities import Quantity
from symplyphysics.laws.dynamics import mechanical_work_from_force_and_move as work_law

# Description
## Force of 100N is applied to heavy object lying oh horizontal table. Force is directed 60 degrees up.
## Object slides on the table with no friction and is been moved to 3m far.
## As cosine of 60 degrees is 1/2, work should be 100 * 3 / 2 = 150 Joules.

@fixture
def test_args():
    F = Quantity(100 * units.newton)
    S = Quantity(3 * units.meter) 
    angle_between_force_and_axis = Quantity(60 * units.degree)
    angle_between_movement_and_axis = Quantity(0 * units.degree)
    Args = namedtuple("Args", ["F", "S", "Fa", "Sa"])
    return Args(F=F, S=S, Fa=angle_between_force_and_axis, Sa=angle_between_movement_and_axis)

def test_basic_work(test_args):
    result = work_law.calculate_work(test_args.F, test_args.S, test_args.Fa, test_args.Sa)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).subs(units.joule, 1).evalf(4)
    assert result_work == approx(150, 0.01)

def test_work_with_number(test_args):
    work_law.calculate_work(test_args.F, test_args.S, 100, test_args.Sa)
    work_law.calculate_work(test_args.F, test_args.S, test_args.Fa, 100)

def test_bad_force(test_args):
    Fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_work(Fb, test_args.S, test_args.Fa, test_args.Sa)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.S, test_args.Fa, test_args.Sa)

def test_bad_move(test_args):
    Sb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.F, Sb, test_args.Fa, test_args.Sa)
    with raises(TypeError):
        work_law.calculate_work(test_args.F, 100, test_args.Fa, test_args.Sa)

def test_bad_angle(test_args):
    ab = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.F, test_args.S, ab, test_args.Sa)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.F, test_args.S, test_args.Fa, ab)
