from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import work_done_by_general_force_in_one_dimension as work_law

# Description
## The force acting on a particle ranges linearly from 3 N to 5 N when
## the particle travels from a point at 1 m to a point at 2 m. The work
## of the force should amount to 4 J.


@fixture(name="test_args")
def test_args_fixture():
    f0 = Quantity(3.0 * units.newton)
    f1 = Quantity(5.0 * units.newton)
    x0 = Quantity(1.0 * units.meter)
    x1 = Quantity(2.0 * units.meter)
    Args = namedtuple("Args", "f0 f1 x0 x1")
    return Args(f0=f0, f1=f1, x0=x0, x1=x1)


def test_basic_law(test_args):
    result = work_law.calculate_work(test_args.f0, test_args.f1, test_args.x0, test_args.x1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_work = convert_to(result, units.joule).evalf(3)
    assert result_work == approx(4.0, 1e-3)


def test_bad_force(test_args):
    fb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        work_law.calculate_work(fb, test_args.f1, test_args.x0, test_args.x1)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.f1, test_args.x0, test_args.x1)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.f0, fb, test_args.x0, test_args.x1)
    with raises(TypeError):
        work_law.calculate_work(test_args.f0, 100, test_args.x0, test_args.x1)


def test_bad_position(test_args):
    xb = Quantity(1.0 * units.second)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.f0, test_args.f1, xb, test_args.x1)
    with raises(TypeError):
        work_law.calculate_work(test_args.f0, test_args.f1, 100, test_args.x1)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.f0, test_args.f1, test_args.x0, xb)
    with raises(TypeError):
        work_law.calculate_work(test_args.f0, test_args.f1, test_args.x0, 100)
