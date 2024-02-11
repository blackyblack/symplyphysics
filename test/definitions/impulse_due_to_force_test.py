from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.definitions import impulse_due_to_force as impulse_def

# Description
## During the collision the force was constant and equal to 200 N. The collision started
## at zero time and lasted 3 ms. The impulse due to force in collision was 0.6 kg*m/s.

Args = namedtuple("Args", "F t_i t_f")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    F = Quantity(200.0 * units.newton)
    t_i = Quantity(0)
    t_f = Quantity(3.0 * prefixes.milli * units.second)
    return Args(F=F, t_i=t_i, t_f=t_f)


def test_definition(test_args: Args) -> None:
    result = impulse_def.calculate_impulse(test_args.F, test_args.F, test_args.t_i, test_args.t_f)
    assert_equal(result, 0.6 * units.kilogram * units.meter / units.second)


def test_bad_forces(test_args: Args) -> None:
    Fb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        impulse_def.calculate_impulse(Fb, test_args.F, test_args.t_i, test_args.t_f)
    with raises(TypeError):
        impulse_def.calculate_impulse(100, test_args.F, test_args.t_i, test_args.t_f)
    with raises(errors.UnitsError):
        impulse_def.calculate_impulse(test_args.F, Fb, test_args.t_i, test_args.t_f)
    with raises(TypeError):
        impulse_def.calculate_impulse(test_args.F, 100, test_args.t_i, test_args.t_f)


def test_bad_times(test_args: Args) -> None:
    tb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        impulse_def.calculate_impulse(test_args.F, test_args.F, tb, test_args.t_f)
    with raises(TypeError):
        impulse_def.calculate_impulse(test_args.F, test_args.F, 100, test_args.t_f)
    with raises(errors.UnitsError):
        impulse_def.calculate_impulse(test_args.F, test_args.F, test_args.t_i, tb)
    with raises(TypeError):
        impulse_def.calculate_impulse(test_args.F, test_args.F, test_args.t_i, 100)
