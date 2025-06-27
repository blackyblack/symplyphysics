from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.definitions import electrical_conductance_is_inverse_resistance as conductance_def

# Description
## If the object has 2 Ohm resistance, it should have 0.5 Siemens conductance. No external calculators were used for such computation.

Args = namedtuple("Args", ["R"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    R = Quantity(2 * units.ohm)
    return Args(R=R)


def test_basic_conductance(test_args: Args) -> None:
    result = conductance_def.calculate_conductance(test_args.R)
    assert_equal(result, 0.5 * units.siemens)


def test_bad_resistance() -> None:
    Rb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conductance_def.calculate_conductance(Rb)
    with raises(TypeError):
        conductance_def.calculate_conductance(100)
