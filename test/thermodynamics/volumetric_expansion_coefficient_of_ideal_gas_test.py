from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import (
    volumetric_expansion_coefficient_of_ideal_gas as expansion_law,
)

# Description
## The volumetric expansion coefficient of an ideal gas at temperature T = 100 K is 0.01 1/K.

Args = namedtuple("Args", "t")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    t = Quantity(100 * units.kelvin)
    return Args(t=t)


def test_law(test_args: Args) -> None:
    result = expansion_law.calculate_volumetric_expansion_coefficient(test_args.t)
    assert_equal(result, 0.01 / units.kelvin)


def test_bad_temperature() -> None:
    tb = Quantity(1.0 * units.meter)
    with raises(errors.UnitsError):
        expansion_law.calculate_volumetric_expansion_coefficient(tb)
    with raises(TypeError):
        expansion_law.calculate_volumetric_expansion_coefficient(100)
