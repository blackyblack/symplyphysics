from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.core.symbols.probability import Probability
from symplyphysics.laws.nuclear import thermal_non_leakage_probability_from_diffusion_length as non_leakage_factor

Args = namedtuple("Args", ["Lth", "Bg"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # water diffusion area is 8.1 cm^2
    thermal_diffusion_area = Quantity(8.1 * units.centimeter**2)
    # sphere with radius = 1 meter
    geometric_buckling = Quantity(9.869 / units.meter**2)
    return Args(Lth=thermal_diffusion_area, Bg=geometric_buckling)


def test_basic_non_leakage_factor(test_args: Args) -> None:
    result = non_leakage_factor.calculate_probability(test_args.Lth, test_args.Bg)
    assert isinstance(result, Probability)
    assert_equal(result, 0.9921)


def test_bad_diffusion_area(test_args: Args) -> None:
    Lb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(Lb, test_args.Bg)
    with raises(TypeError):
        non_leakage_factor.calculate_probability(100, test_args.Bg)


def test_bad_buckling(test_args: Args) -> None:
    Bb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        non_leakage_factor.calculate_probability(test_args.Lth, Bb)
    with raises(TypeError):
        non_leakage_factor.calculate_probability(test_args.Lth, 100)
