from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors, prefixes)
from symplyphysics.laws.optics import light_pressure as pressure_law

# Description
## Let the intensity of the incident radiation be 0.6 [watt / meter^2] and the reflection coefficient be 1.
## Then the light pressure will be 4 nanopascal.
## https://portal.tpu.ru/SHARED/m/MTN/academic/Tab1/physics/Tab1/Tab/Tab/11%20Davlenie%20sveta.pdf
## page 10

Args = namedtuple("Args", ["intensity", "reflection_coefficient"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    intensity = Quantity(0.6 * (units.watt / units.meter**2))
    reflection_coefficient = 1
    return Args(intensity=intensity, reflection_coefficient=reflection_coefficient)


def test_basic_pressure(test_args: Args) -> None:
    result = pressure_law.calculate_pressure(test_args.intensity, test_args.reflection_coefficient)
    assert_equal(result, 4 * prefixes.nano * units.pascal)


def test_bad_intensity(test_args: Args) -> None:
    intensity = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure(intensity, test_args.reflection_coefficient)
    with raises(TypeError):
        pressure_law.calculate_pressure(100, test_args.reflection_coefficient)


def test_bad_reflection_coefficient(test_args: Args) -> None:
    reflection_coefficient = Quantity(1 * units.joule)
    with raises(errors.UnitsError):
        pressure_law.calculate_pressure(test_args.intensity, reflection_coefficient)
