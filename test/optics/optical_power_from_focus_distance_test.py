from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity, prefixes,
)
from symplyphysics.laws.optics import optical_power_from_focus_distance as optical_power_law

#  Test example from https://bambookes.ru/stuff/reshenie_zadach/fizika/2-1-0-9505

Args = namedtuple("Args", ["f"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    f = Quantity(20 * prefixes.centi * units.meters)
    return Args(f=f)


def test_basic_power(test_args: Args) -> None:
    result = optical_power_law.calculate_optical_power(test_args.f)
    assert_equal(result, 5 * units.dioptre)


def test_bad_focal_distance(test_args: Args) -> None:
    fb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        optical_power_law.calculate_optical_power(fb)
    with raises(TypeError):
        optical_power_law.calculate_optical_power(100)