from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)

from symplyphysics.laws.dynamics import pressure_from_force_and_area as pressure

Args = namedtuple("Args", ["area", "force"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    force = Quantity(588 * units.newton)
    area = Quantity(0.002 * units.meter**2)
    return Args(
        area=area,
        force=force,
    )


def test_basic_pressure(test_args: Args) -> None:
    result = pressure.calculate_pressure(test_args.force, test_args.area)
    assert_equal(result, 294000 * units.pascal)


def test_bad_force(test_args: Args) -> None:
    ag = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        pressure.calculate_pressure(ag, test_args.area)
    with raises(errors.UnitsError):
        pressure.calculate_pressure(test_args.force, ag)
    with raises(TypeError):
        pressure.calculate_pressure(100, test_args.area)
    with raises(TypeError):
        pressure.calculate_pressure(test_args.force, 100)
