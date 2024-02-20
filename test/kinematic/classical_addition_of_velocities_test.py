from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.kinematic import classical_addition_of_velocities as galilean_law

# Description
## A kayak is moving down a river with a speed of 1 m/s relative to the water in the river.
## The water moves with a speed of 4 m/s relative to the riverbank. Therefore the kayak's speed
## relative to the riverbank is 5 m/s.

Args = namedtuple("Args", "kayak_speed_relative_to_water water_speed_relative_to_bank")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    kayak_speed_relative_to_water = Quantity(1.0 * units.meter / units.second)
    water_speed_relative_to_bank = Quantity(4.0 * units.meter / units.second)
    return Args(
        kayak_speed_relative_to_water=kayak_speed_relative_to_water,
        water_speed_relative_to_bank=water_speed_relative_to_bank,
    )


def test_law(test_args: Args) -> None:
    result = galilean_law.calculate_speed_relative_to_first_frame(
        test_args.kayak_speed_relative_to_water,
        test_args.water_speed_relative_to_bank,
    )
    assert_equal(result, 5.0 * units.meter / units.second)


def test_bad_velocities(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        galilean_law.calculate_speed_relative_to_first_frame(
            vb,
            test_args.water_speed_relative_to_bank,
        )
    with raises(TypeError):
        galilean_law.calculate_speed_relative_to_first_frame(
            100,
            test_args.water_speed_relative_to_bank,
        )
    with raises(errors.UnitsError):
        galilean_law.calculate_speed_relative_to_first_frame(
            test_args.kayak_speed_relative_to_water,
            vb,
        )
    with raises(TypeError):
        galilean_law.calculate_speed_relative_to_first_frame(
            test_args.kayak_speed_relative_to_water,
            100,
        )
