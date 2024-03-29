from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import volume_of_heated_body_from_initial_volume_and_heating_temperature as heating

Args = namedtuple("Args",
    ["start_volume", "expansion_coefficient", "finish_temperature", "start_temperature"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    start_volume = Quantity(20 * units.liter)
    expansion_coefficient = Quantity(0.001 * (1 / units.kelvin))
    start_temperature = Quantity(283.15 * units.kelvin)
    finish_temperature = Quantity(303.15 * units.kelvin)
    return Args(start_volume=start_volume,
        expansion_coefficient=expansion_coefficient,
        finish_temperature=finish_temperature,
        start_temperature=start_temperature)


def test_basic_volume(test_args: Args) -> None:
    result = heating.calculate_final_volume(test_args.start_volume, test_args.expansion_coefficient,
        test_args.finish_temperature, test_args.start_temperature)
    assert_equal(result, 20.4 * units.liter)


def test_bad_coefficient(test_args: Args) -> None:
    bad_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heating.calculate_final_volume(test_args.start_volume, bad_coefficient,
            test_args.finish_temperature, test_args.start_temperature)
    with raises(TypeError):
        heating.calculate_final_volume(test_args.start_volume, 100, test_args.finish_temperature,
            test_args.start_temperature)


def test_bad_temperature(test_args: Args) -> None:
    bad_temperature = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heating.calculate_final_volume(test_args.start_volume, test_args.expansion_coefficient,
            bad_temperature, test_args.start_temperature)
    with raises(errors.UnitsError):
        heating.calculate_final_volume(test_args.start_volume, test_args.expansion_coefficient,
            test_args.finish_temperature, bad_temperature)
    with raises(TypeError):
        heating.calculate_final_volume(test_args.start_volume, test_args.expansion_coefficient, 100,
            test_args.start_temperature)
    with raises(TypeError):
        heating.calculate_final_volume(test_args.start_volume, test_args.expansion_coefficient,
            test_args.finish_temperature, 100)


def test_bad_volume(test_args: Args) -> None:
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heating.calculate_final_volume(bad_volume, test_args.expansion_coefficient,
            test_args.finish_temperature, test_args.start_temperature)
    with raises(TypeError):
        heating.calculate_final_volume(100, test_args.expansion_coefficient,
            test_args.finish_temperature, test_args.start_temperature)
