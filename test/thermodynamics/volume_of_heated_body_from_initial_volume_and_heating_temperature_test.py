from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
)

from symplyphysics.laws.thermodynamics import volume_of_heated_body_from_initial_volume_and_heating_temperature as heating


@fixture(name="test_args")
def test_args_fixture():
    start_volume = Quantity(20 * units.liter)
    expansion_coefficient = Quantity(0.001 * (1 / units.kelvin))
    start_temperature = Quantity(283.15 * units.kelvin)
    finish_temperature = Quantity(303.15 * units.kelvin)
    Args = namedtuple("Args",
        ["start_volume", "expansion_coefficient", "finish_temperature", "start_temperature"])
    return Args(start_volume=start_volume,
        expansion_coefficient=expansion_coefficient,
        finish_temperature=finish_temperature,
        start_temperature=start_temperature)


def test_basic_volume(test_args):
    result = heating.calculate_final_volume(test_args.start_volume, test_args.expansion_coefficient,
        test_args.finish_temperature, test_args.start_temperature)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.volume)
    result_volume = convert_to(result, units.liter).evalf(5)
    assert result_volume == approx(20.4, 0.001)


def test_bad_coefficient(test_args):
    bad_coefficient = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heating.calculate_final_volume(test_args.start_volume, bad_coefficient,
            test_args.finish_temperature, test_args.start_temperature)
    with raises(TypeError):
        heating.calculate_final_volume(test_args.start_volume, 100, test_args.finish_temperature,
            test_args.start_temperature)


def test_bad_temperature(test_args):
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


def test_bad_volume(test_args):
    bad_volume = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        heating.calculate_final_volume(bad_volume, test_args.expansion_coefficient,
            test_args.finish_temperature, test_args.start_temperature)
    with raises(TypeError):
        heating.calculate_final_volume(100, test_args.expansion_coefficient,
            test_args.finish_temperature, test_args.start_temperature)
