from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (units, SI, convert_to, Quantity, errors)
from symplyphysics.laws.electricity import magnetic_field_intensity_of_infinite_wire as intensity_law

# Description
## With a current value of 1 ampere and a distance of 2 meters from the wire, the magnetic
## field intensity will be equal to 39.79e-3 amperes per meter.
## https://www.fxyz.ru/формулы_по_физике/электричество/магнитное_поле/правило_буравчика_электромагнетизм/напряженность_магнитного_поля/


@fixture(name="test_args")
def test_args_fixture():
    current = Quantity(0.5 * units.ampere)
    distance = Quantity(2 * units.meter)

    Args = namedtuple("Args", ["current", "distance"])
    return Args(current=current, distance=distance)


def test_basic_magnetic_intensity(test_args):
    result = intensity_law.calculate_magnetic_intensity(test_args.current, test_args.distance)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.current / units.length)
    result = convert_to(result, units.ampere / units.meter).evalf(5)
    assert result == approx(39.79e-3, rel=0.01)


def test_bad_current(test_args):
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_magnetic_intensity(current, test_args.distance)
    with raises(TypeError):
        intensity_law.calculate_magnetic_intensity(100, test_args.distance)


def test_bad_distance(test_args):
    distance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        intensity_law.calculate_magnetic_intensity(test_args.current, distance)
    with raises(TypeError):
        intensity_law.calculate_magnetic_intensity(test_args.current, 100)
