from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    Quantity,
    SI,
    convert_to,
    dimensionless
)

from symplyphysics.definitions import period_of_the_diffraction_grating_from_number_of_strokes as period_of_the_diffraction_grating


@fixture(name="test_args")
def test_args_fixture():
    number_of_strokes = 1000
    Args = namedtuple("Args", ["number_of_strokes"])
    return Args(
        number_of_strokes=number_of_strokes,
               )


def test_basic_period_of_the_diffraction_grating(test_args):
    result = period_of_the_diffraction_grating.calculate_diffraction_grating_period(test_args.number_of_strokes)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.length)
    result_diffraction_grating_period = convert_to(result, units.millimeter).evalf(5)
    assert result_diffraction_grating_period == approx(0.001, 0.0001)


def test_bad_number(test_args):
    bad_number = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        period_of_the_diffraction_grating.calculate_diffraction_grating_period(bad_number)
