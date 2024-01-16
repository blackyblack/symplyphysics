from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.conservation import mass_after_equals_to_mass_before as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    ms = Quantity(5 * units.kilograms)
    Args = namedtuple("Args", ["ms"])
    return Args(ms=ms)


def test_basic_conservation(test_args):
    result = conservation_law.calculate_mass_after(test_args.ms)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.mass)
    result_ = convert_to(result, units.kilograms).evalf(2)
    assert result_ == approx(5.0, 0.01)


def test_bad_mass():
    mb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_mass_after(mb)
    with raises(TypeError):
        conservation_law.calculate_mass_after(100)
