from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.conservation import mechanical_energy_is_constant as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    Es = Quantity(5 * units.joule)
    Args = namedtuple("Args", ["Es"])
    return Args(Es=Es)


def test_basic_conservation(test_args):
    result = conservation_law.calculate_energy_after(test_args.Es)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_ = convert_to(result, units.joule).subs(units.joule, 1).evalf(2)
    assert result_ == approx(5.0, 0.01)


def test_bad_energy():
    Eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_energy_after(Eb)
    with raises(AttributeError):
        conservation_law.calculate_energy_after(100)
