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


@fixture
def test_args():
    E_before = Quantity(5 * units.joule)
    Args = namedtuple("Args", ["E_before"])
    return Args(E_before=E_before)


def test_basic_conservation(test_args):
    result = conservation_law.calculate_energy_after(test_args.E_before)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_ = convert_to(result, units.joule).subs(units.joule, 1).evalf(2)
    assert result_ == approx(5.0, 0.01)


def test_bad_energy():
    Eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_energy_after(Eb)
    with raises(TypeError):
        conservation_law.calculate_energy_after(100)
