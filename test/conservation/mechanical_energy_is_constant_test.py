from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.conservation import mechanical_energy_is_constant as conservation_law


@fixture(name="test_args")
def test_args_fixture():
    Es = Quantity(5 * units.joule)
    Args = namedtuple("Args", ["Es"])
    return Args(Es=Es)


def test_basic_conservation(test_args):
    result = conservation_law.calculate_energy_after(test_args.Es)
    assert_equal(result, 5 * units.joule)


def test_bad_energy():
    Eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_energy_after(Eb)
    with raises(TypeError):
        conservation_law.calculate_energy_after(100)
