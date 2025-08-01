from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.conservation import initial_mechanical_energy_equals_final_mechanical_energy as conservation_law

Args = namedtuple("Args", ["Es"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    Es = Quantity(5 * units.joule)
    return Args(Es=Es)


def test_basic_conservation(test_args: Args) -> None:
    result = conservation_law.calculate_energy_after(test_args.Es)
    assert_equal(result, 5 * units.joule)


def test_bad_energy() -> None:
    Eb = Quantity(1 * units.meter)
    with raises(errors.UnitsError):
        conservation_law.calculate_energy_after(Eb)
    with raises(TypeError):
        conservation_law.calculate_energy_after(100)
