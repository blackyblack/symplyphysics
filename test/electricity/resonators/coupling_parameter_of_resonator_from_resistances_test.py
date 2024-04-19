from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.resonators import coupling_parameter_of_resonator_from_resistances as parameter_law

# Description
## Resonator's resistance is 80 ohms, the load resistance is 100 ohms.
## Then the coupling parameter will be 0.8.

Args = namedtuple("Args", ["resonator_resistance", "load_resistance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resonator_resistance = Quantity(80 * units.ohm)
    load_resistance = Quantity(100 * units.ohm)

    return Args(resonator_resistance=resonator_resistance, load_resistance=load_resistance)


def test_basic_coupling_parameter(test_args: Args) -> None:
    result = parameter_law.calculate_coupling_parameter(test_args.resonator_resistance,
        test_args.load_resistance)
    assert_equal(result, 0.8)


def test_bad_resistances(test_args: Args) -> None:
    bad_resistance = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        parameter_law.calculate_coupling_parameter(bad_resistance,
            test_args.load_resistance)
    with raises(errors.UnitsError):
        parameter_law.calculate_coupling_parameter(test_args.resonator_resistance,
            bad_resistance)
