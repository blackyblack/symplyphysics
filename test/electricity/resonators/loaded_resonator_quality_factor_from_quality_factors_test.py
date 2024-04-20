from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity.resonators import loaded_resonator_quality_factor_from_quality_factors as factor_law

# Description
## Resonator's own quality factor is 8000, the quality factor of the external circuit is 10000.
## Then the loaded quality factor of the resonator will be 4444.

Args = namedtuple("Args", ["resonator_quality_factor", "external_circuit_quality_factor"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    resonator_quality_factor = 8000
    external_circuit_quality_factor = 10000

    return Args(resonator_quality_factor=resonator_quality_factor,
        external_circuit_quality_factor=external_circuit_quality_factor)


def test_basic_loaded_resonator_quality_factor(test_args: Args) -> None:
    result = factor_law.calculate_quality_factor(test_args.resonator_quality_factor,
        test_args.external_circuit_quality_factor)
    assert_equal(result, 4444)


def test_bad_quality_factors(test_args: Args) -> None:
    bad_quality_factor = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(bad_quality_factor,
            test_args.external_circuit_quality_factor)
    with raises(errors.UnitsError):
        factor_law.calculate_quality_factor(test_args.resonator_quality_factor, bad_quality_factor)
