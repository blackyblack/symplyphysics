from collections import namedtuple
from pytest import fixture
from symplyphysics import (units, Quantity, prefixes)
from symplyphysics.core.fields.vector_field import VectorField
from symplyphysics.conditions.electricity.maxwell_equations import divergence_of_magnetic_induction_field_is_zero as divergence_cond

## The magnetic induction vector field is given. For its distribution, the magnetic amplitude is known, which is equal to 1 kilotesla.
## Then the divergence of magnetic induction is equal to zero.

Args = namedtuple("Args", ["magnetic_induction"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    magnetic_amplitude = Quantity(1 * prefixes.kilo * units.tesla)
    magnetic_induction = VectorField(
        lambda point: [(point.x / Quantity(1 * units.meter)) * magnetic_amplitude,
        (-point.y / Quantity(1 * units.meter)) * magnetic_amplitude, magnetic_amplitude])
    return Args(magnetic_induction=magnetic_induction)


def test_basic_magnetic_field_divergence_condition(test_args: Args) -> None:
    result = divergence_cond.magnetic_field_divergence_condition(test_args.magnetic_induction)
    assert result is True


def test_bad_condition() -> None:
    magnetic_amplitude = Quantity(1 * prefixes.kilo * units.tesla)
    magnetic_induction = VectorField(
        lambda point: [(point.x / Quantity(1 * units.meter)) * magnetic_amplitude,
        (point.y / Quantity(1 * units.meter)) * magnetic_amplitude, magnetic_amplitude])
    result = divergence_cond.magnetic_field_divergence_condition(magnetic_induction)
    assert result is False
