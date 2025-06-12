from collections import namedtuple
from pytest import fixture
from sympy import true, Eq
from sympy.core.numbers import Zero
from symplyphysics import units, Quantity, prefixes
from symplyphysics.conditions.electricity.maxwell_equations import (
    divergence_of_magnetic_induction_field_is_zero as divergence_cond,)

from symplyphysics.core.experimental.coordinate_systems import (CARTESIAN, CoordinateVector,
    CoordinateScalar)

## The magnetic induction vector field is given. For its distribution, the magnetic amplitude is
# known, which is equal to 1 kilotesla. Then the divergence of magnetic induction is equal to zero.

Args = namedtuple("Args", ["magnetic_induction"])


def quantity_coordinate_scalar(scalar: CoordinateScalar) -> CoordinateScalar | Zero:
    new_value = Quantity(scalar.scalar)

    if new_value.scale_factor == 0:
        return Zero()

    return CoordinateScalar(new_value, scalar.system, scalar.point)


@fixture(name="test_args")
def test_args_fixture() -> Args:
    magnetic_amplitude = Quantity(1 * prefixes.kilo * units.tesla)

    x, y, _ = CARTESIAN.base_scalars

    magnetic_induction = CoordinateVector([
        (x / Quantity(1 * units.meter)) * magnetic_amplitude,
        (-y / Quantity(1 * units.meter)) * magnetic_amplitude,
        magnetic_amplitude,
    ], CARTESIAN)

    return Args(magnetic_induction=magnetic_induction)


def test_basic_magnetic_field_divergence_condition(test_args: Args) -> None:
    result = divergence_cond.law.subs(
        divergence_cond.magnetic_flux_density(divergence_cond.position_vector),
        test_args.magnetic_induction,
    ).doit()

    for scalar in result.atoms(CoordinateScalar):
        result = result.subs(scalar, quantity_coordinate_scalar(scalar))

    assert result is true


def test_bad_condition() -> None:
    magnetic_amplitude = Quantity(1 * prefixes.kilo * units.tesla)

    x, y, _ = CARTESIAN.base_scalars

    magnetic_induction = CoordinateVector([
        (x / Quantity(1 * units.meter)) * magnetic_amplitude,
        (y / Quantity(1 * units.meter)) * magnetic_amplitude,
        magnetic_amplitude,
    ], CARTESIAN)

    result = divergence_cond.law.subs(
        divergence_cond.magnetic_flux_density(divergence_cond.position_vector),
        magnetic_induction,
    ).doit()

    for scalar in result.atoms(CoordinateScalar):
        result = result.subs(scalar, quantity_coordinate_scalar(scalar))

    assert isinstance(result, Eq)
