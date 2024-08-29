from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal_vectors,
    units,
    errors,
    Quantity,
    QuantityVector,
)
from symplyphysics.laws.electricity.vector import electric_dipole_moment_is_charge_times_displacement as dipole_moment

# Description
## The vector of the dipole moment of two point charges of magnitude Q = 2e-13 C, the displacement
## vector between which is [1, -1, 0.5] m, amounts to [2e-13, -2e-13, 1e-13] C*m

Args = namedtuple("Args", "q l")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    q = Quantity(2e-13 * units.coulomb)
    l = QuantityVector([
        Quantity(1.0 * units.meter),
        Quantity(-1.0 * units.meter),
        Quantity(0.5 * units.meter),
    ])
    return Args(q=q, l=l)


def test_basic_law(test_args: Args) -> None:
    result = dipole_moment.calculate_dipole_moment(test_args.q, test_args.l)
    unit = units.coulomb * units.meter
    assert_equal_vectors(
        result,
        QuantityVector([2e-13 * unit, -2e-13 * unit, 1e-13 * unit]),
    )


def test_bad_charge(test_args: Args) -> None:
    qb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        dipole_moment.calculate_dipole_moment(qb, test_args.l)
    with raises(TypeError):
        dipole_moment.calculate_dipole_moment(100, test_args.l)


def test_bad_displacement(test_args: Args) -> None:
    lb = Quantity(1 * units.second)
    with raises(errors.UnitsError):
        dipole_moment.calculate_dipole_moment(test_args.q, lb)
    with raises(TypeError):
        dipole_moment.calculate_dipole_moment(test_args.q, 100)
