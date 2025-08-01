from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import errors, units, Quantity
from symplyphysics.laws.kinematics.vector import center_of_mass_for_system_of_particles as com_def

from symplyphysics.core.experimental.coordinate_systems import CARTESIAN, QuantityCoordinateVector
from symplyphysics.core.experimental.approx import assert_equal_vectors

# Description
## Particle 1 weighs 1 kg and has position (0, -1, 2) m.
## Particle 2 weighs 4 kg and has position (2, 3, 1) m.
## Particle 3 weighs 3 kg and has position (-4, 1, -2) m.
## The com of the system consisting of particles 1 and 2 has position (1.6, 2.2, 1.2) m.
## The com of the system consisting of particles 1, 2, and 3 has position (-0.5, 1.75, 0.0) m.

Args = namedtuple("Args", "m1 m2 m3 r1 r2 r3")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    m1 = Quantity(1.0 * units.kilogram)
    m2 = Quantity(4.0 * units.kilogram)
    m3 = Quantity(3.0 * units.kilogram)
    r1 = QuantityCoordinateVector([
        Quantity(0.0 * units.meter),
        Quantity(-1.0 * units.meter),
        Quantity(2.0 * units.meter),
    ], CARTESIAN)
    r2 = QuantityCoordinateVector([
        Quantity(2.0 * units.meter),
        Quantity(3.0 * units.meter),
        Quantity(1.0 * units.meter),
    ], CARTESIAN)
    r3 = QuantityCoordinateVector([
        Quantity(-4.0 * units.meter),
        Quantity(1.0 * units.meter),
        Quantity(-2.0 * units.meter),
    ], CARTESIAN)
    return Args(m1=m1, m2=m2, m3=m3, r1=r1, r2=r2, r3=r3)


def test_basic_law_two_particles(test_args: Args) -> None:
    result = com_def.calculate_center_of_mass(
        [test_args.m1, test_args.m2],
        [test_args.r1, test_args.r2],
    )
    expected = QuantityCoordinateVector([
        1.6 * units.meter,
        2.2 * units.meter,
        1.2 * units.meter,
    ], CARTESIAN)
    assert_equal_vectors(result, expected)


def test_basic_law_three_particles(test_args: Args) -> None:
    result = com_def.calculate_center_of_mass(
        [test_args.m1, test_args.m2, test_args.m3],
        [test_args.r1, test_args.r2, test_args.r3],
    )
    expected = QuantityCoordinateVector([
        -0.5 * units.meter,
        1.75 * units.meter,
        0,
    ], CARTESIAN)
    assert_equal_vectors(result, expected)


def test_bad_masses(test_args: Args) -> None:
    mb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        com_def.calculate_center_of_mass([mb], [test_args.r1])
    with raises(errors.UnitsError):
        com_def.calculate_center_of_mass(mb, [test_args.r1])

    with raises(errors.UnitsError):
        com_def.calculate_center_of_mass([mb, test_args.r2], [test_args.r1, test_args.r2])
    with raises(ValueError):
        com_def.calculate_center_of_mass([test_args.m1, mb], [test_args.r1, test_args.r2])

    with raises(TypeError):
        com_def.calculate_center_of_mass(test_args.m1, [test_args.r1])
    with raises(TypeError):
        com_def.calculate_center_of_mass(100, [test_args.r1])
    with raises(TypeError):
        com_def.calculate_center_of_mass([100], [test_args.r1])
    with raises(ValueError):
        com_def.calculate_center_of_mass([], [test_args.r1])
    with raises(ValueError):
        com_def.calculate_center_of_mass([test_args.m1, test_args.m1], [test_args.r1])


def test_bad_positions(test_args: Args) -> None:
    rb = QuantityCoordinateVector([
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
        Quantity(1.0 * units.coulomb),
    ], CARTESIAN)
    with raises(errors.UnitsError):
        com_def.calculate_center_of_mass([test_args.m1], [rb])
    with raises(errors.UnitsError):
        com_def.calculate_center_of_mass([test_args.m1], rb)

    with raises(ValueError):
        com_def.calculate_center_of_mass([test_args.m1, test_args.m2], [rb, test_args.r1])
    with raises(ValueError):
        com_def.calculate_center_of_mass([test_args.m1, test_args.m2], [test_args.r1, rb])

    with raises(TypeError):
        com_def.calculate_center_of_mass([test_args.m1], test_args.r1)
    with raises(TypeError):
        com_def.calculate_center_of_mass([test_args.m1], 100)
    with raises(TypeError):
        com_def.calculate_center_of_mass([test_args.m1], [100])
    with raises(ValueError):
        com_def.calculate_center_of_mass([test_args.m1], [])
