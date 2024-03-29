from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import magnetic_induction_of_wire as induction_law

# Description
## For a distance of 0.5 meters from the conductor and a current of 20 A, the magnetic field induction is 8e-6 tesla.
## https://pressbooks.online.ucf.edu/osuniversityphysics2/chapter/magnetic-field-due-to-a-thin-straight-wire/#:~:text=Summary,the%20permeability%20of%20free%20space.

Args = namedtuple("Args", ["relative_permeability", "current", "distance"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    relative_permeability = 1
    current = Quantity(20 * units.ampere)
    distance = Quantity(0.5 * units.meter)
    return Args(relative_permeability=relative_permeability, current=current, distance=distance)


def test_basic_induction(test_args: Args) -> None:
    result = induction_law.calculate_induction(test_args.relative_permeability, test_args.current,
        test_args.distance)
    assert_equal(result, 8e-6 * units.tesla)


def test_bad_relative_permeability(test_args: Args) -> None:
    relative_permeability = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(relative_permeability, test_args.current,
            test_args.distance)


def test_bad_current(test_args: Args) -> None:
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.relative_permeability, current,
            test_args.distance)
    with raises(TypeError):
        induction_law.calculate_induction(test_args.relative_permeability, 100, test_args.distance)


def test_bad_distance(test_args: Args) -> None:
    distance = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        induction_law.calculate_induction(test_args.relative_permeability, test_args.current,
            distance)
    with raises(TypeError):
        induction_law.calculate_induction(test_args.relative_permeability, test_args.current, 100)
