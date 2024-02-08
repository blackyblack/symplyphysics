from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.electricity import mass_of_the_substance_deposited_on_electrode as mass_law

# Description
## Consider copper with an electrochemical equivalent of 1.19 [gram / (ampere * hour)]. Let the current be 2 ampere and the time
## be 1 hour. Then the mass of the substance deposited on the electrode will be 2.38 gram.
## https://www.indigomath.ru//raschety/2JNq0p.html


@fixture(name="test_args")
def test_args_fixture():
    equivalent = Quantity(1.19 * units.gram / (units.ampere * units.hour))
    current = Quantity(2 * units.ampere)
    time = Quantity(1 * units.hour)

    Args = namedtuple("Args", ["equivalent", "current", "time"])
    return Args(equivalent=equivalent, current=current, time=time)


def test_basic_mass(test_args):
    result = mass_law.calculate_mass(test_args.equivalent, test_args.current, test_args.time)
    assert_equal(result, 2.38 * units.gram)


def test_bad_equivalent(test_args):
    equivalent = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mass_law.calculate_mass(equivalent, test_args.current, test_args.time)
    with raises(TypeError):
        mass_law.calculate_mass(100, test_args.current, test_args.time)


def test_bad_current(test_args):
    current = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        mass_law.calculate_mass(test_args.equivalent, current, test_args.time)
    with raises(TypeError):
        mass_law.calculate_mass(test_args.equivalent, 100, test_args.time)


def test_bad_time(test_args):
    time = Quantity(1 * units.kelvin)
    with raises(errors.UnitsError):
        mass_law.calculate_mass(test_args.equivalent, test_args.current, time)
    with raises(TypeError):
        mass_law.calculate_mass(test_args.equivalent, test_args.current, 100)
