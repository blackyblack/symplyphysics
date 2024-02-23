from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (assert_equal, units, Quantity, errors)
from symplyphysics.laws.dynamics import coefficient_of_stiffness_from_area_length as coefficient_law

# Description
## Consider an aluminum bar. The Young's modulus for aluminum is 70e9 pascal. With a bar cross-section of 0.1 [meter^2] and a length of 3 meter,
## the stiffness coefficient is 2333333333 [newton / meter].
## https://educon.by/index.php/formuly/phystheory
## https://ru.wikipedia.org/wiki/%D0%9C%D0%BE%D0%B4%D1%83%D0%BB%D1%8C_%D0%AE%D0%BD%D0%B3%D0%B0

Args = namedtuple("Args", ["module_of_young", "area", "length"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    module_of_young = Quantity(70e9 * units.pascal)
    area = Quantity(0.1 * units.meter**2)
    length = Quantity(3 * units.meter)

    return Args(module_of_young=module_of_young,
        area=area,
        length=length)


def test_basic_resistance(test_args: Args) -> None:
    result = coefficient_law.calculate_coefficient_of_stiffness(test_args.module_of_young,
        test_args.area, test_args.length)
    assert_equal(result, 2333333333 * (units.newton / units.meter))


def test_bad_module_of_young(test_args: Args) -> None:
    module_of_young = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_of_stiffness(module_of_young, test_args.area,
            test_args.length)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_of_stiffness(100, test_args.area,
            test_args.length)


def test_bad_area(test_args: Args) -> None:
    area = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_of_stiffness(test_args.module_of_young, area,
            test_args.length)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_of_stiffness(test_args.module_of_young, 100,
            test_args.length)


def test_bad_length(test_args: Args) -> None:
    length = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        coefficient_law.calculate_coefficient_of_stiffness(test_args.module_of_young,
            test_args.area, length)
    with raises(TypeError):
        coefficient_law.calculate_coefficient_of_stiffness(test_args.module_of_young,
            test_args.area, 100)
