from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, errors, units, Quantity
from symplyphysics.laws.electricity import voltage_is_line_integral_of_electric_field as law

Args = namedtuple("Args", "e0 e1 s0 s1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e0 = Quantity(3 * units.volt / units.meter)
    e1 = Quantity(7 * units.volt / units.meter)
    s0 = Quantity(0.1 * units.meter)
    s1 = Quantity(0.2 * units.meter)
    return Args(e0=e0, e1=e1, s0=s0, s1=s1)


def test_law(test_args: Args) -> None:
    result = law.calculate_voltage(test_args.e0, test_args.e1, test_args.s0, test_args.s1)
    assert_equal(result, -0.5 * units.volt)


def test_bad_electric_field(test_args: Args) -> None:
    eb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_voltage(eb, test_args.e1, test_args.s0, test_args.s1)
    with raises(TypeError):
        law.calculate_voltage(100, test_args.e1, test_args.s0, test_args.s1)
    with raises(errors.UnitsError):
        law.calculate_voltage(test_args.e0, eb, test_args.s0, test_args.s1)
    with raises(TypeError):
        law.calculate_voltage(test_args.e0, 100, test_args.s0, test_args.s1)


def test_bad_distance(test_args: Args) -> None:
    sb = Quantity(units.coulomb)
    with raises(errors.UnitsError):
        law.calculate_voltage(test_args.e0, test_args.e1, sb, test_args.s1)
    with raises(TypeError):
        law.calculate_voltage(test_args.e0, test_args.e1, 100, test_args.s1)
    with raises(errors.UnitsError):
        law.calculate_voltage(test_args.e0, test_args.e1, test_args.s0, sb)
    with raises(TypeError):
        law.calculate_voltage(test_args.e0, test_args.e1, test_args.s0, 100)
