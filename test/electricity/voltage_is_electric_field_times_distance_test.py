from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import assert_equal, errors, units, Quantity
from symplyphysics.laws.electricity import voltage_is_electric_field_times_distance as law

Args = namedtuple("Args", "e d")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    e = Quantity(0.5 * units.volt / units.meter)
    d = Quantity(0.1 * units.meter)
    return Args(e=e, d=d)


def test_law(test_args: Args) -> None:
    result = law.calculate_voltage(test_args.e, test_args.d)
    assert_equal(result, 0.05 * units.volt)


def test_bad_electric_field(test_args: Args) -> None:
    eb = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_voltage(eb, test_args.d)
    with raises(TypeError):
        law.calculate_voltage(100, test_args.d)


def test_bad_distance(test_args: Args) -> None:
    db = Quantity(units.candela)
    with raises(errors.UnitsError):
        law.calculate_voltage(test_args.e, db)
    with raises(TypeError):
        law.calculate_voltage(test_args.e, 100)
