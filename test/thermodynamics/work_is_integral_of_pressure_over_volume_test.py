from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    Quantity,
)
from symplyphysics.laws.thermodynamics import work_is_integral_of_pressure_over_volume as work_law

# Description
## At the start of a quasi-static process, the volume of the system was 1 L and the pressure in the system
## was 1 atm. At the end, the volume was 1.5 L and the pressure was 2 atm. Assuming the pressure changed
## linearly with volume during the expansion, the work done by the system amounts to ...

Args = namedtuple("Args", "v0 v1 p0 p1")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    v0 = Quantity(1 * units.liter)
    v1 = Quantity(1.5 * units.liter)
    p0 = Quantity(1 * units.atmosphere)
    p1 = Quantity(2 * units.atmosphere)
    return Args(v0=v0, v1=v1, p0=p0, p1=p1)


def test_law(test_args: Args) -> None:
    result = work_law.calculate_work(test_args.v0, test_args.v1, test_args.p0, test_args.p1)
    assert_equal(result, 76 * units.joule)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_work(vb, test_args.v1, test_args.p0, test_args.p1)
    with raises(TypeError):
        work_law.calculate_work(100, test_args.v1, test_args.p0, test_args.p1)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.v0, vb, test_args.p0, test_args.p1)
    with raises(TypeError):
        work_law.calculate_work(test_args.v0, 100, test_args.p0, test_args.p1)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.v0, test_args.v1, pb, test_args.p1)
    with raises(TypeError):
        work_law.calculate_work(test_args.v0, test_args.v1, 100, test_args.p1)
    with raises(errors.UnitsError):
        work_law.calculate_work(test_args.v0, test_args.v1, test_args.p0, pb)
    with raises(TypeError):
        work_law.calculate_work(test_args.v0, test_args.v1, test_args.p0, 100)
