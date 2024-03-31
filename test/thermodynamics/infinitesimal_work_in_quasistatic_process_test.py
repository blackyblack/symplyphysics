from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.thermodynamics import infinitesimal_work_in_quasistatic_process as work_law

# Description
## In an infinitesimal quasi-static process the pressure of the system was 1 Pa and the infinitesimal
## change of volume was 1 mm**3. The infinitesimal work done by the system was 1 nJ.

Args = namedtuple("Args", "p dv")


@fixture(name="test_args")
def test_args_fixture() -> Args:
    p = Quantity(1 * units.pascal)
    dv = Quantity(1 * units.millimeter**3)
    return Args(p=p, dv=dv)


def test_law(test_args: Args) -> None:
    result = work_law.calculate_infinitesimal_work_done(test_args.p, test_args.dv)
    assert_equal(result, 1 * prefixes.nano * units.joule)


def test_bad_pressure(test_args: Args) -> None:
    pb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_infinitesimal_work_done(pb, test_args.dv)
    with raises(TypeError):
        work_law.calculate_infinitesimal_work_done(100, test_args.dv)


def test_bad_volume(test_args: Args) -> None:
    vb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        work_law.calculate_infinitesimal_work_done(test_args.p, vb)
    with raises(TypeError):
        work_law.calculate_infinitesimal_work_done(test_args.p, 100)
