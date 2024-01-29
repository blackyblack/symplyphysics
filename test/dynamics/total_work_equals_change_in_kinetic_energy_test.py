from collections import namedtuple
from pytest import approx, fixture, raises
from symplyphysics import (
    errors,
    units,
    convert_to,
    Quantity,
    SI,
)
from symplyphysics.laws.dynamics import (
    total_work_equals_change_in_kinetic_energy as work_energy_principle,)

# Description
## Several forces act on a particle, so that at first its kinetic energy amounted to 5 J
## and after the interaction to 9 J. The total work done by the forces on the particle is 4 J.


@fixture(name="test_args")
def test_args_fixture():
    K0 = Quantity(5.0 * units.joule)
    K1 = Quantity(9.0 * units.joule)
    Args = namedtuple("Args", "K0 K1")
    return Args(K0=K0, K1=K1)


def test_basic_law(test_args):
    result = work_energy_principle.calculate_total_work(test_args.K0, test_args.K1)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, units.energy)
    result_value = convert_to(result, units.joule).evalf(3)
    assert result_value == approx(4.0, 1e-3)


def test_bad_energies(test_args):
    Kb = Quantity(1.0 * units.coulomb)
    with raises(errors.UnitsError):
        work_energy_principle.calculate_total_work(Kb, test_args.K1)
    with raises(errors.UnitsError):
        work_energy_principle.calculate_total_work(test_args.K0, Kb)
    with raises(TypeError):
        work_energy_principle.calculate_total_work(100, test_args.K1)
    with raises(TypeError):
        work_energy_principle.calculate_total_work(test_args.K0, 100)
