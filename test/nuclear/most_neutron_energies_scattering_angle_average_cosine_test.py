from collections import namedtuple
from pytest import approx, fixture
from symplyphysics.laws.nuclear import most_neutron_energies_scattering_angle_average_cosine


@fixture(name="test_args")
def test_args_fixture():
    # carbon mass number is 12
    target_nucleus_mass_number = 12
    Args = namedtuple("Args", ["A"])
    return Args(A=target_nucleus_mass_number)


def test_basic_scattering_angle(test_args):
    result = most_neutron_energies_scattering_angle_average_cosine.calculate_average_scattering_angle_cosine(
        test_args.A)
    assert result == approx(0.0555, 0.01)
