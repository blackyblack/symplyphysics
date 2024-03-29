from collections import namedtuple
from pytest import fixture
from symplyphysics import assert_equal
from symplyphysics.laws.nuclear import most_neutron_energies_scattering_angle_average_cosine

Args = namedtuple("Args", ["A"])


@fixture(name="test_args")
def test_args_fixture() -> Args:
    # carbon mass number is 12
    target_nucleus_mass_number = 12
    return Args(A=target_nucleus_mass_number)


def test_basic_scattering_angle(test_args: Args) -> None:
    result = most_neutron_energies_scattering_angle_average_cosine.calculate_average_scattering_angle_cosine(
        test_args.A)
    assert_equal(result, 0.0555)
