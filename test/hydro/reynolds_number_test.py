from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import Quantity, units
from symplyphysics.laws.hydro import reynolds_number


# Description
## There is pipe with diameter 0.1m, density 1000kg/m3, velocity 1m/s
## and dynamic viscosity 0.000894Pa*s. The reynolds number should be 111856.823

@fixture(name="test_args")
def test_args_fixture():
    rho = Quantity(1000 * units.kilogram / units.meter**3)
    d = Quantity(0.1 * units.meter)
    v = Quantity(1 * units.meter / units.second)
    mu = Quantity(0.000894 * units.pascal * units.second)
    Args = namedtuple("Args", ["d", "rho", "v", "mu"])
    return Args(d=d, rho=rho, v=v, mu=mu)


def test_reynolds_number(test_args):
    result = reynolds_number.calculate_reynolds_number(
        test_args.d, test_args.rho, test_args.v, test_args.mu)
    assert result == approx(111856.823, 0.001)
