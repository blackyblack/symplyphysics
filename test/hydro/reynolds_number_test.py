from collections import namedtuple
from pytest import approx, fixture
from symplyphysics import (
    convert_to,
    Quantity,
    SI,
    dimensionless,
    units
)
from symplyphysics.laws.hydro import reynolds_number

@fixture(name="test_args")
def test_args_fixture():
    rho = Quantity(1000 * units.kilogram / units.meter**3)
    d = Quantity(0.1 * units.meter)
    v = Quantity(1 * units.meter / units.second)
    mu = Quantity(0.000894 * units.pressure * units.time)
    Args = namedtuple("Args", ["rho", "d", "v", "mu"])
    return Args(rho=rho, d=d, v=v, mu=mu)

def test_reynolds_number(test_args):
    result = reynolds_number.calculate_reynolds_number(test_args.d, test_args.rho, test_args.v, test_args.mu)
    assert SI.get_dimension_system().equivalent_dims(result.dimension, dimensionless)
    result_reynolds_number = convert_to(result, dimensionless).evalf(3)
    assert result_reynolds_number == approx(100000, 0.001)