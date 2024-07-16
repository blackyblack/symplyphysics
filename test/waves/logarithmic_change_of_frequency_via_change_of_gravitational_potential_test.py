from collections import namedtuple
from pytest import fixture, raises
from symplyphysics import (
    assert_equal,
    errors,
    units,
    prefixes,
    Quantity,
)
from symplyphysics.laws.waves import (
    logarithmic_change_of_frequency_via_change_of_gravitational_potential as law,
)

Args = namedtuple("Args", "nu1 dphi")



