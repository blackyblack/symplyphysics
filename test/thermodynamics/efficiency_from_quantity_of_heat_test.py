import pytest
from sympy import solve
from sympy.physics.units import temperature
from symplyphysics import units, Quantity, Symbol, validate_input, validate_output, LessEq
from symplyphysics.laws.thermodynamics.efficiency_from_quantity_of_heat import calculate_eta  # Replace "your_module" with the name of your module

# Define epsilon (tolerance) for floating-point number comparison
EPSILON = 1e-6

def test_calculate_eta():
    # Testing the main function calculate_eta

    # Testing with equal temperature values
    result = calculate_eta(Quantity(100, units.joule), Quantity(100, units.joule))
    assert result.evalf() == 0.0

    # Testing with final temperature greater than initial temperature
    result = calculate_eta(Quantity(200, units.joule), Quantity(100, units.joule))
    assert result.evalf() == 0.5

    # Testing with initial temperature equal to 0
    result = calculate_eta(Quantity(0, units.joule), Quantity(100, units.joule))
    assert result.evalf() == 1.0


# Additional tests for validate_input and validate_output
def test_validate_input():
    # Testing input data validation

    # Correct temperature values
    assert validate_input(temperature_start_=Quantity(100, units.joule), temperature_end_=Quantity(50, units.joule))

    # Error: final temperature > initial temperature
    with pytest.raises(ValueError):
        validate_input(temperature_start_=Quantity(50, units.joule), temperature_end_=Quantity(100, units.joule))


def test_validate_output():
    # Testing output data validation

    # Correct efficiency value
    assert validate_output(Quantity(0.5, units.joule / units.joule))

    # Error: efficiency > 1
    with pytest.raises(ValueError):
        validate_output(Quantity(1.5, units.joule / units.joule))

