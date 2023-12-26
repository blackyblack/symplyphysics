import pytest
from sympy import sin, rad
from your_optics_module import calculate_refraction_angle

# Replace 'your_optics_module' with the actual name of your code file

def test_calculate_refraction_angle():
    # Test with known values
    theta_i_value = rad(30)  # incident angle of 30 degrees
    n1_value = 1.5  # refractive index of the first medium
    n2_value = 1.2  # refractive index of the second medium

    # Calculate the expected result
    expected_result = sin(theta_i_value) * n1_value / n2_value

    # Call the function to calculate the refraction angle
    result = calculate_refraction_angle(theta_i_value, n1_value, n2_value)

    # Check that the result is close to the expected value
    assert abs(result - expected_result) < 1e-6

    # Test with different values
    theta_i_value = rad(45)
    n1_value = 1.4
    n2_value = 1.6

    expected_result = sin(theta_i_value) * n1_value / n2_value
    result = calculate_refraction_angle(theta_i_value, n1_value, n2_value)

    assert abs(result - expected_result) < 1e-6

if __name__ == "__main__":
    pytest.main()
