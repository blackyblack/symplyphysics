import pytest
from symplyphysics.laws.relativistic.moving_particle_energy import calculate_particle_energy
from sympy.physics.units import (speed_of_light, 
                                 meter, 
                                 second, 
                                 kg, 
                                 convert_to)
import math

def test_basic():
    e = calculate_particle_energy(1 * kg, 0.5 * speed_of_light)
    check = 1 * kg * convert_to(
        speed_of_light, [meter, second]) ** 2 * \
        math.sqrt(1 / (1 - 0.5 ** 2))
    assert (convert_to(e, [meter, kg, second])) == check
    

def test_exception():
    with pytest.raises(ValueError) as e:
        res = calculate_particle_energy(1 * kg, 1.0 * speed_of_light)
    with pytest.raises(ValueError) as e:
        res = calculate_particle_energy(1 * kg, 1.1 * speed_of_light)


# run from ... simplyphysics/symplyphysics dir
# python -m pytest test/relativistic/moving_particle_energy_test.py        

        
if __name__ == "__main__":
    test_basic()
    test_exception()