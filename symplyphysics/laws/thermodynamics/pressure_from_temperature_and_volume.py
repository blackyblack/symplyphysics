from sympy import symbols, Eq, pretty, solve
from sympy.physics.units import molar_gas_constant

# Description
## Ideal gas law: P * V = n * R * T
## Where:
## P is pressure,
## V is volume,
## n is number of moles,
## R is ideal gas constant,
## T is temperature

pressure, volume, mole_count, temperature, R_constant = symbols('pressure volume mole_count temperature R_constant')
law = Eq(pressure, mole_count * temperature * R_constant / volume)

def print():
    return pretty(law, use_unicode=False)

def calculate_pressure(volume_, temperature_, mole_count_):
    return solve(law
      .subs(R_constant, molar_gas_constant)
      .subs(volume, volume_)
      .subs(temperature, temperature_)
      .subs(mole_count, mole_count_))[0].evalf()