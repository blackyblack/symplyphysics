from sympy import symbols, Eq
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from examples.thermodynamics import energy_from_combustion_of_gasoline as combustion_energy_example

# Note: fuel_consumption = volume_of_fuel / distance
fuel_consumption = symbols("fuel_consumption")

# Create efficiency factor values from 15% to 75% in increments of 15%
efficiency_factors = (round(0.15 * (i + 1), 2) for i in range(5))

fuel_consumption_value = (combustion_energy_example.answer.rhs / combustion_energy_example.volume_of_gasoline) ** (-1)
fuel_consumption_equation = Eq(fuel_consumption, fuel_consumption_value)
gasoline_consumption_equation = fuel_consumption_equation.subs({
    combustion_energy_example.velocity_of_car: 15,      # 1 kilometers / hour = (1 / 3.6) meters / second,
    combustion_energy_example.density_of_gasoline: 700,     # kilogram / (meter ** 3),
    combustion_energy_example.gasoline_specific_heat_combustion: 46 * 1E6       # 1 megajoules / kilogram = 1E6 joules / kilogram
})

print(f"Formula is:\n{print_expression(gasoline_consumption_equation)}")

# 1 m^3 / m = 10^(-3) liters / m = (10^(-3) / 10^(-3)) liters / km = 1 liters / km
base_plot = plot(title="Coulomb Law",
                 xlabel="$N, W$",
                 ylabel="$V/S, l/km$",
                 backend=MatplotlibBackend,
                 legend=True,
                 show=False)

# Create plots for every efficient factor in sequence and add plot to base plot
for efficiency_factor in efficiency_factors:
    gasoline_compustion_value = gasoline_consumption_equation.subs({
        combustion_energy_example.efficiency_factor: efficiency_factor
    }).rhs
    subplot = plot(gasoline_compustion_value, (combustion_energy_example.power_of_car, 1_000, 100_000),
                   label="$\eta_{engine}=" + f"{efficiency_factor}$",
                   show=False)
    base_plot.append(subplot[0])

base_plot.show()
