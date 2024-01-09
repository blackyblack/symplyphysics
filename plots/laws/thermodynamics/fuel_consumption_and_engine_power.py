from sympy import symbols, Eq
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics import print_expression
from examples.thermodynamics import \
    energy_from_combustion_of_gasoline as combustion_energy_example

# Note: fuel_consumption = volume_of_fuel / distance
fuel_consumption = symbols("fuel_consumption")

# Create efficiency factor values from 15% to 75% in increments of 15%
EFFICIENCY_FACTOR_1 = 0.15
EFFICIENCY_FACTOR_2 = 0.3
EFFICIENCY_FACTOR_3 = 0.45
EFFICIENCY_FACTOR_4 = 0.6
EFFICIENCY_FACTOR_5 = 0.75

fuel_consumption_value = combustion_energy_example.volume_of_gasoline / combustion_energy_example.answer.rhs
fuel_consumption_equation = Eq(fuel_consumption, fuel_consumption_value)
gasoline_consumption_equation = fuel_consumption_equation.subs({
    combustion_energy_example.velocity_of_car: 15,
    # 1 kilometers / hour = (1 / 3.6) meters / second,
    combustion_energy_example.density_of_gasoline: 700,
    # kilogram / (meter ** 3),
    combustion_energy_example.gasoline_specific_heat_combustion: 46 * 1E6
    # 1 megajoules / kilogram = 1E6 joules / kilogram
})

print(f"Formula is:\n{print_expression(gasoline_consumption_equation)}")

gasoline_compustion_value_1 = gasoline_consumption_equation.subs({
    combustion_energy_example.efficiency_factor: EFFICIENCY_FACTOR_1
}).rhs
gasoline_compustion_value_2 = gasoline_consumption_equation.subs({
    combustion_energy_example.efficiency_factor: EFFICIENCY_FACTOR_2
}).rhs
gasoline_compustion_value_3 = gasoline_consumption_equation.subs({
    combustion_energy_example.efficiency_factor: EFFICIENCY_FACTOR_3
}).rhs
gasoline_compustion_value_4 = gasoline_consumption_equation.subs({
    combustion_energy_example.efficiency_factor: EFFICIENCY_FACTOR_4
}).rhs
gasoline_compustion_value_5 = gasoline_consumption_equation.subs({
    combustion_energy_example.efficiency_factor: EFFICIENCY_FACTOR_5
}).rhs

# 1 m^3 / m = 10^(-3) liters / m = (10^(-3) / 10^(-3)) liters / km = 1 liters / km
base_plot = plot(title="Coulomb Law",
                 xlabel="$N, W$",
                 ylabel="$V/S, liters/km$",
                 backend=MatplotlibBackend,
                 legend=True,
                 show=False)

p1 = plot(gasoline_compustion_value_1,
          (combustion_energy_example.power_of_car, 1_000, 100_000),
          label="$\eta_{engine}=" + f"{EFFICIENCY_FACTOR_1}$",
          show=False)
p2 = plot(gasoline_compustion_value_2,
          (combustion_energy_example.power_of_car, 1_000, 100_000),
          label="$\eta_{engine}=" + f"{EFFICIENCY_FACTOR_2}$",
          show=False)
p3 = plot(gasoline_compustion_value_3,
          (combustion_energy_example.power_of_car, 1_000, 100_000),
          label="$\eta_{engine}=" + f"{EFFICIENCY_FACTOR_3}$",
          show=False)
p4 = plot(gasoline_compustion_value_4,
          (combustion_energy_example.power_of_car, 1_000, 100_000),
          label="$\eta_{engine}=" + f"{EFFICIENCY_FACTOR_4}$",
          show=False)
p5 = plot(gasoline_compustion_value_5,
          (combustion_energy_example.power_of_car, 1_000, 100_000),
          label="$\eta_{engine}=" + f"{EFFICIENCY_FACTOR_5}$",
          show=False)

base_plot.append(p1[0])
base_plot.append(p2[0])
base_plot.append(p3[0])
base_plot.append(p4[0])
base_plot.append(p5[0])

base_plot.show()
