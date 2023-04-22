from sympy import solve
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.laws.dynamics import kinetic_energy_from_mass_and_velocity as kinetic_energy

print("Formula is:\n{}".format(kinetic_energy.print()))

solved = solve(kinetic_energy.law,
               kinetic_energy.kinetic_energy_of_body,
               dict=True)[0][kinetic_energy.kinetic_energy_of_body]
result_energy = solved.subs(kinetic_energy.body_mass, 1)

p1 = plot(result_energy, (kinetic_energy.body_velocity, 0, 6),
          ylim=(0, 3),
          axis_center=(0.0, 0.0),
          line_color="red",
          title="Kinetic energy of body (velocity)",
          xlabel="velocity",
          ylabel="energy",
          label="Energy(velocity)",
          legend=True,
          backend=MatplotlibBackend,
          show=False)
p1.show()
