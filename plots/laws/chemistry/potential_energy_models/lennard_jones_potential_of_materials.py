from dataclasses import dataclass
from sympy import symbols
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend
from symplyphysics.laws.chemistry.potential_energy_models import lennard_jones_potential

print(f"Formula of Lennard-Jones potential:\n{lennard_jones_potential.print_law()}\n")


@dataclass
class PlotData:
    label: str
    dispersion_energy: float  # eV
    particle_size: float  # nm
    line_color: str


# Source: https://www.researchgate.net/figure/Lennard-Jones-LJ-potential-parameters-of-different-materials-considered-in-thepresent_tbl2_319412425
datas_ = [
    PlotData("Ar-Ar", 0.0104, 0.3400, "red"),
    PlotData("Pt-Pt", 0.5200, 0.2475, "blue"),
    PlotData("Ag-Ag", 0.3510, 0.2574, "green"),
    PlotData("Al-Al", 0.4080, 0.2551, "orange"),
]

distance = symbols("distance", positive=True)
law = lennard_jones_potential.law.rhs.subs(
    lennard_jones_potential.distance, distance
)

potentials_plot = plot(
    ylim=(-1.0, 2.0),
    axis_center=(0.0, 0.0),
    line_color="red",
    title="Lennard-Jones potential of different materials",
    xlabel="r, nm",
    ylabel="U, eV",
    backend=MatplotlibBackend,
    show=False,
    annotations=False,
    legend=True,
)

for data_ in datas_:
    applied_law = law.subs({
        lennard_jones_potential.dispersion_energy: data_.dispersion_energy,
        lennard_jones_potential.particle_size: data_.particle_size,
    })
    subplot = plot(
        applied_law, (distance, 0.0, 0.8),
        label=data_.label,
        line_color=data_.line_color,
        show=False,
    )
    potentials_plot.append(subplot[0])

potentials_plot.show()
