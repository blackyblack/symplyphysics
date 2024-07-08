#!/usr/bin/env python3

from dataclasses import dataclass
from pytest import approx
from sympy import symbols, solve, pi, sin
from sympy.plotting import plot
from sympy.plotting.plot import MatplotlibBackend, Plot
from symplyphysics.laws.optics import refraction_angle_from_environments as refraction_law

# Description
## Plot the dependency of refraction angle on incidence angle for different environments

# all indices were measured at the wavelength of 587.6 nm
AIR_REFRACTIVE_INDEX = 1.0003  # Source: https://refractiveindex.info/?shelf=other&book=air&page=Ciddor
WATER_REFRACTIVE_INDEX = 1.3325  # Source: https://refractiveindex.info/?shelf=main&book=H2O&page=Hale
BENZENE_REFRACTIVE_INDEX = 1.4957  # Source: https://refractiveindex.info/?shelf=organic&book=benzene&page=Moutzouris
DIAMOND_REFRACTIVE_INDEX = 2.4168  # Source: https://refractiveindex.info/?shelf=3d&book=crystals&page=diamond

incident_index, refracted_index = symbols("incident_index refracted_index", positive=True)
angle_of_incidence, angle_of_refraction = symbols("angle_of_incidence angle_of_refraction",
    real=True)

law = refraction_law.law.subs({
    refraction_law.incidence_refractive_index: incident_index,
    refraction_law.resulting_refractive_index: refracted_index,
    refraction_law.incidence_angle: angle_of_incidence,
    refraction_law.refraction_angle: angle_of_refraction,
})

# first solution is (pi - angle), ignore it
solved_angle_of_refraction = solve(law, angle_of_refraction)[1]


def maximum_angle_of_incidence(incident_index_: float, refracted_index_: float) -> float:
    # sin(angle_of_incidence) <= refracted_index / incident_index
    eqn = law.subs({
        sin(angle_of_refraction): 1,
        incident_index: incident_index_,
        refracted_index: refracted_index_,
    })
    solved = solve(eqn, angle_of_incidence)

    # if the equation has no solution, it means that the maximum angle is pi/2
    return solved[0] if solved else pi / 2


assert maximum_angle_of_incidence(AIR_REFRACTIVE_INDEX,
    WATER_REFRACTIVE_INDEX) == approx(pi / 2, 1e-5)
assert maximum_angle_of_incidence(WATER_REFRACTIVE_INDEX,
    AIR_REFRACTIVE_INDEX) == approx(0.84911, 1e-5)


@dataclass
class SubplotData:
    incident_index: float
    refracted_index: float
    label: str
    color: str


def make_subplot(data_: SubplotData) -> Plot:
    maximum_angle = maximum_angle_of_incidence(data_.incident_index, data_.refracted_index)
    angle_of_refraction_value = solved_angle_of_refraction.subs({
        incident_index: data_.incident_index,
        refracted_index: data_.refracted_index,
        # convert back, from degrees to radians
        angle_of_incidence: angle_of_incidence * pi / 180
    })
    return plot(
        # convert Y result to degrees
        angle_of_refraction_value * 180 / pi,
        # convert X values to degrees (requires to convert back to degrees in the formula)
        (angle_of_incidence, 0, maximum_angle * 180 / pi),
        label=data_.label,
        line_color=data.color,
        show=False,
    )


datas = [
    SubplotData(
    incident_index=DIAMOND_REFRACTIVE_INDEX,
    refracted_index=AIR_REFRACTIVE_INDEX,
    label="diamond to air",
    color="black",
    ),
    SubplotData(
    incident_index=BENZENE_REFRACTIVE_INDEX,
    refracted_index=WATER_REFRACTIVE_INDEX,
    label="benzene to water",
    color="blue",
    ),
    SubplotData(
    incident_index=WATER_REFRACTIVE_INDEX,
    refracted_index=BENZENE_REFRACTIVE_INDEX,
    label="water to benzene",
    color="green",
    ),
    SubplotData(
    incident_index=AIR_REFRACTIVE_INDEX,
    refracted_index=WATER_REFRACTIVE_INDEX,
    label="air to water",
    color="red",
    ),
]

p = plot(
    title="Angle of refraction against angle of incidence",
    xlabel="Angle of incidence, deg",
    ylabel="Angle of refraction, deg",
    backend=MatplotlibBackend,
    legend=True,
    show=False,
)
for data in datas:
    sub_p = make_subplot(data)
    p.append(sub_p[0])
p.show()
