from collections import namedtuple
from pytest import fixture
from sympy import I, Integral, evaluate, exp, log, pi, sin, sqrt
from symplyphysics import Quantity, SymbolNew, clone_as_function, clone_as_symbol, units
from symplyphysics.docs.printer_latex import latex_str

Args = namedtuple(
    "Args",
    [
    "mass",
    "temperature",
    "boltzmann_constant",
    "energy",
    "time",
    "force",
    "speed",
    "speed_of_light",
    "charge",
    "distance",
    "vacuum_permittivity",
    "electric_dipole_moment",
    "intensity",
    ],
)


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass = SymbolNew("m")
    temperature = SymbolNew("T")
    boltzmann_constant = Quantity(display_latex="k_\\text{B}")
    energy = SymbolNew("E")
    time = SymbolNew("t")
    force = SymbolNew("F")
    speed = SymbolNew("v")
    speed_of_light = Quantity(display_symbol="c")
    charge = SymbolNew("q")
    distance = SymbolNew("d")
    vacuum_permittivity = Quantity(display_latex="\\varepsilon_0")
    electric_dipole_moment = SymbolNew("p")
    intensity = SymbolNew("I")
    return Args(
        mass=mass,
        temperature=temperature,
        boltzmann_constant=boltzmann_constant,
        energy=energy,
        time=time,
        force=force,
        speed=speed,
        speed_of_light=speed_of_light,
        charge=charge,
        distance=distance,
        vacuum_permittivity=vacuum_permittivity,
        electric_dipole_moment=electric_dipole_moment,
        intensity=intensity,
    )


def test_add() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a + b
    assert latex_str(expr) == "a + b"


def test_add_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    with evaluate(False):
        expr = a + b + c
    assert latex_str(expr) == "a + b + c"


def test_mul() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a * b
    assert latex_str(expr) == "a b"


def test_div() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a / b
    assert latex_str(expr) == "\\frac{a}{b}"


def test_div_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = a / b + c / d
    assert latex_str(expr) == "\\frac{a}{b} + \\frac{c}{d}"


def test_div_no_braces_common_frac() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = (a + b) / (c + d)
    assert latex_str(expr) == "\\frac{a + b}{c + d}"


def test_sub() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a - b
    assert latex_str(expr) == "a - b"


def test_sub_no_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    with evaluate(False):
        expr = a - b + c
    assert latex_str(expr) == "a - b + c"


def test_mul_with_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    with evaluate(False):
        expr = a * (b + a)
    assert latex_str(expr) == "a \\left(b + a\\right)"
    with evaluate(False):
        expr = (a + b) * a
    assert latex_str(expr) == "\\left(a + b\\right) a"


def test_mul_with_double_braces() -> None:
    a = SymbolNew("a")
    b = SymbolNew("b")
    c = SymbolNew("c")
    d = SymbolNew("d")
    with evaluate(False):
        expr = (a + b) * (c + d)
    assert latex_str(expr) == "\\left(a + b\\right) \\left(c + d\\right)"


def test_mul_with_sqrt_and_quantity(test_args: Args) -> None:
    with evaluate(False):
        expr = sqrt(2 * pi /
            (test_args.mass * test_args.boltzmann_constant * test_args.temperature))
    assert latex_str(expr) == "\\sqrt{\\frac{2 \\pi}{m k_\\text{B} T}}"


def test_mul_minus_one(test_args: Args) -> None:
    with evaluate(False):
        expr = exp(-1 * test_args.energy / (test_args.boltzmann_constant * test_args.temperature))
    assert latex_str(expr) == "\\exp{\\left(- \\frac{E}{k_\\text{B} T} \\right)}"


def test_mul_integral(test_args: Args) -> None:
    time = test_args.time
    time_before = clone_as_symbol(time, subscript="0")
    time_after = clone_as_symbol(time, subscript="1")
    force = clone_as_function(test_args.force)
    with evaluate(False):
        expr = Integral(force(time), (time, time_before, time_after))
    assert latex_str(expr) == "\\int\\limits_{t_{0}}^{t_{1}} F{\\left(t \\right)}\\, dt"


def test_mul_sqrt(test_args: Args) -> None:
    with evaluate(False):
        expr = 1 / sqrt(1 - test_args.speed**2 / test_args.speed_of_light**2)
    assert latex_str(expr) == "\\frac{1}{\\sqrt{1 - \\frac{v^{2}}{c^{2}}}}"


def test_mul_fraction_and_sin(test_args: Args) -> None:
    natural_angular_frequency = SymbolNew("w_0", display_latex="\\omega_0")
    time = test_args.time
    with evaluate(False):
        expr = (test_args.force / (2 * test_args.mass * natural_angular_frequency)) * time * sin(
            natural_angular_frequency * time)
    assert latex_str(expr) == "\\frac{F}{2 m \\omega_{0}} t \\sin{\\left(\\omega_{0} t \\right)}"


def test_fraction_mul_inner_fraction(test_args: Args) -> None:
    first_charge = clone_as_symbol(test_args.charge, subscript="1")
    second_charge = clone_as_symbol(test_args.charge, subscript="2")
    with evaluate(False):
        expr = first_charge * second_charge / (4 * pi *
            test_args.vacuum_permittivity) / test_args.distance**2
    assert latex_str(expr) == "\\frac{q_{1} q_{2}}{4 \\pi \\varepsilon_0} \\frac{1}{d^{2}}"


def test_fraction_mul_fraction(test_args: Args) -> None:
    with evaluate(False):
        expr = 1 / (2 * pi * test_args.vacuum_permittivity) * (test_args.electric_dipole_moment /
            test_args.distance**3)
    assert latex_str(expr) == "\\frac{1}{2 \\pi \\varepsilon_0} \\frac{p}{d^{3}}"


def test_log10(test_args: Args) -> None:
    reference_intensity = Quantity(1e-12 * units.watt / units.meter**2, display_symbol="I_0")
    with evaluate(False):
        expr = log(test_args.intensity / reference_intensity, 10)
    assert latex_str(expr) == "\\log_{10} \\left( \\frac{I}{I_0} \\right)"


def test_minus_imaginary(test_args: Args) -> None:
    with evaluate(False):
        expr = exp(-1 * (I / test_args.boltzmann_constant))
    assert latex_str(expr) == "\\exp{\\left(- \\frac{i}{k_\\text{B}} \\right)}"
