from collections import namedtuple
from pytest import fixture
from sympy import I, Integral, evaluate, exp, log, pi, sin, sqrt, Matrix, ImmutableMatrix
from symplyphysics import (
    Quantity,
    Symbol,
    clone_as_function,
    clone_as_symbol,
    units,
    IndexedSymbol,
    IndexedSum,
    IndexedProduct,
    global_index,
)
from symplyphysics.docs.printer_latex import latex_str
from symplyphysics.core.operations import symbolic

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
    mass = Symbol("m")
    temperature = Symbol("T")
    boltzmann_constant = Quantity(display_latex="k_\\text{B}")
    energy = Symbol("E")
    time = Symbol("t")
    force = Symbol("F")
    speed = Symbol("v")
    speed_of_light = Quantity(display_symbol="c")
    charge = Symbol("q")
    distance = Symbol("d")
    vacuum_permittivity = Quantity(display_latex="\\varepsilon_0")
    electric_dipole_moment = Symbol("p")
    intensity = Symbol("I")
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
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a + b
    assert latex_str(expr) == "a + b"


def test_add_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    with evaluate(False):
        expr = a + b + c
    assert latex_str(expr) == "a + b + c"


def test_mul() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a * b
    assert latex_str(expr) == "a b"


def test_div() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a / b
    assert latex_str(expr) == "\\frac{a}{b}"


def test_div_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = a / b + c / d
    assert latex_str(expr) == "\\frac{a}{b} + \\frac{c}{d}"


def test_div_no_braces_common_frac() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = (a + b) / (c + d)
    assert latex_str(expr) == "\\frac{a + b}{c + d}"


def test_sub() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a - b
    assert latex_str(expr) == "a - b"


def test_sub_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    with evaluate(False):
        expr = a - b + c
    assert latex_str(expr) == "a - b + c"


def test_mul_with_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a * (b + a)
    assert latex_str(expr) == "a \\left(b + a\\right)"
    with evaluate(False):
        expr = (a + b) * a
    assert latex_str(expr) == "\\left(a + b\\right) a"


def test_mul_with_double_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
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
    natural_angular_frequency = Symbol("w_0", display_latex="\\omega_0")
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


def test_log_squared(test_args: Args) -> None:
    with evaluate(False):
        expr = log(test_args.boltzmann_constant)**2
    assert latex_str(expr) == "\\log \\left( k_\\text{B} \\right)^{2}"


def test_mutable_matrix() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")

    # 2-by-2 matrix
    with evaluate(False):
        expr = Matrix([[a, b], [c, d]])
    assert latex_str(expr) == "\\begin{pmatrix} a & b \\\\ c & d \\end{pmatrix}"

    # 4-row vector
    with evaluate(False):
        expr = Matrix([a, b, c, d])
    assert latex_str(expr) == "\\begin{pmatrix} a \\\\ b \\\\ c \\\\ d \\end{pmatrix}"

    # 4-col vector
    with evaluate(False):
        expr = Matrix([a, b, c, d]).T
    assert latex_str(expr) == "\\begin{pmatrix} a & b & c & d \\end{pmatrix}"


def test_immutable_matrix() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")

    # 2-by-2 matrix
    with evaluate(False):
        expr = ImmutableMatrix([[a, b], [c, d]])
    assert latex_str(expr) == "\\begin{pmatrix} a & b \\\\ c & d \\end{pmatrix}"

    # 4-row vector
    with evaluate(False):
        expr = ImmutableMatrix([a, b, c, d])
    assert latex_str(expr) == "\\begin{pmatrix} a \\\\ b \\\\ c \\\\ d \\end{pmatrix}"

    # 4-col vector
    with evaluate(False):
        expr = ImmutableMatrix([a, b, c, d]).T
    assert latex_str(expr) == "\\begin{pmatrix} a & b & c & d \\end{pmatrix}"


def test_indexed_symbol() -> None:
    f = IndexedSymbol("F")
    expr = f[global_index]
    assert latex_str(expr) == "{F}_{i}"


def test_indexed_sum() -> None:
    f = IndexedSymbol("F")
    expr = IndexedSum(f[global_index], global_index)
    assert latex_str(expr) == "\\sum_i {F}_{i}"

    with evaluate(False):
        expr = IndexedSum(f[global_index]**2, global_index)
    assert latex_str(expr) == "\\sum_i {F}_{i}^{2}"


def test_indexed_product() -> None:
    f = IndexedSymbol("F")

    expr = IndexedProduct(f[global_index], global_index)
    assert latex_str(expr) == "\\prod_i {F}_{i}"

    with evaluate(False):
        expr = IndexedProduct(f[global_index]**2, global_index)
    assert latex_str(expr) == "\\prod_i {F}_{i}^{2}"


def test_average() -> None:
    a = Symbol("a")
    expr = symbolic.Average(a)
    assert latex_str(expr) == "\\langle a \\rangle"


def test_finite_difference() -> None:
    a = Symbol("a")

    expr = symbolic.FiniteDifference(a)
    assert latex_str(expr) == "\\Delta a"

    expr = symbolic.FiniteDifference(a, wrap_latex=True)
    assert latex_str(expr) == "\\Delta \\left( a \\right)"


def test_exact_differential() -> None:
    a = Symbol("a")

    expr = symbolic.ExactDifferential(a)
    assert latex_str(expr) == "d a"

    expr = symbolic.ExactDifferential(a, wrap_latex=True)
    assert latex_str(expr) == "d \\left( a \\right)"
