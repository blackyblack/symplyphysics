from collections import namedtuple
from pytest import fixture
from sympy import I, Integral, evaluate, exp, log, pi, sin, sqrt, Matrix, ImmutableMatrix
from symplyphysics import (
    Quantity,
    Symbol,
    clone_as_symbol,
    clone_as_function,
    units,
    SymbolIndexedNew,
    SumIndexed,
    ProductIndexed,
    global_index,
)
from symplyphysics.docs.printer_code import code_str

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
    "intensity",
    ],
)


@fixture(name="test_args")
def test_args_fixture() -> Args:
    mass = Symbol("m")
    temperature = Symbol("T")
    boltzmann_constant = Quantity(display_symbol="k_B")
    energy = Symbol("E")
    time = Symbol("t")
    force = Symbol("F")
    speed = Symbol("v")
    speed_of_light = Quantity(display_symbol="c")
    charge = Symbol("q")
    distance = Symbol("d")
    vacuum_permittivity = Quantity(display_symbol="epsilon_0")
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
        intensity=intensity,
    )


def test_add() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a + b
    assert code_str(expr) == "a + b"


def test_add_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    with evaluate(False):
        expr = a + b + c
    assert code_str(expr) == "a + b + c"


def test_mul() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a * b
    assert code_str(expr) == "a * b"


def test_div() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a / b
    assert code_str(expr) == "a / b"


def test_div_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = a / b + c / d
    assert code_str(expr) == "a / b + c / d"


def test_div_braces_common_frac() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = (a + b) / (c + d)
    assert code_str(expr) == "(a + b) / (c + d)"


def test_sub() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a - b
    assert code_str(expr) == "a - b"


def test_sub_no_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    with evaluate(False):
        expr = a - b + c
    assert code_str(expr) == "a - b + c"


def test_mul_with_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    with evaluate(False):
        expr = a * (b + a)
    assert code_str(expr) == "a * (b + a)"
    with evaluate(False):
        expr = (a + b) * a
    assert code_str(expr) == "(a + b) * a"


def test_mul_with_double_braces() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")
    with evaluate(False):
        expr = (a + b) * (c + d)
    assert code_str(expr) == "(a + b) * (c + d)"


def test_mul_with_sqrt_and_quantity(test_args: Args) -> None:
    with evaluate(False):
        expr = sqrt(2 * pi /
            (test_args.mass * test_args.boltzmann_constant * test_args.temperature))
    assert code_str(expr) == "sqrt(2 * pi / (m * k_B * T))"


def test_mul_minus_one(test_args: Args) -> None:
    with evaluate(False):
        expr = exp(-1 * test_args.energy / (test_args.boltzmann_constant * test_args.temperature))
    assert code_str(expr) == "exp(-E / (k_B * T))"


def test_mul_integral(test_args: Args) -> None:
    time = test_args.time
    time_before = clone_as_symbol(time, subscript="0")
    time_after = clone_as_symbol(time, subscript="1")
    force = clone_as_function(test_args.force)
    with evaluate(False):
        expr = Integral(force(time), (time, time_before, time_after))
    assert code_str(expr) == "Integral(F(t), (t, t_0, t_1))"


def test_mul_sqrt(test_args: Args) -> None:
    with evaluate(False):
        expr = 1 / sqrt(1 - test_args.speed**2 / test_args.speed_of_light**2)
    assert code_str(expr) == "1 / sqrt(1 - v^2 / c^2)"


def test_mul_fraction_and_sin(test_args: Args) -> None:
    natural_angular_frequency = Symbol("w_0", display_latex="\\omega_0")
    time = test_args.time
    with evaluate(False):
        expr = (test_args.force / (2 * test_args.mass * natural_angular_frequency)) * time * sin(
            natural_angular_frequency * time)
    assert code_str(expr) == "F / (2 * m * w_0) * t * sin(w_0 * t)"


def test_fraction_mul_inner_fraction(test_args: Args) -> None:
    first_charge = clone_as_symbol(test_args.charge, subscript="1")
    second_charge = clone_as_symbol(test_args.charge, subscript="2")
    with evaluate(False):
        expr = first_charge * second_charge / (4 * pi *
            test_args.vacuum_permittivity) / test_args.distance**2
    assert code_str(expr) == "q_1 * q_2 / (4 * pi * epsilon_0) / d^2"


def test_fraction_mul_inner_fraction_with_braces(test_args: Args) -> None:
    first_charge = clone_as_symbol(test_args.charge, subscript="1")
    second_charge = clone_as_symbol(test_args.charge, subscript="2")
    with evaluate(False):
        expr = 1 / (4 * pi * test_args.vacuum_permittivity) * (first_charge * second_charge /
            test_args.distance**2)
    assert code_str(expr) == "q_1 * q_2 / d^2 / (4 * pi * epsilon_0)"


def test_log10(test_args: Args) -> None:
    reference_intensity = Quantity(1e-12 * units.watt / units.meter**2, display_symbol="I_0")
    with evaluate(False):
        expr = log(test_args.intensity / reference_intensity, 10)
    assert code_str(expr) == "log(I / I_0, 10)"


def test_minus_imaginary(test_args: Args) -> None:
    with evaluate(False):
        expr = exp(-1 * (I / test_args.boltzmann_constant))
    assert code_str(expr) == "exp(-I / k_B)"


def test_log_squared(test_args: Args) -> None:
    with evaluate(False):
        expr = log(test_args.boltzmann_constant)**2
    assert code_str(expr) == "log(k_B)^2"


def test_mutable_matrix() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")

    # 2-by-2 matrix
    with evaluate(False):
        expr = Matrix([[a, b], [c, d]])
    assert code_str(expr) == "[[a, b], [c, d]]"

    # 4-row vector
    with evaluate(False):
        expr = Matrix([a, b, c, d])
    assert code_str(expr) == "[a, b, c, d]"

    # 4-column vector
    with evaluate(False):
        expr = Matrix([a, b, c, d]).T
    assert code_str(expr) == "[a, b, c, d].T"


def test_immutable_matrix() -> None:
    a = Symbol("a")
    b = Symbol("b")
    c = Symbol("c")
    d = Symbol("d")

    # 2-by-2 matrix
    with evaluate(False):
        expr = ImmutableMatrix([[a, b], [c, d]])
    assert code_str(expr) == "[[a, b], [c, d]]"

    # 4-row vector
    with evaluate(False):
        expr = ImmutableMatrix([a, b, c, d])
    assert code_str(expr) == "[a, b, c, d]"

    # 4-column vector
    with evaluate(False):
        expr = ImmutableMatrix([a, b, c, d]).T
    assert code_str(expr) == "[a, b, c, d].T"


def test_indexed_symbol() -> None:
    f = SymbolIndexedNew("F")
    expr = f[global_index]
    assert code_str(expr) == "F[i]"


def test_indexed_sum() -> None:
    f = SymbolIndexedNew("F")

    expr = SumIndexed(f[global_index], global_index)
    assert code_str(expr) == "Sum(F[i], i)"

    with evaluate(False):
        expr = SumIndexed(f[global_index]**2, global_index)
    assert code_str(expr) == "Sum(F[i]^2, i)"


def test_indexed_product() -> None:
    f = SymbolIndexedNew("F")

    expr = ProductIndexed(f[global_index], global_index)
    assert code_str(expr) == "Product(F[i], i)"

    with evaluate(False):
        expr = ProductIndexed(f[global_index]**2, global_index)
    assert code_str(expr) == "Product(F[i]^2, i)"
