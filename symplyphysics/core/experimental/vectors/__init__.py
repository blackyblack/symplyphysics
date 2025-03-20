from __future__ import annotations

from typing import Any, Optional

from sympy import Atom, Basic, Expr, S, sympify, ask, Q, simplify
from sympy.core.parameters import global_parameters
from sympy.physics.units import Dimension
from sympy.physics.units.systems.si import dimsys_SI
from sympy.printing.printer import Printer

from symplyphysics.core.dimensions import collect_expression_and_dimension
from symplyphysics.core.errors import UnitsError
from symplyphysics.core.symbols.id_generator import next_id
from symplyphysics.core.symbols.symbols import DimensionSymbol
from symplyphysics.docs.miscellaneous import needs_mul_brackets


class VectorExpr(Basic):  # type: ignore[misc]
    """
    Base class for all vector expressions.
    """

    def doit(self, **_hints: Any) -> VectorExpr:
        return self

    def __add__(self, other: VectorExpr) -> VectorExpr:
        return VectorAdd(self, other).doit()

    def __sub__(self, other: VectorExpr) -> VectorExpr:
        # Refer to the consequence in VectorAdd
        return VectorAdd(self, -other).doit()

    def __mul__(self, other: Any) -> VectorExpr:
        return VectorScale(self, other).doit()

    def __neg__(self) -> VectorExpr:
        # Refer to consequence #4 in VectorScale
        return VectorScale(self, -1).doit()

    def __pos__(self) -> VectorExpr:
        return self.doit()

    def __truediv__(self, other: Any) -> VectorExpr:
        # NOTE: probably need to check if `other` is not `0` to return a special "NaN" vector

        return VectorScale(self, 1 / other).doit()

    @property
    def is_zero(self) -> bool:
        # check with `isinstance` in case the user instantiates their own zero vector.
        return isinstance(self, _VectorZero)

    def subs(self, *args: Any, **kwargs: Any) -> VectorExpr:
        return super().subs(*args, **kwargs)  # type: ignore[no-any-return]

    def as_symbol_combination(self) -> tuple[tuple[VectorSymbol, Expr], ...]:
        """
        Express `self` as a linear combination of `VectorSymbol`. Each term is represented by a
        `(VectorSymbol, Expr)` tuple, and each symbol must appear in the combination only once.

        Examples:
        =========

        >>> from symplyphysics.core.experimental.coordinate_systems import CartesianCoordinateSystem
        >>> c = CartesianCoordinateSystem()
        >>> VectorAdd(c.i, VectorScale(c.j / S(2), S(4)), ZERO, -c.i * 2).as_symbol_combination()
        ((i, -1), (j, 2))
        """

        raise NotImplementedError(f"Implement this method in {type(self).__name__}.")


class _VectorZero(VectorExpr):
    """
    Class expressing the notion of a zero vector. This class isn't intended to be instantiated
    except for the constant `ZERO` since under the `definition of vector spaces
    <https://en.wikipedia.org/wiki/Vector_space#Definition_and_basic_properties>` there can only
    be one zero.

    Note that the zero vector can be considered having arbitrary (physical) dimension as it can be
    added to all vectors.
    """

    def _sympystr(self, _p: Printer) -> str:
        return "0"

    def as_symbol_combination(self) -> tuple[tuple[VectorSymbol, Expr], ...]:
        return ()


ZERO = _VectorZero()


# NOTE: Instead of marking the `norm` of the `VectorSymbol` on the symbol itself, another
#       possibility is to create a separate class for `UnitVector`s. This way, when we add support
#       for vector components, we wouldn't need to check for the fact that the norm calculated
#       using the supplied components equals the norm given at the instantiation of the symbol.
#       But perhaps this simply gives us additional information about the vector and there's no
#       need to worry.
class VectorSymbol(DimensionSymbol, VectorExpr, Atom):  # type: ignore[misc]
    """
    Class representing a symbolic vector.

    The ``display_symbol``, ``display_latex``, and ``dimension`` parameters are used to instantiate
    the ``DimensionSymbol`` class.

    The ``norm`` parameter represents the **norm**, in other words **length** or **magnitude**, of
    the vector. It can be a number, an expression, a quantity, or a symbol, but its dimension must
    match the ``dimension`` of the vector symbol. Note that only ``dimensionless`` vectors with
    ``norm=1`` can be considered as **unit vectors**. As a counterexample, a force vector of
    magnitude `1 N` is not a unit vector since its norm contains a dimensionful quantity `N`.
    """

    _norm: Optional[Expr]
    _id: int

    is_symbol = True

    def __new__(
        cls,
        display_symbol: Optional[str] = None,
        dimension: Dimension = Dimension(1),
        *,
        norm: Optional[Any] = None,  # pylint: disable=redefined-outer-name
        display_latex: Optional[str] = None,
    ) -> VectorExpr:
        if norm is not None:
            norm = simplify(norm)

            if not isinstance(norm, Expr):
                raise TypeError(f"Norm {norm} must be an Expr, got {type(norm).__name__}.")

            if not ask(Q.nonnegative(norm)):  # pylint: disable=too-many-function-args
                raise ValueError(f"Norm must be non-negative, got {norm}.")

            if norm == 0:
                return ZERO

        return VectorExpr.__new__(cls)

    def __init__(
        self,
        display_symbol: Optional[str] = None,
        dimension: Dimension = Dimension(1),
        *,
        norm: Optional[Any] = None,  # pylint: disable=redefined-outer-name
        display_latex: Optional[str] = None,
    ) -> None:
        self._id = next_id("VEC")

        if not display_symbol:
            display_symbol = f"VEC{self._id}"
            if not display_latex:
                display_latex = f"\\mathbf{{v}}_{{{self._id}}}"
        elif not display_latex:
            display_latex = f"\\mathbf{{{display_symbol}}}"

        DimensionSymbol.__init__(
            self,
            display_name=display_symbol,
            dimension=dimension,
            display_latex=display_latex,
        )

        if norm is not None:
            norm_dim = collect_expression_and_dimension(norm)[1]
            if not dimsys_SI.equivalent_dims(norm_dim, self.dimension):
                raise UnitsError(f"The norm must be {self.dimension}, got {norm_dim}.")

        self._norm = norm

    @property
    def norm(self) -> Optional[Expr]:
        return self._norm

    def _hashable_content(self) -> tuple[Any, ...]:
        return (self._id,)

    def as_symbol_combination(self) -> tuple[tuple[VectorSymbol, Expr], ...]:
        return ((self, S.One),)


class VectorNorm(Expr):  # type: ignore[misc]
    """
    Class representing the Euclidean norm (see link 1, *Euclidean norm*) of a vector expression.

    The vector argument is stored in position `0` of `self.args`.

    The vector norm has the following **properties**:

    1. **Subadditivity**: for all vectors `a` and `b`, `norm(a + b) <= norm(a) + norm(b)`.

    2. **Absolute homogeneity**: for all scalars `k` and vectors `a`, `norm(k * a) = abs(k) * norm(a)`.

    3. **Positive definiteness**: for all vectors `a`, `norm(a) = 0` if and only if `a = 0`.

    **Links:**

    1. `Wikipedia <https://en.wikipedia.org/wiki/Norm_(mathematics)>`__.
    """

    is_nonnegative = True

    @property
    def argument(self) -> VectorExpr:
        return self.args[0]  # type: ignore[no-any-return]

    # NOTE: Add __new__ that would dispatch the code execution depending on the value of `vector`
    # For now, this is handled by the function `norm` below.
    def __init__(self, vector: VectorExpr) -> None:
        self._args = (vector,)

    def doit(self, **hints: Any) -> Expr:
        vector = self.argument

        if hints.get("deep", True):
            vector = vector.doit()

        # Refer to property #3
        if vector.is_zero:
            return S.Zero

        if isinstance(vector, VectorSymbol) and vector.norm is not None:
            return vector.norm

        # Refer to property #2
        if isinstance(vector, VectorScale):
            return VectorNorm(vector.args[0]) * abs(vector.args[1])

        # TODO: add support for the following relation:
        # for all vectors `a, b` and scalars `k`, `norm(a * k + b * k) = norm(a + b) * abs(k)`

        return self

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return f"norm({p.doprint(self.argument)})"


def norm(vector: VectorExpr) -> Expr:
    return VectorNorm(vector).doit()


class VectorScale(VectorExpr):
    """
    Class representing the notion of scalar multiplication as a property of vectors.

    This operation has the following **properties**:

    1. Field and vector multiplications are compatible: for all scalars `k, l` and vectors `a`,
       `k * (s * a) = (k * s) * a`.

    2. Identity element exists: for all vectors `a`, `1 * a = 1`.

    3. Distributivity of scalar multiplication w.r.t. vector addition: for all scalars `k` and
       vectors `a, b`, `k * (a + b) = k * a + k * b`.

    The last property, distributivity of scalar multiplication w.r.t. field addition, is not
    represented within the functionality of `VectorScale`.

    As a **consequence** of the properties of the vector field, one has

    1. For all vectors `a`, `0 * k = 0`.

    2. For all scalars `k`, `k * 0 = 0` where `0` is the zero vector.

    3. For all scalars `k` and vectors `a`, `k * a = 0` implies `k = 0` or `a = 0`.

    4. For all vectors `a`, `(-1) * a = -a` where `-a` is the additive inverse of `a`.

    **Links:**

    1. `Wikipedia <https://en.wikipedia.org/wiki/Vector_space#Definition_and_basic_properties>`__.
    """

    @property
    def vector(self) -> VectorExpr:
        return self.args[0]  # type: ignore[no-any-return]

    @property
    def scale(self) -> Expr:
        return self.args[1]

    def __new__(cls, vector: VectorExpr, scale: Any, **kwargs: Any) -> VectorScale:
        # TODO: Add dispatch depending on the value of `vector` and `scale`?

        return super().__new__(cls)  # type: ignore[no-any-return]

    def __init__(self, vector: VectorExpr, scale: Any, **kwargs: Any) -> None:
        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            scale = sympify(scale, strict=True)
            if not isinstance(scale, Expr):
                raise TypeError(f"Scale {scale} must be an Expr, got {type(scale).__name__}.")

        self._args = (vector, scale)

    def doit(self, **_hints: Any) -> VectorExpr:
        combination = self.as_symbol_combination()

        match combination:
            case ():
                # Refer to consequence #2
                return ZERO
            case ((v, s),):
                # Refer to property #2
                if s == 1:
                    return v

                # Refer to consequence #1
                if s == 0:
                    return ZERO

                return VectorScale(v, s, evaluate=False)
            case _:
                # Refer to property #3
                return VectorAdd(*(v * s for v, s in combination))

    def _sympystr(self, p: Printer) -> str:
        vector, value = self.args

        if needs_mul_brackets(value, last=True):
            return f"{p.doprint(vector)}*({p.doprint(value)})"

        return f"{p.doprint(vector)}*{p.doprint(value)}"

    def as_symbol_combination(self) -> tuple[tuple[VectorSymbol, Expr], ...]:
        scale = self.scale
        vector = self.vector

        return tuple((v, s * scale) for v, s in vector.as_symbol_combination())


class VectorAdd(VectorExpr):
    """
    Class representing the notion of vector *addition* as a property of vectors.

    Note that the addends must have the same (physical) dimension to be added together.

    This operation has the following properties:

    1. **Associativity**: for all vectors `a, b, c`, `a + (b + c) = (a + b) + c`.

    2. **Commutativity**: for all vectors `a, b`, `a + b = b + a`.

    3. Existence of **identity vector** `0`: for all vectors `a`, `a + 0 = a`.

    4. Existence of **inverse vector**: for all vectors `a`, there exists a vector `-a` s.t.
    `a + (-a) = 0`.

    As a **consequence**, the *subtraction* of two vectors can be defined as such: for all vectors
    `a, b`, `a - b = a + (-b)`.

    **Links:**

    1. `Wikipedia <https://en.wikipedia.org/wiki/Vector_space#Definition_and_basic_properties>`__.
    """

    @property
    def addends(self) -> tuple[VectorExpr]:
        return self.args  # type: ignore[no-any-return]

    def __init__(self, *vectors: VectorExpr) -> None:
        # TODO: Add dispatch depending on the value of `vector` and `scale`?
        # TODO: if `vectors` is empty, return `ZERO` (in __new__?)

        # TODO: add dimension check for the arguments

        for vector in vectors:
            if not isinstance(vector, VectorExpr):
                raise TypeError(f"All addends must be VectorExpr, got {type(vector).__name__}")

        self._args = vectors

    def doit(self, **_hints: Any) -> VectorExpr:
        combination = self.as_symbol_combination()

        match combination:
            case ():
                return ZERO
            case ((v, s),):
                return v * s  # type: ignore[no-any-return]
            case _:
                return VectorAdd(*(v * s for v, s in combination))

    # TODO: order the addends in __init__? so that we could return a consistent `_hashable_contents`
    # tuple and get rid of custom __eq__
    def __eq__(self, other: object) -> bool:
        return type(other) is type(self) and self.args == other.args  # type: ignore[attr-defined]

    def __hash__(self) -> int:
        return hash(self.args)

    def _sympystr(self, p: Printer) -> str:
        return " + ".join(map(p.doprint, self.args))

    def as_symbol_combination(self) -> tuple[tuple[VectorSymbol, Expr], ...]:
        mapping: dict[VectorSymbol, Expr] = {}

        for addend in self.addends:
            for v, s in addend.as_symbol_combination():
                if v in mapping:
                    mapping[v] += s
                else:
                    mapping[v] = s

        for v, s in mapping.items():
            mapping[v] = simplify(s)

        combination = tuple((v, s) for v, s in mapping.items() if s != 0)
        return combination


__all__ = [
    "ZERO",
    "VectorAdd",
    "VectorExpr",
    "VectorNorm",
    "VectorScale",
    "VectorSymbol",
    "norm",
]
