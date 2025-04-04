# pylint: disable=too-many-lines

from __future__ import annotations

from typing import Any, Optional, TypeAlias, assert_never, Sequence, Self
from collections import defaultdict

from sympy import Atom, Basic, Expr, S, sympify
from sympy.core import function as sym_fn
from sympy.core.parameters import global_parameters
from sympy.physics.units import Dimension
from sympy.printing.printer import Printer

from symplyphysics.core.symbols.symbols import DimensionSymbol, next_name
from symplyphysics.core.symbols.id_generator import last_id
from symplyphysics.docs.miscellaneous import needs_mul_brackets
from ..miscellaneous import sort_with_sign, Registry, cacheit


class _AtomicRegistry:
    """
    Helper class that associates all `AtomicVectorExpr` with a number. This is needed to order
    the arguments in such operations as `VectorDot` or `VectorCross`.
    """

    _symbol_registry: Registry[VectorSymbol]
    _cross_registry: Registry[_VectorSymbolCross]
    _function_registry: Registry[AppliedVectorFunction]
    _derivative_registry: Registry[VectorDerivative]

    def __init__(self) -> None:
        self._symbol_registry = Registry()
        self._cross_registry = Registry()
        self._function_registry = Registry()
        self._derivative_registry = Registry()

    def add(self, value: AtomicVectorExpr) -> None:
        if isinstance(value, VectorSymbol):
            self._symbol_registry.add(value)
            return

        if isinstance(value, _VectorSymbolCross):
            self._cross_registry.add(value)
            return

        if isinstance(value, AppliedVectorFunction):
            self._function_registry.add(value)
            return

        if isinstance(value, VectorDerivative):
            self._derivative_registry.add(value)
            return

        assert_never(value)

    def get(self, value: AtomicVectorExpr) -> int:
        offset = self._offset(value)

        if isinstance(value, VectorSymbol):
            return offset + self._symbol_registry.get(value)

        if isinstance(value, _VectorSymbolCross):
            return offset + self._cross_registry.get(value)

        if isinstance(value, AppliedVectorFunction):
            return offset + self._function_registry.get(value)

        if isinstance(value, VectorDerivative):
            return offset + self._derivative_registry.get(value)

        assert_never(value)

    def _offset(self, value: AtomicVectorExpr) -> int:
        if isinstance(value, VectorSymbol):
            return 0

        # `_VectorSymbolCross` should come after `VectorSymbol`
        if isinstance(value, _VectorSymbolCross):
            return len(self._symbol_registry)

        if isinstance(value, AppliedVectorFunction):
            return len(self._symbol_registry) + len(self._cross_registry)

        if isinstance(value, VectorDerivative):
            return (len(self._symbol_registry) + len(self._cross_registry) +
                len(self._function_registry))

        assert_never(value)


_atomic_registry = _AtomicRegistry()


class VectorExpr(Basic):  # type: ignore[misc]
    """
    Base class for all vector expressions.
    """

    def doit(self, **_hints: Any) -> VectorExpr:
        return self

    def __add__(self, other: VectorExpr) -> VectorExpr:
        return VectorAdd(self, other)

    def __sub__(self, other: VectorExpr) -> VectorExpr:
        # Refer to the consequence in VectorAdd
        return VectorAdd(self, -other)

    def __mul__(self, other: Any) -> VectorExpr:
        return VectorScale(self, other)

    def __neg__(self) -> VectorExpr:
        # Refer to consequence #4 in VectorScale
        return VectorScale(self, -1)

    def __pos__(self) -> VectorExpr:
        return self

    def __truediv__(self, other: Any) -> VectorExpr:
        # NOTE: probably need to check if `other` is not `0` to return a special "NaN" vector

        return VectorScale(self, 1 / other)

    def subs(self, *args: Any, **kwargs: Any) -> VectorExpr:
        return Basic.subs(self, *args, **kwargs)  # type: ignore[no-any-return]

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        """
        Express `self` as a linear combination of `AtomicVectorExpr`. Each term is represented by a
        `(AtomicVectorExpr, Expr)` tuple, and each vector must appear in the combination only once.

        Examples:
        =========

        >>> from symplyphysics.core.experimental.coordinate_systems import CartesianCoordinateSystem
        >>> c = CartesianCoordinateSystem()
        >>> VectorAdd(c.i, VectorScale(c.j / S(2), S(4)), ZERO, -c.i * 2).as_symbol_combination()
        ((i, -1), (j, 2))
        """

        raise NotImplementedError(f"Implement this method in {type(self).__name__}.")

    def _eval_vector_derivative(self, _symbol: Expr) -> Optional[VectorExpr]:
        return None

    def _eval_vector_norm(self) -> Optional[Expr]:
        return None


class _VectorZero(VectorExpr, Atom):  # type: ignore[misc]
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

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        return ()

    def _eval_vector_derivative(self, _symbol: Expr) -> VectorExpr:
        return ZERO

    def _eval_vector_norm(self) -> Expr:
        # Refer to property #3 in VectorNorm
        return S.Zero


ZERO = _VectorZero()


def _process_vector_names(
    code: Optional[str] = None,
    latex: Optional[str] = None,
    *,
    base: str = "VEC",
    i: Optional[int] = None,
) -> tuple[str, str]:
    if not code:
        code = f"{base}{i}"
        if not latex:
            base = base[0].lower()
            latex = f"\\mathbf{{{base}}}_{{{i}}}"
    elif not latex:
        latex = f"\\mathbf{{{code}}}"

    return code, latex


# TODO: Add support for axial vectors.
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

    is_symbol = True

    def __new__(
            cls,
            display_symbol: Optional[str] = None,
            dimension: Dimension = Dimension(1),
            *,
            display_latex: Optional[str] = None,
    ) -> VectorExpr:
        obj = super().__new__(cls)

        _atomic_registry.add(obj)

        return obj  # type: ignore[no-any-return]

    def __init__(
            self,
            display_symbol: Optional[str] = None,
            dimension: Dimension = Dimension(1),
            *,
            display_latex: Optional[str] = None,
    ):
        id_ = _atomic_registry.get(self)

        display_symbol, display_latex = _process_vector_names(display_symbol, display_latex, i=id_)

        DimensionSymbol.__init__(
            self,
            display_name=display_symbol,
            dimension=dimension,
            display_latex=display_latex,
        )
        VectorExpr.__init__(self)
        Atom.__init__(self)

    def _hashable_content(self) -> tuple[Any, ...]:
        return (id(self),)

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        return ((self, S.One),)

    def _eval_vector_derivative(self, _symbol: Expr) -> VectorExpr:
        return ZERO


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

    @cacheit
    def __new__(cls, vector: VectorExpr, **kwargs: Any) -> Expr:
        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            obj = cls.from_vector(vector)
            if obj is not None:
                return obj

        obj = super().__new__(cls)
        obj._args = (vector,)
        return obj

    @classmethod
    def from_vector(cls, vector: VectorExpr) -> Optional[Expr]:
        result = vector._eval_vector_norm()  # pylint: disable=protected-access
        if result is not None:
            return result

        if isinstance(vector, VectorScale):
            return cls(vector.vector) * abs(vector.scale)

        # TODO: add support for the following relation:
        # for all vectors `a, b` and scalars `k`, `norm(a * k + b * k) = norm(a + b) * abs(k)`

        return None

    def doit(self, **hints: Any) -> Expr:
        vector = self.argument

        if hints.get("deep", True):
            vector = vector.doit(**hints)

        result = self.from_vector(vector)

        if result is not None:
            return result

        return self

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return f"norm({p.doprint(self.argument)})"

    def _eval_derivative(self, symbol: Expr) -> Expr:
        done = self.doit()

        if isinstance(done, VectorNorm):
            # (d/dx)[norm(v(x))]
            # = (d/dx)[sqrt(dot(v(x), v(x)))]
            # = (d/dx)[dot(v(x), v(x))] / (2 * sqrt(v(x), v(x)))
            # = (2 * dot(v(x), (d/dx)[v(x)])) / (2 * norm(v(x)))
            # = dot(v(x), (d/dx)[v(x)]) / norm(v(x))

            vector = done.argument
            return VectorDot(vector, VectorDerivative(vector, symbol)) / done

        return done.diff(symbol)


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

    @cacheit
    def __new__(cls, vector: VectorExpr, scale: Any, **kwargs: Any) -> VectorScale:
        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            return cls.from_arguments(vector, scale)

        obj = super().__new__(cls)
        obj._args = (vector, scale)
        return obj  # type: ignore[no-any-return]

    @classmethod
    def from_arguments(cls, vector: VectorExpr, scale: Any) -> VectorExpr:
        scale = sympify(scale, strict=True)
        if not isinstance(scale, Expr):
            raise TypeError(f"Scale {scale} must be an Expr, got {type(scale).__name__}.")

        combination = [(v, s * scale) for v, s in vector.as_symbol_combination()]

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

                return cls(v, s, evaluate=False)
            case _:
                # Refer to property #3
                return VectorAdd(
                    *(cls(v, s, evaluate=False) for v, s in combination),
                    evaluate=False,
                )

    def doit(self, **_hints: Any) -> VectorExpr:
        return self.from_arguments(self.vector, self.scale)

    def _sympystr(self, p: Printer) -> str:
        vector, value = self.args

        if needs_mul_brackets(value, last=True):
            return f"{p.doprint(vector)}*({p.doprint(value)})"

        return f"{p.doprint(vector)}*{p.doprint(value)}"

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        scale = self.scale
        vector = self.vector

        result: list[tuple[AtomicVectorExpr, Expr]] = []
        for v, s0 in vector.as_symbol_combination():
            s1 = s0 * scale
            if s1 != 0:
                result.append((v, s1))
        return tuple(result)

    def _eval_vector_derivative(self, symbol: Expr) -> VectorExpr:
        vector = self.vector
        scale = self.scale

        derived_vector = VectorDerivative(vector, symbol) * scale
        derived_scaled = vector * scale.diff(symbol)

        return derived_vector + derived_scaled  # type: ignore[no-any-return]

    def _eval_vector_norm(self) -> Expr:
        # Refer to property #2 in VectorNorm
        return VectorNorm(self.vector) * abs(self.scale)


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

    @cacheit
    def __new__(cls, *vectors: VectorExpr, **kwargs: Any) -> VectorExpr:
        # NOTE: add dimension check for the arguments?

        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            for vector in vectors:
                if not isinstance(vector, VectorExpr):
                    raise TypeError(f"All addends must be VectorExpr, got {type(vector).__name__}")

            return cls.from_vectors(*vectors)

        obj = super().__new__(cls)
        obj._args = vectors
        return obj  # type: ignore[no-any-return]

    @classmethod
    def from_vectors(cls, *vectors: VectorExpr) -> VectorExpr:
        if not vectors:
            return ZERO

        combination = cls._as_symbol_combination(*vectors)

        match combination:
            case ():
                return ZERO
            case ((v, s),):
                return v * s  # type: ignore[no-any-return]
            case _:
                return cls(*(v * s for v, s in combination), evaluate=False)

    def doit(self, **_hints: Any) -> VectorExpr:
        return self.from_vectors(*self.addends)

    def _sympystr(self, p: Printer) -> str:
        return " + ".join(map(p.doprint, self.args))

    @staticmethod
    def _as_symbol_combination(*addends: VectorExpr) -> tuple[tuple[VectorExpr, Expr], ...]:
        mapping: dict[AtomicVectorExpr, Expr] = defaultdict(lambda: S.Zero)

        for addend in addends:
            for v, s in addend.as_symbol_combination():
                mapping[v] += s

        return tuple((v, s) for v, s in mapping.items() if s != 0)

    def as_symbol_combination(self) -> tuple[tuple[VectorSymbol, Expr], ...]:
        return self._as_symbol_combination(*self.addends)

    def _eval_vector_derivative(self, symbol: Expr) -> VectorExpr:
        derived_args = [VectorDerivative(arg, symbol) for arg in self.args]
        return VectorAdd(*derived_args)


class VectorDot(Expr):  # type: ignore[misc]
    """
    The **dot product**, or **scalar product**, is a binary operation that takes two vectors and
    returns a single number.

    Geometrically, the dot product can be expressed using the length (`norm`) of the vectors and
    the (non-directional) `angle` between them: `dot(a, b) = norm(a) * norm(b) * cos(angle(a, b))`
    where `a` and `b` are vectors.

    In particular,

    1. `dot(a, b) = 0` if and only if `a` and `b` are orthogonal.

    2. `dot(a, b) = norm(a) * norm(b)` if and only if `a` and `b` are codirectional.

    3. `dot(a, a) = norm(a)^2` as a result of (2).

    The dot product has the following **properties**:

    1. **Commutativity**: for all vectors `a, b`, `dot(a, b) = dot(b, a)`.

    2. **Linearity** in the **first** argument: for all vectors `a, b, c` and scalars `k, l`,
       `dot(a * k + b * l, c) = k * dot(a, c) + l * dot(b, c)`.

    3. **Linearity** in the **second** argument, which follows from (1) and (2).

    4. **Absense of cancellation**: for all vectors `a, b, c` s.t. `a ≠ 0`, `dot(a, b) = dot(a, c)`
       does not imply `b = c`.

    5. Applicability of the **product rule**: for all vector-valued differentiable functions
       `a, b`, `d[dot(a, b)] = dot(d[a], b) + dot(a, d[b])` where `d[v]` denotes the derivative of
       vector `v`.

    Note that the properties and relations mentioned only apply to **real-valued vectors**.

    The dot product is a *true scalar* in a sense that it is unchanged if the orientation of the
    frame is reversed.
    """

    @property
    def lhs(self) -> VectorExpr:
        return self.args[0]  # type: ignore[no-any-return]

    @property
    def rhs(self) -> VectorExpr:
        return self.args[1]  # type: ignore[no-any-return]

    @cacheit
    def __new__(cls, lhs: VectorExpr, rhs: VectorExpr, **kwargs: Any) -> Expr:
        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            return cls.from_vectors(lhs, rhs)

        obj = super().__new__(cls)
        obj._args = (lhs, rhs)
        return obj

    def doit(self, **_hints: Any) -> Expr:
        return self.from_vectors(self.lhs, self.rhs)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        return f"dot({p.doprint(self.lhs)}, {p.doprint(self.rhs)})"

    @classmethod
    def from_vectors(cls, lhs: VectorExpr, rhs: VectorExpr) -> Expr:
        if isinstance(lhs, AtomicVectorExpr) and isinstance(rhs, AtomicVectorExpr):
            return cls.from_atomic(lhs, rhs)

        result = S.Zero

        for lhs_v, lhs_s in lhs.as_symbol_combination():
            for rhs_v, rhs_s in rhs.as_symbol_combination():
                result += cls.from_atomic(lhs_v, rhs_v) * lhs_s * rhs_s

        return result

    @classmethod
    def from_symbols(
        cls,
        lhs: VectorSymbol | AppliedVectorFunction,
        rhs: VectorSymbol | AppliedVectorFunction,
    ) -> Expr:
        sign, args = sort_with_sign((lhs, rhs), key=_atomic_registry.get)

        if sign == 0:
            return VectorNorm(lhs)**2

        return cls(*args, evaluate=False)

    @classmethod
    def from_atomic(cls, lhs: AtomicVectorExpr, rhs: AtomicVectorExpr) -> Expr:
        """
        1. `dot(v, w)` is left unchanged unless `v = w`, in which case `dot(v, v) = norm(v)**2`.

        2. `dot(v, cross(c, d)) = mixed(v, c, d)`.

        3. `dot(cross(a, b), w) = mixed(w, a, b)`.

        4. `dot(cross(a, b), cross(c, d)) = dot(a, b) * dot(c, d) - dot(b, c) * dot(a, d)`.
        """

        if isinstance(lhs, _VectorSymbolCross):
            if isinstance(rhs, _VectorSymbolCross):
                a = lhs.lhs
                b = lhs.rhs
                c = rhs.lhs
                d = rhs.rhs

                return (cls.from_symbols(a, b) * cls.from_symbols(b, d) -
                    cls.from_symbols(b, c) * cls.from_symbols(a, d))

            return VectorMixedProduct.from_symbols(rhs, lhs.lhs, lhs.rhs)

        if isinstance(rhs, _VectorSymbolCross):
            return VectorMixedProduct.from_symbols(lhs, rhs.lhs, rhs.rhs)

        return cls.from_symbols(lhs, rhs)

    def _eval_derivative(self, symbol: Expr) -> Expr:
        lhs = self.lhs
        rhs = self.rhs

        derived_lhs = VectorDot(VectorDerivative(lhs, symbol), rhs)
        derived_rhs = VectorDot(lhs, VectorDerivative(rhs, symbol))

        return derived_lhs + derived_rhs


class VectorCross(VectorExpr):
    """
    The **cross product**, or **vector product**, is a binary operation that takes two vectors and
    returns another vector. The cross product is only defined in a *3-dimensional space* (although
    its construction is also possible in a 7-dimensional space, the following properties do not
    hold there).

    Geometrically, the cross product between vectors `a` and `b` can be defined as `cross(a, b) =
    norm(a) * norm(b) * sin(angle(a, b)) * n` where `norm` is the length operator and `n` is a unit
    vector orthogonal to the plane containing `a` and `b` with such direction that the ordered set
    `(a, b, n)` is positively oriented.

    The cross product has the following **properties**:

    1. **Anticommutativity**: for all vectors `a, b`, `cross(a, b) = -cross(b, a)`.

    2. For any vector `a`, `cross(a, a) = 0`, which follows from (1).

    3. **Distributivity over addition**: for all vectors `a, b, c`, `cross(a, b + c) =
       cross(a, b) + cross(a, c)`.

    4. **Absense of associativity**: for all vectors `a, b, c`, `cross(a, cross(b, c)) ≠
       cross(cross(a, b), c)`.

    5. However, the **Jacobi identity** is satisfied: for all vectors `a, b, c`,
       `cross(a, cross(b, c)) + cross(b, cross(c, a)) + cross(c, cross(b, a)) = 0`.

    6. **Absense of cancellation**: for all vectors `a, b, c` s.t. `a ≠ 0`, `cross(a, b) =
       cross(a, c)` does not imply `b = c`. This only happens if `dot(a, b) = dot(a, c)` holds.

    7. Applicability of the **product rule**: for all vector-valued differentiable functions
       `a, b`, `d[cross(a, b)] = cross(d[a], b) + cross(a, d[b])` where `d[v]` denotes the
       derivative of vector `v`.

    It is related to the dot product by the following relation: for all vectors `a, b`,
    `norm(cross(a, b))^2 + dot(a, b)^2 = norm(a)^2 * norm(b)^2`.

    The cross product is a *pseudovector*, i.e. it is negated if the orientation of the frame is
    reversed.
    """

    @property
    def lhs(self) -> VectorExpr:
        return self.args[0]  # type: ignore[no-any-return]

    @property
    def rhs(self) -> VectorExpr:
        return self.args[1]  # type: ignore[no-any-return]

    @cacheit
    def __new__(cls, lhs: VectorExpr, rhs: VectorExpr, **kwargs: Any) -> VectorExpr:
        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            return cls.from_vectors(lhs, rhs)

        obj = super().__new__(cls)
        obj._args = (lhs, rhs)
        return obj  # type: ignore[no-any-return]

    @classmethod
    def from_vectors(cls, lhs: VectorExpr, rhs: VectorExpr) -> VectorExpr:
        if isinstance(lhs, VectorSymbol) and isinstance(rhs, VectorSymbol):
            return _VectorSymbolCross(lhs, rhs)

        combination = cls._as_symbol_combination(lhs, rhs)

        return VectorAdd(
            *(VectorScale(v, s, evaluate=False) for v, s in combination),
            evaluate=False,
        )

    def doit(self, **_hints: Any) -> VectorExpr:
        return self.from_vectors(self.lhs, self.rhs)

    @classmethod
    def _as_symbol_combination(
        cls,
        lhs: VectorExpr,
        rhs: VectorExpr,
    ) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        mapping: dict[AtomicVectorExpr, Expr] = defaultdict(lambda: S.Zero)

        for lhs_v, lhs_s in lhs.as_symbol_combination():
            for rhs_v, rhs_s in rhs.as_symbol_combination():
                new_factor = lhs_s * rhs_s

                for v, s in cls.from_atomic(lhs_v, rhs_v).as_symbol_combination():
                    mapping[v] += s * new_factor

        return tuple((v, s) for v, s in mapping.items() if s != 0)

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        return self._as_symbol_combination(self.lhs, self.rhs)

    def _sympystr(self, p: Printer) -> str:
        return f"cross({p.doprint(self.lhs)}, {p.doprint(self.rhs)})"

    @classmethod
    def from_atomic(cls, lhs: AtomicVectorExpr, rhs: AtomicVectorExpr) -> VectorExpr:
        """
        1. `cross(v, w)` is left unchanged.

        2. `cross(v, cross(c, d)) = c * dot(v, d) - d * dot(v, c)`, referred to as the
            Lagrange's formula.

        3. `cross(cross(a, b), w) = b * dot(w, a) - a * dot(w, b)` due to anticommutativity of
            the cross product.

        4. `cross(cross(a, b), cross(c, d)) = c * mixed(d, a, b) - d * mixed(c, a, b)` from the
            Lagrange's formula.

        **Links:**

        1. `Lagrange's formula <https://en.wikipedia.org/wiki/Triple_product#Vector_triple_product>`__.
        """

        if isinstance(lhs, _VectorSymbolCross):
            a = lhs.lhs
            b = lhs.rhs

            if isinstance(rhs, _VectorSymbolCross):
                c = rhs.lhs
                d = rhs.rhs

                # Refer to formula #4
                return (  # type: ignore[no-any-return]
                    c * VectorMixedProduct.from_symbols(d, a, b) -
                    d * VectorMixedProduct.from_symbols(c, a, b))

            # Refer to formula #3
            return (  # type: ignore[no-any-return]
                b * VectorDot.from_atomic(rhs, a) - a * VectorDot.from_atomic(rhs, b))

        if isinstance(rhs, _VectorSymbolCross):
            c, d = rhs.args

            # Refer to formula #2
            return (  # type: ignore[no-any-return]
                c * VectorDot.from_atomic(lhs, d) - d * VectorDot.from_atomic(lhs, c))

        # Refer to formula #1
        return _VectorSymbolCross.from_symbols(lhs, rhs)

    def _eval_vector_derivative(self, symbol: Expr) -> VectorExpr:
        lhs = self.lhs
        rhs = self.rhs

        derived_lhs = VectorCross(VectorDerivative(lhs, symbol), rhs)
        derived_rhs = VectorCross(lhs, VectorDerivative(rhs, symbol))
        return derived_lhs + derived_rhs


class _VectorSymbolCross(VectorCross):
    """
    A helper class for the cross product between two symbolic vectors.
    """

    @property
    def lhs(self) -> VectorSymbol | AppliedVectorFunction | VectorDerivative:
        return self.args[0]  # type: ignore[no-any-return]

    @property
    def rhs(self) -> VectorSymbol | AppliedVectorFunction | VectorDerivative:
        return self.args[1]  # type: ignore[no-any-return]

    @cacheit
    def __new__(
        cls,
        lhs: VectorSymbol | AppliedVectorFunction | VectorDerivative,
        rhs: VectorSymbol | AppliedVectorFunction | VectorDerivative,
        **kwargs: Any,
    ) -> _VectorSymbolCross:
        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            return cls.from_symbols(lhs, rhs)

        obj = super().__new__(cls, lhs, rhs, evaluate=False)
        _atomic_registry.add(obj)
        return obj

    def doit(self, **_hints: Any) -> VectorExpr:
        return self.from_symbols(self.lhs, self.rhs)

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        sign, args = sort_with_sign((self.lhs, self.rhs), key=_atomic_registry.get)

        if sign == 0:
            return ()

        if sign == -1:
            result = type(self)(*args, evaluate=False)
        else:
            result = self

        return ((result, S(sign)),)

    @classmethod
    def from_symbols(
        cls,
        lhs: VectorSymbol | AppliedVectorFunction | VectorDerivative,
        rhs: VectorSymbol | AppliedVectorFunction | VectorDerivative,
    ) -> VectorExpr:
        sign, args = sort_with_sign((lhs, rhs), key=_atomic_registry.get)

        if sign == 0:
            return ZERO

        return cls(*args, evaluate=False) * sign


class VectorMixedProduct(Expr):  # type: ignore[misc]
    """
    The **scalar triple product**, or **mixed product**, is defined as the dot product of one
    vector with the cross product with the other two.

    Geometrically, given vectors `a, b, c` the scalar triple product `dot(a, cross(b, c))` can be
    interpreted as the signed volume of the parallelepiped defined by these vectors.

    The mixed product has the following **properties**:

    1. It does not change under a positive permutation of the arguments: for all vectors `a, b, c`,
       `mixed(a, b, c) = mixed(b, c, a) = mixed(c, a, b)`.

    2. It is negated under a negative permutation of the arguments: for all vectors `a, b, c`,
       `mixed(a, c, b) = mixed(b, a, c) = mixed(c, b, a) = -mixed(a, b, c)`.

    3. The scalar triple product is zero if and only if the three vectors `a, b, c` are coplanar.

    4. It is linear in all arguments.

    The scalar triple product is a *pseudoscalar*, i.e. it is negated if the orientation of the
    frame is reversed.
    """

    is_real = True

    @property
    def vectors(self) -> tuple[VectorExpr, VectorExpr, VectorExpr]:
        a, b, c = self.args
        return a, b, c

    @cacheit
    def __new__(cls, *vectors: VectorExpr, **kwargs: Any) -> Expr:
        a, b, c = vectors

        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            return cls.from_vectors(*vectors)

        obj = super().__new__(cls)
        obj._args = (a, b, c)
        return obj

    @classmethod
    def from_vectors(cls, *vectors: VectorExpr) -> Expr:
        if all(isinstance(v, VectorSymbol) for v in vectors):
            return cls.from_symbols(*vectors)

        a, b, c = vectors
        return VectorDot(a, VectorCross(b, c))

    def doit(self, **_hints: Any) -> Expr:
        return self.from_vectors(*self.vectors)

    def _eval_nseries(self, x: Any, n: Any, logx: Any, cdir: Any) -> Any:
        pass

    def _sympystr(self, p: Printer) -> str:
        a, b, c = self.args
        return f"mixed({p.doprint(a)}, {p.doprint(b)}, {p.doprint(c)})"

    @classmethod
    def from_symbols(cls, *vectors: VectorSymbol | AppliedVectorFunction) -> Expr:
        sign, sorted_args = sort_with_sign(vectors, key=_atomic_registry.get)

        if sign == 0:
            return S.Zero

        return sign * cls(*sorted_args, evaluate=False)

    def _eval_derivative(self, symbol: Expr) -> Expr:
        return self.doit().diff(symbol)


class AppliedVectorFunction(sym_fn.Application, VectorExpr):  # type: ignore[misc]
    """This class represents the result of applying a vector-valued function to some arguments."""

    @cacheit
    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        if cls is AppliedVectorFunction:
            raise TypeError("Call `VectorFunction` instead to instantiate a new function.")

        n = len(args)

        if not cls._valid_nargs(n):
            template = "{name} takes {qual} {args} argument{plural} ({given} given)"
            arguments = getattr(cls, "arguments", None)
            nargs = min(cls.nargs) if arguments is None else len(arguments)
            message = template.format(
                name=cls,
                qual="exactly" if len(cls.nargs) == 1 else "at least",
                args=nargs,
                plural="s" * (nargs != 1),
                given=n,
            )
            raise TypeError(message)

        args = tuple(sympify(arg, strict=True) for arg in args)
        undefineds = [arg.name for arg in args if isinstance(arg, sym_fn.FunctionClass)]
        if undefineds:
            template = "Invalid argument: expecting an expression, not undefined function{plural}: {names}"
            message = template.format(
                plural="s" * (len(undefineds) > 1),
                names=", ".join(undefineds),
            )
            raise TypeError(message)

        result = super().__new__(cls, *args, **kwargs)
        return result  # type: ignore[no-any-return]

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        return ((self, S.One),)

    def _eval_vector_derivative(self, symbol: Expr) -> VectorExpr:
        # For now, the derivative of a function application is only implemented when no more than
        # one argument of the function is equal to `symbol`. All other cases require a vector
        # alternative of `sympy.Subs` or a Jacobian matrix.

        derived_args = []

        found = False

        for arg in self.args:
            if isinstance(arg, Expr):
                scale = arg.diff(symbol)

                if scale == 0:
                    continue

                if scale != 1 or found:
                    raise NotImplementedError("VectorSubs has not been implemented yet.")

                found = True

                partial = VectorDerivative(self, arg, evaluate=False)
                derived_args.append(partial * scale)
            elif isinstance(arg, VectorExpr):
                # Assuming `v` and `w` are vector-valued functions:
                #   (d/dx)[v(w(x))] = jacobian[v](w(x)) * (d/dx)[w(x)]
                # Here, `jacobian[v](w(x))` denotes the Jacobian matrix of `v` evaluated at `w(x)`.

                raise NotImplementedError("Jacobian matrix has not been implemented yet.")
            else:
                raise TypeError(
                    f"The argument must be `Expr` or `VectorExpr`, got {type(arg).__name__}.")

        if not found:
            return ZERO

        return VectorAdd(*derived_args)


class UndefinedVectorFunction(sym_fn.FunctionClass):  # type: ignore[misc]
    """The (meta)class of undefined vector functions."""

    def __new__(
        mcs,
        name: str,
        bases: Optional[Sequence[type]] = None,
        __dict__: Optional[dict[str, Any]] = None,
        **kwargs: Any,
    ) -> Self:
        bases = bases or (AppliedVectorFunction,)

        if __dict__ is None:
            __dict__ = {}
        __dict__ |= kwargs
        __dict__["_kwargs"] = kwargs
        __dict__["__module__"] = None

        obj = super().__new__(mcs, name, bases, __dict__)
        obj.name = name
        return obj  # type: ignore[no-any-return]


class VectorFunction(DimensionSymbol, UndefinedVectorFunction):
    _arguments: Optional[tuple[Basic]]

    @property
    def arguments(cls) -> Optional[tuple[Basic]]:
        return cls._arguments

    def __new__(  # pylint: disable=signature-differs
        mcs,
        display_name: Optional[str] = None,
        arguments: Optional[Sequence[Basic]] = None,
        *,
        dimension: Dimension = Dimension(1),
        display_latex: Optional[str] = None,
        **kwargs: Any,
    ) -> Self:
        name = next_name("FUN")
        obj = UndefinedVectorFunction.__new__(mcs, name, **kwargs)

        return obj

    def __init__(
        cls,
        display_name: Optional[str] = None,
        arguments: Optional[Sequence[Basic]] = None,
        *,
        dimension: Dimension = Dimension(1),
        display_latex: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        if arguments is not None and not isinstance(arguments, tuple):
            arguments = tuple(arguments)
        cls._arguments = arguments

        display_name, display_latex = _process_vector_names(
            display_name,
            display_latex,
            base="FUN",
            i=last_id("FUN"),
        )
        DimensionSymbol.__init__(
            cls,
            display_name=display_name,
            dimension=dimension,
            display_latex=display_latex,
        )

        if arguments is not None:
            kwargs["nargs"] = len(arguments)
        UndefinedVectorFunction.__init__(cls, **kwargs)

    def __repr__(cls) -> str:
        return str(cls.display_name)


class VectorDerivative(VectorExpr):

    @cacheit
    def __new__(cls, vector: VectorExpr, symbol: Expr, **kwargs: Any) -> VectorExpr:
        if not isinstance(vector, VectorExpr):
            raise TypeError("VectorDerivative only accepts a vector expression.")

        symbol = sympify(symbol, strict=True)

        if not isinstance(symbol, Expr):
            raise TypeError(
                f"VectorDerivative only accepts a single differentiation symbol, got {type(symbol).__name__}."
            )

        # Check that `symbol` is not a number since differentiation is only defined for
        # non-constant expressions

        is_number = True

        try:
            _ = complex(symbol)
        except TypeError:
            is_number = False

        if is_number:
            raise ValueError("The differentiation symbol must not be a number.")

        evaluate = kwargs.get("evaluate", global_parameters.evaluate)

        if evaluate:
            result = vector._eval_vector_derivative(symbol)
            if result is not None:
                return result

        obj = super().__new__(cls)
        obj._args = (vector, symbol)

        _atomic_registry.add(obj)
        return obj  # type: ignore[no-any-return]

    def _eval_vector_derivative(self, symbol: Expr) -> VectorExpr:
        # TODO: probably check for inner derivatives

        return self.func(self, symbol)  # type: ignore[no-any-return]

    def doit(self, **_hints: Any) -> VectorExpr:
        return VectorDerivative(*self.args)

    def as_symbol_combination(self) -> tuple[tuple[AtomicVectorExpr, Expr], ...]:
        # TODO: since `self` can contain an unevaluated derivative, re-evaluate it without making
        # an infinite recursion loop

        return ((self, S.One),)


AtomicVectorExpr: TypeAlias = VectorSymbol | _VectorSymbolCross | AppliedVectorFunction | VectorDerivative

__all__ = [
    "ZERO",
    "AtomicVectorExpr",
    "VectorAdd",
    "VectorCross",
    "VectorDot",
    "VectorExpr",
    "VectorMixedProduct",
    "VectorNorm",
    "VectorScale",
    "VectorSymbol",
    "AppliedVectorFunction",
    "VectorFunction",
]
