from pytest import approx

APPROX_RELATIVE_TOLERANCE = 0.001


def assert_approx(lhs: float, rhs: float, *, tolerance: float = APPROX_RELATIVE_TOLERANCE):
    rhs_approx = approx(rhs, rel=tolerance, abs=abs(lhs * tolerance))
    assert lhs == rhs_approx
