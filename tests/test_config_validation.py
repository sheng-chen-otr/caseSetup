"""Unit tests for configValidation (suspension hardpoint sanity checks)."""
import configValidation


def _corner(center, scale):
    """Build a plausible corner dict at a given unit scale.

    ``scale`` multiplies a nominal metric geometry (~0.2 m spans), so scale=1.0
    is metres and scale=1000.0 is millimetres.
    """
    def s(vec):
        return [v * scale for v in vec]

    return {
        "wheel_center_static": s(center),
        "uca_f_inner": s([center[0] + 0.05, center[1] - 0.20, center[2] + 0.15]),
        "uca_r_inner": s([center[0] - 0.05, center[1] - 0.20, center[2] + 0.15]),
        "lca_f_inner": s([center[0] + 0.06, center[1] - 0.22, center[2] - 0.12]),
        "lca_r_inner": s([center[0] - 0.06, center[1] - 0.22, center[2] - 0.12]),
        "tie_inner": s([center[0] + 0.08, center[1] - 0.20, center[2] - 0.02]),
        "uca_outer_static": s([center[0], center[1] - 0.01, center[2] + 0.10]),
        "lca_outer_static": s([center[0], center[1] - 0.01, center[2] - 0.10]),
        "tie_outer_static": s([center[0] + 0.05, center[1] - 0.01, center[2]]),
        "pushrod_outer_static": s([center[0], center[1] - 0.01, center[2] + 0.08]),
    }


class TestCornerCharacteristicSpan:
    def test_returns_positive_span(self):
        c = _corner([1.5, 0.75, 0.30], scale=1.0)
        span = configValidation.cornerCharacteristicSpan(c)
        assert span is not None and span > 0.0

    def test_missing_center_returns_none(self):
        c = _corner([1.5, 0.75, 0.30], scale=1.0)
        del c["wheel_center_static"]
        assert configValidation.cornerCharacteristicSpan(c) is None

    def test_scale_is_proportional(self):
        c_m = _corner([1.5, 0.75, 0.30], scale=1.0)
        c_mm = _corner([1.5, 0.75, 0.30], scale=1000.0)
        span_m = configValidation.cornerCharacteristicSpan(c_m)
        span_mm = configValidation.cornerCharacteristicSpan(c_mm)
        assert span_mm / span_m == __import__("pytest").approx(1000.0, rel=1e-9)


class TestCheckHardpointUnitConsistency:
    def test_consistent_metres_ok(self):
        corners = {
            "fl": _corner([1.5, 0.75, 0.30], 1.0),
            "fr": _corner([1.5, -0.75, 0.30], 1.0),
            "rl": _corner([-1.5, 0.75, 0.30], 1.0),
            "rr": _corner([-1.5, -0.75, 0.30], 1.0),
        }
        ok, msg = configValidation.checkHardpointUnitConsistency(corners)
        assert ok
        assert msg == ""

    def test_consistent_millimetres_ok(self):
        corners = {
            "fl": _corner([1.5, 0.75, 0.30], 1000.0),
            "fr": _corner([1.5, -0.75, 0.30], 1000.0),
            "rl": _corner([-1.5, 0.75, 0.30], 1000.0),
            "rr": _corner([-1.5, -0.75, 0.30], 1000.0),
        }
        ok, msg = configValidation.checkHardpointUnitConsistency(corners)
        assert ok

    def test_mixed_units_flagged(self):
        """FL in mm, others in m -> ~1000x span ratio -> flagged.

        Regresses the recurring template bug (FL mm, FR/RL/RR m, uniform scale).
        """
        corners = {
            "fl": _corner([1.5, 0.75, 0.30], 1000.0),
            "fr": _corner([1.5, -0.75, 0.30], 1.0),
            "rl": _corner([-1.5, 0.75, 0.30], 1.0),
            "rr": _corner([-1.5, -0.75, 0.30], 1.0),
        }
        ok, msg = configValidation.checkHardpointUnitConsistency(corners)
        assert not ok
        assert "mismatched units" in msg
        assert "HARDPOINT_SCALE" in msg

    def test_single_corner_cannot_judge(self):
        corners = {"fl": _corner([1.5, 0.75, 0.30], 1.0)}
        ok, msg = configValidation.checkHardpointUnitConsistency(corners)
        assert ok
