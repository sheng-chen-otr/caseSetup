"""Unit tests for the pure geometry / kinematics / cornering math in utilities.py.

These functions are self-contained (numpy + dict inputs) and require no
OpenFOAM installation. Each test group maps to a behaviour that has previously
regressed and been validated by hand (see repo memory notes), so these serve as
regression guards.
"""
import math

import numpy as np
import pytest

import utilities


# ---------------------------------------------------------------------------
# _build_block_rigid_transform / _transform_vertex_block
# ---------------------------------------------------------------------------
class TestBlockRigidTransform:
    def test_none_cfg_and_pivot_is_identity_no_nan(self):
        """cfg=None, pivot=None must yield identity with NO NaN.

        Regresses the static-PID NaN corruption: pivot=None previously became
        np.array(None) -> NaN and destroyed unmatched (static) blocks.
        """
        R, p, d = utilities._build_block_rigid_transform(None, None)
        assert np.allclose(R, np.eye(3))
        assert np.allclose(p, [0.0, 0.0, 0.0])
        assert np.allclose(d, [0.0, 0.0, 0.0])
        assert not np.isnan(p).any()
        assert not np.isnan(d).any()

    def test_none_cfg_leaves_block_unchanged(self):
        block = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        out = utilities._transform_vertex_block(None, None, block, np.arange(2))
        assert np.allclose(out, block)
        assert not np.isnan(out).any()

    def test_pure_translation(self):
        cfg = {"translation": [10.0, -2.0, 0.5]}
        block = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        out = utilities._transform_vertex_block(cfg, [0, 0, 0], block, np.arange(2))
        assert np.allclose(out, block + np.array([10.0, -2.0, 0.5]))

    def test_rotation_90_about_z(self):
        cfg = {"rotation": {"x": 0.0, "y": 0.0, "z": 90.0}}
        block = np.array([[1.0, 0.0, 0.0]])
        out = utilities._transform_vertex_block(cfg, [0, 0, 0], block, np.arange(1))
        # +90 deg about z takes (1,0,0) -> (0,1,0)
        assert np.allclose(out[0], [0.0, 1.0, 0.0], atol=1e-12)

    def test_rotation_about_pivot(self):
        """Rotation is applied about the pivot, not the world origin."""
        cfg = {"rotation": {"x": 0.0, "y": 0.0, "z": 90.0}, "pivot": [1.0, 0.0, 0.0]}
        block = np.array([[2.0, 0.0, 0.0]])
        out = utilities._transform_vertex_block(cfg, [1.0, 0.0, 0.0], block, np.arange(1))
        # point 1 unit +x of pivot, rotate +90 about z -> 1 unit +y of pivot
        assert np.allclose(out[0], [1.0, 1.0, 0.0], atol=1e-12)

    def test_block_matches_per_vertex(self):
        """Vectorized block transform equals per-vertex transform (rigid)."""
        cfg = {
            "translation": [0.3, -1.2, 4.0],
            "rotation": {"x": 12.0, "y": -30.0, "z": 45.0},
            "pivot": [0.5, 0.5, 0.5],
        }
        rng = np.random.default_rng(0)
        block = rng.normal(size=(50, 3))
        idx = np.arange(block.shape[0])
        blockOut = utilities._transform_vertex_block(cfg, [0, 0, 0], block, idx)

        R, p, d = utilities._build_block_rigid_transform(cfg, [0, 0, 0])
        perVertex = np.vstack([(v - p) @ R.T + p + d for v in block])
        assert np.allclose(blockOut, perVertex, atol=1e-12)


# ---------------------------------------------------------------------------
# calcLoadedRadius
# ---------------------------------------------------------------------------
def _bc_dict(refcor=(0.0, 0.0, 0.0), pitch=0.0, roll=0.0):
    return {
        "BC_SETUP": {
            "REFCOR": [str(refcor[0]), str(refcor[1]), str(refcor[2])],
            "DOMAIN_PITCH": [str(pitch)],
            "DOMAIN_ROLL": [str(roll)],
        }
    }


class TestCalcLoadedRadius:
    def test_zero_tilt_reduces_to_zcenter_minus_refz(self):
        d = _bc_dict(refcor=(1.0, 2.0, 0.25), pitch=0.0, roll=0.0)
        r = utilities.calcLoadedRadius(5.0, -3.0, 0.6, d)
        assert r == pytest.approx(0.6 - 0.25)

    def test_matches_plane_formula_with_tilt(self):
        refcor = (0.5, -0.5, 0.1)
        pitch_deg, roll_deg = 3.0, -2.0
        d = _bc_dict(refcor=refcor, pitch=pitch_deg, roll=roll_deg)
        xc, yc, zc = 1.2, 0.3, 0.55

        pitch = math.radians(pitch_deg)
        roll = math.radians(roll_deg)
        nx = math.sin(pitch) * math.cos(roll)
        ny = -math.sin(roll)
        nz = math.cos(pitch) * math.cos(roll)
        groundZ = refcor[2] - (nx * (xc - refcor[0]) + ny * (yc - refcor[1])) / nz
        expected = zc - groundZ

        assert utilities.calcLoadedRadius(xc, yc, zc, d) == pytest.approx(expected)


# ---------------------------------------------------------------------------
# corneringAxis
# ---------------------------------------------------------------------------
class TestCorneringAxis:
    def test_zero_tilt_is_vertical(self):
        d = _bc_dict(pitch=0.0, roll=0.0)
        axis = utilities.corneringAxis(d)
        assert np.allclose(axis, [0.0, 0.0, 1.0])

    def test_is_unit_vector(self):
        d = _bc_dict(pitch=4.0, roll=3.0)
        axis = utilities.corneringAxis(d)
        assert np.linalg.norm(axis) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# corneringFrame
# ---------------------------------------------------------------------------
def _corner_dict(inlet_mag=30.0, radius=100.0, direction="left",
                 refcor=(0.0, 0.0, 0.0), center=("default",), factor=3.0,
                 domain_size=(60.0, 40.0, 20.0)):
    return {
        "BC_SETUP": {
            "REFCOR": [str(refcor[0]), str(refcor[1]), str(refcor[2])],
            "DOMAIN_PITCH": ["0.0"],
            "DOMAIN_ROLL": ["0.0"],
            "INLET_MAG": [str(inlet_mag)],
            "DOMAIN_SIZE": [str(domain_size[0]), str(domain_size[1]), str(domain_size[2])],
        },
        "CORNERING_SETUP": {
            "CORNER_RADIUS": [str(radius)],
            "CORNER_DIR": [direction],
            "CORNER_CENTER": list(center),
            "CORNER_CLEARANCE_FACTOR": [str(factor)],
        },
    }


class TestCorneringFrame:
    def test_left_turn_centre_and_sign(self):
        d = _corner_dict(inlet_mag=30.0, radius=100.0, direction="left", refcor=(0.0, 0.0, 0.0))
        omega, axis, centre = utilities.corneringFrame(d)
        assert omega == pytest.approx(30.0 / 100.0)  # positive
        assert centre[1] == pytest.approx(-100.0)     # centre to the left (-y)

    def test_right_turn_centre_and_sign(self):
        d = _corner_dict(inlet_mag=30.0, radius=80.0, direction="right", refcor=(0.0, 0.0, 0.0))
        omega, axis, centre = utilities.corneringFrame(d)
        assert omega == pytest.approx(-30.0 / 80.0)   # negative
        assert centre[1] == pytest.approx(80.0)        # centre to the right (+y)

    def test_explicit_centre_override(self):
        d = _corner_dict(center=("1.0", "2.0", "3.0"))
        _, _, centre = utilities.corneringFrame(d)
        assert np.allclose(centre, [1.0, 2.0, 3.0])

    def test_invalid_direction_exits(self):
        d = _corner_dict(direction="sideways")
        with pytest.raises(SystemExit):
            utilities.corneringFrame(d)

    def test_nonpositive_radius_exits(self):
        d = _corner_dict(radius=0.0)
        with pytest.raises(SystemExit):
            utilities.corneringFrame(d)


# ---------------------------------------------------------------------------
# checkCorneringDomain
# ---------------------------------------------------------------------------
class TestCheckCorneringDomain:
    def test_returns_min_radius_and_passes_when_large_enough(self):
        # domain y = 40 -> halfWidth 20, factor 3 -> rMin 60; radius 100 passes
        d = _corner_dict(radius=100.0, factor=3.0, domain_size=(60.0, 40.0, 20.0))
        rMin = utilities.checkCorneringDomain(d)
        assert rMin == pytest.approx(60.0)

    def test_too_small_radius_exits(self):
        d = _corner_dict(radius=50.0, factor=3.0, domain_size=(60.0, 40.0, 20.0))
        with pytest.raises(SystemExit):
            utilities.checkCorneringDomain(d)

    def test_empty_radius_exits_cleanly(self):
        d = _corner_dict()
        d["CORNERING_SETUP"]["CORNER_RADIUS"] = [""]
        with pytest.raises(SystemExit):
            utilities.checkCorneringDomain(d)
