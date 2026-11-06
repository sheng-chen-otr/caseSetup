"""
Copyright (c) 2026 Zhi Sheng Chen

Licensed under the PolyForm Noncommercial License 1.0.0 (the "License").
You may not use this file except in compliance with the License. Any
commercial use requires a separate commercial license from the copyright
holder.

You may obtain a copy of the License at:
    https://polyformproject.org/licenses/noncommercial/1.0.0

Or see the LICENSE file in the root of this repository.

This software is provided "as is", without warranty of any kind, express or
implied. See the License for the specific language governing permissions
and limitations.
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np


ArrayLike = List[float]


@dataclass
class RockerConfig:
    pivot: np.ndarray
    axis: np.ndarray
    pushrod_joint_ref: np.ndarray
    damper_joint_ref: np.ndarray
    damper_chassis: np.ndarray


@dataclass
class KinematicState:
    rvec: np.ndarray
    wheel_center: np.ndarray
    rocker_theta: float
    rack_travel: float = 0.0


def _as_vec3(v: ArrayLike, name: str) -> np.ndarray:
    arr = np.array(v, dtype=np.float64)
    if arr.shape != (3,):
        raise ValueError(f"{name} must be a 3-vector")
    return arr


def _normalize(v: np.ndarray, name: str) -> np.ndarray:
    n = np.linalg.norm(v)
    if n == 0.0:
        raise ValueError(f"{name} cannot be zero vector")
    return v / n


def _rodrigues(rvec: np.ndarray) -> np.ndarray:
    theta = np.linalg.norm(rvec)
    if theta < 1e-15:
        return np.eye(3)

    k = rvec / theta
    kx, ky, kz = k
    K = np.array(
        [[0.0, -kz, ky], [kz, 0.0, -kx], [-ky, kx, 0.0]],
        dtype=np.float64,
    )
    return np.eye(3) + math.sin(theta) * K + (1.0 - math.cos(theta)) * (K @ K)


def _rotate_about_axis(point: np.ndarray, pivot: np.ndarray, axis: np.ndarray, theta: float) -> np.ndarray:
    axis = _normalize(axis, "rocker axis")
    v = point - pivot
    ct = math.cos(theta)
    st = math.sin(theta)
    return pivot + v * ct + np.cross(axis, v) * st + axis * np.dot(axis, v) * (1.0 - ct)


def _distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


# solver result is a failure when the leftover constraint error is above this (m)
RESIDUAL_TOL = 1e-6

# names of the 8 solver constraints, in residual order, for error messages
CONSTRAINT_NAMES = (
    "UCA front leg length",
    "UCA rear leg length",
    "LCA front leg length",
    "LCA rear leg length",
    "tie rod length",
    "pushrod length",
    "wheel centre height",
    "steer / rack travel",
)

# pushrod mount detection thresholds (m)
_MOUNT_ON_ARM_TOL = 0.005      # closer than this to an arm triangle -> mounted on it
_MOUNT_AMBIGUOUS_RATIO = 0.2   # distances within 20% of each other -> can't tell


def _closest_point_on_triangle(p: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    # closest point to p on triangle abc (Ericson, Real-Time Collision Detection 5.1.5)
    ab = b - a
    ac = c - a
    ap = p - a
    d1 = np.dot(ab, ap)
    d2 = np.dot(ac, ap)
    if d1 <= 0.0 and d2 <= 0.0:
        return a
    bp = p - b
    d3 = np.dot(ab, bp)
    d4 = np.dot(ac, bp)
    if d3 >= 0.0 and d4 <= d3:
        return b
    vc = d1 * d4 - d3 * d2
    if vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0:
        return a + (d1 / (d1 - d3)) * ab
    cp = p - c
    d5 = np.dot(ab, cp)
    d6 = np.dot(ac, cp)
    if d6 >= 0.0 and d5 <= d6:
        return c
    vb = d5 * d2 - d1 * d6
    if vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0:
        return a + (d2 / (d2 - d6)) * ac
    va = d3 * d6 - d5 * d4
    if va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0:
        return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b)
    denom = 1.0 / (va + vb + vc)
    return a + ab * (vb * denom) + ac * (vc * denom)


def resolve_pushrod_mount(hardpoints: Dict[str, object]) -> str:
    """
    Which wishbone ('uca' or 'lca') the pushrod outer end is bolted to, from
    PUSHROD_OUTER_STATIC. Each wishbone is the triangle of its two inner mounts and outer
    joint; the nearest triangle wins. Raises ValueError when the point is about equally far
    from both (not a real mounting, usually a coordinate typo). Hardpoints in metres.
    """
    p = _as_vec3(hardpoints["pushrod_outer_static"], "pushrod_outer_static")
    dist = {}
    for arm in ("uca", "lca"):
        a = _as_vec3(hardpoints["%s_f_inner" % arm], "%s_f_inner" % arm)
        b = _as_vec3(hardpoints["%s_r_inner" % arm], "%s_r_inner" % arm)
        c = _as_vec3(hardpoints["%s_outer_static" % arm], "%s_outer_static" % arm)
        dist[arm] = float(np.linalg.norm(p - _closest_point_on_triangle(p, a, b, c)))

    nearest = min(dist, key=dist.get)
    other = "lca" if nearest == "uca" else "uca"
    if dist[nearest] < _MOUNT_ON_ARM_TOL:
        return nearest
    if dist[other] - dist[nearest] <= _MOUNT_AMBIGUOUS_RATIO * dist[other]:
        raise ValueError(
            "PUSHROD_OUTER_STATIC (%.1f, %.1f, %.1f) mm is %.1f mm from the UCA and %.1f mm from the LCA, "
            "can't tell which wishbone it is mounted on. Check the coordinates."
            % (p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0, dist["uca"] * 1000.0, dist["lca"] * 1000.0)
        )
    return nearest


def _signed_angle_about_axis(axis_point: np.ndarray, axis_unit: np.ndarray,
                             p_from: np.ndarray, p_to: np.ndarray) -> float:
    # angle about the hinge line that swings p_from onto p_to (perpendicular parts only)
    a = p_from - axis_point
    b = p_to - axis_point
    a = a - np.dot(a, axis_unit) * axis_unit
    b = b - np.dot(b, axis_unit) * axis_unit
    na = np.linalg.norm(a)
    nb = np.linalg.norm(b)
    if na <= 1e-15 or nb <= 1e-15:
        return 0.0
    a = a / na
    b = b / nb
    return math.atan2(float(np.dot(np.cross(a, b), axis_unit)), float(np.clip(np.dot(a, b), -1.0, 1.0)))


class DoubleWishbonePushrodSolver:
    """
    solve wheel/upright kinematics for a double wishbone + pushrod/rocker corner.
    convention: x forward, y lateral, z vertical. the upright is a rigid body with
    unknown pose (rvec + wheel-center translation), rocker angle solved together to
    hold pushrod length.

    required hardpoint keys: wheel_center_static, uca/lca/tie inner+outer, pushrod_outer_static,
    rocker {pivot, axis, pushrod_joint_ref, damper_joint_ref, damper_chassis}.
    optional: wheel_axis_local ([0,1,0]), wheel_forward_local ([1,0,0]),
    pushrod_mount ('uca'/'lca', detected from pushrod_outer_static when missing).
    the pushrod outer end is bolted to a wishbone, so it swings with that arm about the
    arm's inner-mount axis. the vehicle centreline is y = 0.
    """

    def __init__(self, hardpoints: Dict[str, object]):
        self.hp = hardpoints
        self._load_geometry()
        self._build_static_lengths()

    def _load_geometry(self) -> None:
        wc0 = _as_vec3(self.hp["wheel_center_static"], "wheel_center_static")
        self.wc0 = wc0

        self.inner = {
            "uca_f": _as_vec3(self.hp["uca_f_inner"], "uca_f_inner"),
            "uca_r": _as_vec3(self.hp["uca_r_inner"], "uca_r_inner"),
            "lca_f": _as_vec3(self.hp["lca_f_inner"], "lca_f_inner"),
            "lca_r": _as_vec3(self.hp["lca_r_inner"], "lca_r_inner"),
            "tie": _as_vec3(self.hp["tie_inner"], "tie_inner"),
        }

        # points rigid on the upright
        self.outer0_global = {
            "uca": _as_vec3(self.hp["uca_outer_static"], "uca_outer_static"),
            "lca": _as_vec3(self.hp["lca_outer_static"], "lca_outer_static"),
            "tie": _as_vec3(self.hp["tie_outer_static"], "tie_outer_static"),
            "wheel_center": wc0,
        }

        # upright local frame starts aligned to global at the static pose
        self.outer_local = {k: (v - wc0) for k, v in self.outer0_global.items()}

        # pushrod outer end rides on a wishbone, not the upright
        self.pushrod_outer_static = _as_vec3(self.hp["pushrod_outer_static"], "pushrod_outer_static")
        mount = self.hp.get("pushrod_mount") or resolve_pushrod_mount(self.hp)
        mount = str(mount).strip().lower()
        if mount not in ("uca", "lca"):
            raise ValueError("pushrod_mount must be 'uca' or 'lca', got %r" % (mount,))
        self.pushrod_mount = mount
        self.mount_hinge_point = self.inner["%s_f" % mount]
        self.mount_hinge_axis = _normalize(self.inner["%s_r" % mount] - self.inner["%s_f" % mount],
                                           "%s hinge axis" % mount)

        # inboard direction: towards the centreline y = 0
        if abs(wc0[1]) < 1e-9:
            raise ValueError("wheel_center_static y is 0, can't tell which side of the car the wheel is on")
        self.inboard_y = -1.0 if wc0[1] > 0.0 else 1.0

        rocker_dict = self.hp["rocker"]
        self.rocker = RockerConfig(
            pivot=_as_vec3(rocker_dict["pivot"], "rocker.pivot"),
            axis=_as_vec3(rocker_dict["axis"], "rocker.axis"),
            pushrod_joint_ref=_as_vec3(rocker_dict["pushrod_joint_ref"], "rocker.pushrod_joint_ref"),
            damper_joint_ref=_as_vec3(rocker_dict["damper_joint_ref"], "rocker.damper_joint_ref"),
            damper_chassis=_as_vec3(rocker_dict["damper_chassis"], "rocker.damper_chassis"),
        )

        self.wheel_axis_local = _normalize(
            _as_vec3(self.hp.get("wheel_axis_local", [0.0, 1.0, 0.0]), "wheel_axis_local"),
            "wheel_axis_local",
        )
        self.wheel_forward_local = _normalize(
            _as_vec3(self.hp.get("wheel_forward_local", [1.0, 0.0, 0.0]), "wheel_forward_local"),
            "wheel_forward_local",
        )

        # tie-rod inner (rack end) static pos + rack direction. steering slides it along rack_axis
        self.tie_inner_static = self.inner["tie"]
        self.rack_axis = _normalize(
            _as_vec3(self.hp.get("rack_axis", [0.0, 1.0, 0.0]), "rack_axis"),
            "rack_axis",
        )

    def _build_static_lengths(self) -> None:
        o = self.outer0_global
        i = self.inner

        self.lengths = {
            "uca_f": _distance(o["uca"], i["uca_f"]),
            "uca_r": _distance(o["uca"], i["uca_r"]),
            "lca_f": _distance(o["lca"], i["lca_f"]),
            "lca_r": _distance(o["lca"], i["lca_r"]),
            "tie": _distance(o["tie"], i["tie"]),
            "pushrod": _distance(self.pushrod_outer_static, self.rocker.pushrod_joint_ref),
            "damper": _distance(self.rocker.damper_joint_ref, self.rocker.damper_chassis),
        }

    def _rocker_points(self, theta: float) -> Dict[str, np.ndarray]:
        push = _rotate_about_axis(
            self.rocker.pushrod_joint_ref,
            self.rocker.pivot,
            self.rocker.axis,
            theta,
        )
        damper = _rotate_about_axis(
            self.rocker.damper_joint_ref,
            self.rocker.pivot,
            self.rocker.axis,
            theta,
        )
        return {"pushrod": push, "damper": damper}

    def _upright_global_points(self, rvec: np.ndarray, wheel_center: np.ndarray) -> Dict[str, np.ndarray]:
        R = _rodrigues(rvec)
        return {k: (R @ v_local) + wheel_center for k, v_local in self.outer_local.items()}

    def _pushrod_outer(self, pg: Dict[str, np.ndarray]) -> np.ndarray:
        # the mount arm's angle is set by where its outer ball joint is; swing the static
        # pushrod point about the arm's inner-mount axis by the same angle
        arm = self.pushrod_mount
        theta = _signed_angle_about_axis(self.mount_hinge_point, self.mount_hinge_axis,
                                         self.outer0_global[arm], pg[arm])
        return _rotate_about_axis(self.pushrod_outer_static, self.mount_hinge_point,
                                  self.mount_hinge_axis, theta)

    def _toe_rad_from_rvec(self, rvec: np.ndarray) -> float:
        R = _rodrigues(rvec)
        fwd = R @ self.wheel_forward_local
        return math.atan2(fwd[1], fwd[0] + 1e-15)

    def _residual(
        self,
        x: np.ndarray,
        target_wheel_z: float,
        steer_target: float = 0.0,
        steer_mode: str = "none",
    ) -> np.ndarray:
        rvec = x[0:3]
        wc = x[3:6]
        theta = float(x[6])
        d = float(x[7])

        tie_inner = self.tie_inner_static + d * self.rack_axis

        pg = self._upright_global_points(rvec, wc)
        rk = self._rocker_points(theta)

        res = np.zeros(8, dtype=np.float64)
        res[0] = _distance(pg["uca"], self.inner["uca_f"]) - self.lengths["uca_f"]
        res[1] = _distance(pg["uca"], self.inner["uca_r"]) - self.lengths["uca_r"]
        res[2] = _distance(pg["lca"], self.inner["lca_f"]) - self.lengths["lca_f"]
        res[3] = _distance(pg["lca"], self.inner["lca_r"]) - self.lengths["lca_r"]
        res[4] = _distance(pg["tie"], tie_inner) - self.lengths["tie"]
        res[5] = _distance(self._pushrod_outer(pg), rk["pushrod"]) - self.lengths["pushrod"]
        res[6] = wc[2] - target_wheel_z

        # 8th constraint closes the rack DOF: prescribe rack travel, prescribe toe angle, or lock at zero
        if steer_mode == "rack":
            res[7] = d - steer_target
        elif steer_mode == "angle":
            toe = self._toe_rad_from_rvec(rvec)
            res[7] = self.lengths["tie"] * (toe - math.radians(steer_target))
        else:
            res[7] = d
        return res

    def _numeric_jacobian(
        self,
        x: np.ndarray,
        target_wheel_z: float,
        steer_target: float = 0.0,
        steer_mode: str = "none",
        eps: float = 1e-6,
    ) -> np.ndarray:
        r0 = self._residual(x, target_wheel_z, steer_target, steer_mode)
        J = np.zeros((r0.size, x.size), dtype=np.float64)
        for i in range(x.size):
            xp = x.copy()
            xm = x.copy()
            xp[i] += eps
            xm[i] -= eps
            rp = self._residual(xp, target_wheel_z, steer_target, steer_mode)
            rm = self._residual(xm, target_wheel_z, steer_target, steer_mode)
            J[:, i] = (rp - rm) / (2.0 * eps)
        return J

    def solve_for_wheel_travel(
        self,
        wheel_dz: float,
        steer: float = 0.0,
        steer_mode: str = "none",
        initial_state: Optional[KinematicState] = None,
        max_iter: int = 80,
        tol: float = 1e-9,
        damping: float = 1e-9,
    ) -> Dict[str, object]:
        """
        solve kinematics for a target wheel-center vertical travel, optionally with steer.
        steer_mode: "none" (rack locked), "angle" (steer = target toe in deg, solve rack),
        "rack" (steer = rack travel directly). initial_state warm-starts a sweep.
        returns solved points, rocker angle, damper length, rack travel and alignment metrics.
        """
        target_z = self.wc0[2] + wheel_dz

        if initial_state is None:
            x = np.zeros(8, dtype=np.float64)
            x[3:6] = self.wc0 + np.array([0.0, 0.0, wheel_dz], dtype=np.float64)
            x[6] = 0.0
            x[7] = steer if steer_mode == "rack" else 0.0
        else:
            x = np.zeros(8, dtype=np.float64)
            x[0:3] = initial_state.rvec
            x[3:6] = initial_state.wheel_center
            x[6] = initial_state.rocker_theta
            x[7] = getattr(initial_state, "rack_travel", 0.0)

        for _ in range(max_iter):
            r = self._residual(x, target_z, steer, steer_mode)
            err = float(np.linalg.norm(r))
            if err < tol:
                break

            J = self._numeric_jacobian(x, target_z, steer, steer_mode)
            JTJ = J.T @ J + damping * np.eye(x.size)
            step = np.linalg.solve(JTJ, -J.T @ r)

            # backtracking line-search for robustness near singular positions
            alpha = 1.0
            accepted = False
            for _ls in range(10):
                x_try = x + alpha * step
                # clamp rocker angle to ±90°
                x_try[6] = np.clip(x_try[6], -np.pi/2, np.pi/2)
                r_try = self._residual(x_try, target_z, steer, steer_mode)
                if np.linalg.norm(r_try) < err:
                    x = x_try
                    accepted = True
                    break
                alpha *= 0.5

            if not accepted:
                x = x + 0.2 * step
                # clamp rocker angle after the fallback step
                x[6] = np.clip(x[6], -np.pi/2, np.pi/2)

        rvec = x[0:3]
        wheel_center = x[3:6]
        theta = float(x[6])
        rack_travel = float(x[7])

        pg = self._upright_global_points(rvec, wheel_center)
        pg["pushrod"] = self._pushrod_outer(pg)
        rk = self._rocker_points(theta)
        damper_len = _distance(rk["damper"], self.rocker.damper_chassis)

        R = _rodrigues(rvec)
        wheel_axis_global = _normalize(R @ self.wheel_axis_local, "wheel_axis_global")
        wheel_forward_global = _normalize(R @ self.wheel_forward_local, "wheel_forward_global")

        # camber relative to the car body: spin axis flipped to point inboard so the sign does not
        # depend on how WHEEL_AXIS_LOCAL was entered. negative = top of the wheel leaning inboard.
        # asin of the vertical part keeps it independent of steer
        axis_inboard = wheel_axis_global if wheel_axis_global[1] * self.inboard_y >= 0.0 else -wheel_axis_global
        camber_deg = math.degrees(math.asin(float(np.clip(axis_inboard[2], -1.0, 1.0))))
        toe_deg = math.degrees(math.atan2(wheel_forward_global[1], wheel_forward_global[0] + 1e-15))

        residual = self._residual(x, target_z, steer, steer_mode)

        return {
            "state": KinematicState(
                rvec=rvec.copy(),
                wheel_center=wheel_center.copy(),
                rocker_theta=theta,
                rack_travel=rack_travel,
            ),
            "wheel_center": wheel_center,
            "upright_points": pg,
            "rocker_points": rk,
            "rocker_theta_rad": theta,
            "rocker_theta_deg": math.degrees(theta),
            "damper_length": damper_len,
            "damper_delta": damper_len - self.lengths["damper"],
            "camber_deg": camber_deg,
            "toe_deg": toe_deg,
            "rack_travel": rack_travel,
            "residual_norm": float(np.linalg.norm(residual)),
            "residual": residual,
            "rocker_at_limit": abs(theta) >= (math.pi / 2.0 - 1e-6),
            "pushrod_mount": self.pushrod_mount,
        }

    @staticmethod
    def describe_failure(solved: Dict[str, object], wheel_dz: float, steer: float, steer_mode: str) -> str:
        """readable explanation of a solve whose residual is above RESIDUAL_TOL (lengths in m)."""
        residual = np.asarray(solved["residual"], dtype=np.float64)
        order = np.argsort(-np.abs(residual))
        lines = []
        for i in order[:3]:
            if abs(residual[i]) <= RESIDUAL_TOL:
                break
            lines.append("%s off by %.3f mm" % (CONSTRAINT_NAMES[i], residual[i] * 1000.0))
        if steer_mode == "angle":
            steerText = "steer %.3f deg" % (steer,)
        elif steer_mode == "rack":
            steerText = "rack travel %.3f mm" % (steer * 1000.0,)
        else:
            steerText = "no steer"
        msg = ("suspension solve did not converge for wheel travel %.3f mm, %s (leftover error %.4f mm).\n"
               "\tLargest constraint errors: %s"
               % (wheel_dz * 1000.0, steerText, solved["residual_norm"] * 1000.0, "; ".join(lines)))
        if solved.get("rocker_at_limit"):
            msg += "\n\tThe rocker hit its +/-90 deg limit."
        msg += ("\n\tLikely causes: the requested travel/steer is beyond what the linkage can reach, a "
                "hardpoint typo, or wrong units (check HARDPOINT_SCALE).")
        return msg

    def sweep_wheel_travel(self, wheel_dz_values: List[float]) -> List[Dict[str, object]]:
        out = []
        state = None
        for dz in wheel_dz_values:
            solved = self.solve_for_wheel_travel(dz, initial_state=state)
            out.append(solved)
            state = solved["state"]
        return out


if __name__ == "__main__":
    # Minimal example with synthetic coordinates.
    sample_hardpoints = {
        "wheel_center_static": [1.45, 0.78, 0.31],
        "uca_f_inner": [1.25, 0.45, 0.50],
        "uca_r_inner": [1.10, 0.46, 0.50],
        "lca_f_inner": [1.25, 0.40, 0.15],
        "lca_r_inner": [1.10, 0.41, 0.15],
        "tie_inner": [1.05, 0.52, 0.26],
        "uca_outer_static": [1.42, 0.77, 0.47],
        "lca_outer_static": [1.43, 0.79, 0.20],
        "tie_outer_static": [1.40, 0.76, 0.27],
        "pushrod_outer_static": [1.43, 0.79, 0.20],  # on the LCA outer joint
        "rocker": {
            "pivot": [1.00, 0.35, 0.55],
            "axis": [1.0, 0.0, 0.0],
            "pushrod_joint_ref": [1.03, 0.44, 0.50],
            "damper_joint_ref": [0.98, 0.29, 0.50],
            "damper_chassis": [0.95, 0.20, 0.62],
        },
    }

    solver = DoubleWishbonePushrodSolver(sample_hardpoints)
    sweep = solver.sweep_wheel_travel([0.0, -0.01, -0.02, 0.01])
    for row in sweep:
        wc = row["wheel_center"]
        print(
            f"z={wc[2]:.6f}, rocker_deg={row['rocker_theta_deg']:.3f}, "
            f"damper_delta={row['damper_delta']:.6f}, camber={row['camber_deg']:.3f}, toe={row['toe_deg']:.3f}, "
            f"res={row['residual_norm']:.3e}"
        )
