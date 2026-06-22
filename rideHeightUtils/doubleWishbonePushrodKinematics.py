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


class DoubleWishbonePushrodSolver:
    """
    Solve wheel/upright kinematics for a double wishbone + pushrod/rocker suspension.

    Coordinate convention used by default:
    - x: longitudinal (vehicle forward)
    - y: lateral
    - z: vertical

    The upright is treated as a rigid body with unknown pose (rotation vector + wheel center translation).
    Rocker angle is solved simultaneously to enforce pushrod length.

    Required hardpoint dictionary keys:
    - wheel_center_static
    - uca_f_inner, uca_r_inner, lca_f_inner, lca_r_inner, tie_inner
    - uca_outer_static, lca_outer_static, tie_outer_static, pushrod_outer_static
    - rocker: {
        pivot, axis, pushrod_joint_ref, damper_joint_ref, damper_chassis
      }

    Optional keys:
    - wheel_axis_local (default [0, 1, 0])
    - wheel_forward_local (default [1, 0, 0])
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

        self.outer0_global = {
            "uca": _as_vec3(self.hp["uca_outer_static"], "uca_outer_static"),
            "lca": _as_vec3(self.hp["lca_outer_static"], "lca_outer_static"),
            "tie": _as_vec3(self.hp["tie_outer_static"], "tie_outer_static"),
            "pushrod": _as_vec3(self.hp["pushrod_outer_static"], "pushrod_outer_static"),
            "wheel_center": wc0,
        }

        # Upright local frame is initialized aligned to global at static pose.
        self.outer_local = {k: (v - wc0) for k, v in self.outer0_global.items()}

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

        # Tie-rod inner (rack end) static position and the rack travel direction.
        # Steering displaces the inner joint along rack_axis by the rack-travel DOF.
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
            "pushrod": _distance(o["pushrod"], self.rocker.pushrod_joint_ref),
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
        res[5] = _distance(pg["pushrod"], rk["pushrod"]) - self.lengths["pushrod"]
        res[6] = wc[2] - target_wheel_z

        # 8th constraint closes the rack DOF: either prescribe rack travel directly,
        # prescribe the road-wheel (toe) angle, or lock the rack at zero (passive).
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
        Solve kinematics for target wheel center vertical travel, optionally with steer.

        :param wheel_dz: Desired wheel center z displacement from static (same units as hardpoints)
        :param steer: Steer command. Interpreted per steer_mode.
        :param steer_mode: One of:
            - "none"  : rack locked at zero (passive ride-height only).
            - "angle" : ``steer`` is the target road-wheel (toe) angle in degrees; the
                        solver finds the rack travel that achieves it.
            - "rack"  : ``steer`` is the rack travel directly (same length units as
                        hardpoints), e.g. shared across an axle for Ackermann.
        :param initial_state: Optional warm-start state; use previous solution for sweep robustness
        :return: Dict with solved points, rocker angle, damper length, rack travel, and alignment metrics
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

            # Backtracking line-search to improve robustness around singular positions.
            alpha = 1.0
            accepted = False
            for _ls in range(10):
                x_try = x + alpha * step
                # Clamp rocker angle to physically reasonable bounds (±90°)
                x_try[6] = np.clip(x_try[6], -np.pi/2, np.pi/2)
                r_try = self._residual(x_try, target_z, steer, steer_mode)
                if np.linalg.norm(r_try) < err:
                    x = x_try
                    accepted = True
                    break
                alpha *= 0.5

            if not accepted:
                x = x + 0.2 * step
                # Also clamp rocker angle after fallback update
                x[6] = np.clip(x[6], -np.pi/2, np.pi/2)

        rvec = x[0:3]
        wheel_center = x[3:6]
        theta = float(x[6])
        rack_travel = float(x[7])

        pg = self._upright_global_points(rvec, wheel_center)
        rk = self._rocker_points(theta)
        damper_len = _distance(rk["damper"], self.rocker.damper_chassis)

        R = _rodrigues(rvec)
        wheel_axis_global = _normalize(R @ self.wheel_axis_local, "wheel_axis_global")
        wheel_forward_global = _normalize(R @ self.wheel_forward_local, "wheel_forward_global")

        camber_deg = math.degrees(math.atan2(wheel_axis_global[2], abs(wheel_axis_global[1]) + 1e-15))
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
        }

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
        "pushrod_outer_static": [1.39, 0.75, 0.25],
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
