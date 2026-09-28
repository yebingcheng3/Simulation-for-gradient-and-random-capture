# -*- coding: utf-8 -*-
"""
Single-window 3D KT-unit microtubule capture simulator

Final hybrid design from discussion:
- GUI/animation style from the user's previous simulator, but one animation window only.
- Calculation/export rules follow the BOSS conservative event-driven model.
- KT objects are single spherical KT units, not sister-pair objects.
- KT positions are generated uniformly at random within the nuclear ellipsoid.
- Exactly 46 KT units are generated: 20 distal + 26 proximal.
- Distal/proximal classification uses the BOSS center-point rule.
- First burst: 23 MTs from each centrosome at t = 0.
- Post-burst emission rate is per centrosome, adjustable from 0.1000 to 2.0000 MT/s.
- Random MT directions: tangent cone + without-replacement angular occupancy.
- Directed MTs: each MT independently has probability f to be directed.
- Random, Directed, and Pivoting share one chronological growth/capture engine;
  conditions differ through their initial direction bias and D_MT updates.
- Pivoting uses D_MT angular diffusion: 0, 3.16e-6, 1e-5, 3.16e-5, 1e-4,
  3.16e-4 rad^2/s.
- Pivoting occurs only during MT growth; nonzero direction changes permit swept-shaft capture.
- One MT can capture multiple KT units along its path.
- Every analyzed MT exports one mean deflection row: captured MTs stop at first
  engagement, while uncaptured MTs stop at true closest approach. Sampling is 0.25 s.

Run GUI:
    python random_directed_microtubule_kinetochore_capture_3D_simulation_unified_engine.py

Run direct BOSS-style output:
    python random_directed_microtubule_kinetochore_capture_3D_simulation_unified_engine.py --direct

Run an exact six-condition D_MT pivoting sweep:
    python random_directed_microtubule_kinetochore_capture_3D_simulation_unified_engine.py --direct --pivoting-sweep dmt

Dependencies:
    pip install PySide6
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import math
import os
import random
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Tuple

try:
    from PySide6.QtCore import Qt, QTimer, QPointF, QPoint, QSize, QRect
    from PySide6.QtGui import QColor, QFont, QPainter, QPen, QBrush, QPainterPath, QImage
    from PySide6.QtSvg import QSvgGenerator
    from PySide6.QtWidgets import (
        QApplication, QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
        QFormLayout, QFrame, QGroupBox, QHBoxLayout, QLabel, QMessageBox,
        QProgressDialog, QPushButton, QRadioButton, QSpinBox, QVBoxLayout, QWidget,
        QTableWidget, QTableWidgetItem, QHeaderView, QSlider,
    )
    QT_AVAILABLE = True
except ModuleNotFoundError:
    QT_AVAILABLE = False

    class _DummyMeta(type):
        def __getattr__(cls, name):
            return cls()

    class _Dummy(metaclass=_DummyMeta):
        def __init__(self, *args, **kwargs):
            pass
        def __call__(self, *args, **kwargs):
            return self
        def __getattr__(self, name):
            return self
        def __bool__(self):
            return False

    Qt = _Dummy()
    QTimer = QPointF = QPoint = QSize = QRect = QColor = QFont = QPainter = QPen = QBrush = QPainterPath = QImage = QSvgGenerator = _Dummy
    QApplication = QButtonGroup = QCheckBox = QComboBox = QDoubleSpinBox = QFileDialog = _Dummy
    QFormLayout = QFrame = QGroupBox = QHBoxLayout = QLabel = QMessageBox = _Dummy
    QProgressDialog = QPushButton = QRadioButton = QSpinBox = QVBoxLayout = QWidget = QSlider = _Dummy

Vec3 = Tuple[float, float, float]


def kt_sphere_volume_um3(diameter_um: float) -> float:
    """Volume of the spherical KT unit used by the simulator."""
    d = max(0.0, float(diameter_um))
    return math.pi / 6.0 * d ** 3


def effective_kt_diameter_um(base_diameter_um: float, volume_multiplier: float) -> float:
    """Convert base diameter and volume multiplier into an effective diameter.

    The simulator collision model uses a spherical KT unit radius. Scaling the
    KT capture volume by M therefore scales the diameter by M^(1/3), not by M.
    """
    base = max(0.001, float(base_diameter_um))
    mult = max(1e-9, float(volume_multiplier))
    return base * (mult ** (1.0 / 3.0))


# -----------------------------------------------------------------------------
# Geometry defaults, um
# -----------------------------------------------------------------------------
NUC_LONG_AXIS = 17.822970831590442
NUC_INTERMEDIATE_AXIS = 13.22539119077223
NUC_SHORT_AXIS = 9.167840990422043
NUC_A = NUC_LONG_AXIS / 2.0
NUC_B = NUC_INTERMEDIATE_AXIS / 2.0
NUC_C = NUC_SHORT_AXIS / 2.0
NUC_RADII = (NUC_A, NUC_B, NUC_C)

CENTROSOME_DISTANCE = 9.402413793
CS_A = (0.0, 0.0, -CENTROSOME_DISTANCE / 2.0)
CS_B = (0.0, 0.0, +CENTROSOME_DISTANCE / 2.0)
CENTROSOMES = (CS_A, CS_B)
D_HALF = CENTROSOME_DISTANCE / 2.0

TOTAL_KT_UNITS = 46
DISTAL_KT_UNITS = 20
PROXIMAL_KT_UNITS = 26
FIRST_BURST_PER_CENTROSOME = 23

# Visual typography
DISPLAY_FONT_FAMILY = "Arial"
VISUAL_FONT_SCALE = 2.0

def visual_font_size(base_pt: float) -> int:
    """Scale publication-view typography without changing simulation geometry."""
    return max(1, int(round(float(base_pt) * VISUAL_FONT_SCALE)))

# Visual colors
COL_BG = QColor(3, 9, 12)
COL_NUC_WIRE = QColor(110, 80, 150, 78)
COL_MT = QColor(0, 170, 82, 155)
COL_MT_DONE = QColor(0, 170, 82, 160)
COL_MT_FAINT = QColor(0, 170, 82, 45)
COL_CS = QColor(0, 175, 85, 220)
COL_KT = QColor(255, 0, 240, 230)
COL_KT_PROX = QColor(255, 0, 240, 145)
COL_KT_DONE = QColor(175, 190, 188, 85)
COL_KT_DONE_LEGEND = QColor(175, 190, 188, 170)
COL_TEXT = QColor(215, 230, 238)

CONDITIONS: List[Tuple[str, float]] = [
    ("Random", 0.0),
    ("20% directed", 0.20),
    ("40% directed", 0.40),
    ("60% directed", 0.60),
    ("80% directed", 0.80),
    ("100% directed", 1.00),
]

# Export-only condition list.  The animation radio buttons intentionally keep
# the original coarse set above, while batch export can use 10% increments.
EXPORT_CONDITIONS: List[Tuple[str, float]] = [
    ("Random", 0.0),
    ("10% directed", 0.10),
    ("20% directed", 0.20),
    ("30% directed", 0.30),
    ("40% directed", 0.40),
    ("50% directed", 0.50),
    ("60% directed", 0.60),
    ("70% directed", 0.70),
    ("80% directed", 0.80),
    ("90% directed", 0.90),
    ("100% directed", 1.00),
]

# Pivoting sweeps are deliberately separate from the unchanged
# Random/Directed conditions above.  Zero is included as one of the six
# conditions in each sweep; this is a user-selected graded/log-style range
# spanning human-to-yeast-relevant values, not six direct human measurements.
# Zero remains a pivoting-mode condition for metadata, but has no swept advantage.
PIVOTING_DMT_CONDITIONS: List[Tuple[str, float]] = [
    ("Pivoting D_MT = 0 rad^2/s", 0.0),
    ("Pivoting D_MT = 3.16e-6 rad^2/s", 3.16e-6),
    ("Pivoting D_MT = 1e-5 rad^2/s", 1e-5),
    ("Pivoting D_MT = 3.16e-5 rad^2/s", 3.16e-5),
    ("Pivoting D_MT = 1e-4 rad^2/s", 1e-4),
    ("Pivoting D_MT = 3.16e-4 rad^2/s", 3.16e-4),
]

@dataclass(frozen=True)
class ConditionSpec:
    """One animation/export condition and its model interpretation."""

    name: str
    directed_prob: float = 0.0
    pivoting_mode: str = "none"  # none or dmt
    pivoting_value: float = 0.0
    family: str = "random_directed"

    @property
    def d_mt_rad2_s(self) -> float:
        return self.pivoting_value if self.pivoting_mode == "dmt" else 0.0

def directed_condition_specs(conditions: Sequence[Tuple[str, float]]) -> List[ConditionSpec]:
    return [ConditionSpec(name, prob) for name, prob in conditions]


ANIMATION_CONDITIONS: List[ConditionSpec] = (
    directed_condition_specs(CONDITIONS)
    + [ConditionSpec(name, 0.0, "dmt", value, "pivoting_dmt") for name, value in PIVOTING_DMT_CONDITIONS]
)

PIVOT_DT_S = 1.0
THETA_SAMPLE_DT_S = 0.25  # Independent of GUI FPS and pivot update interval.

PROGRESS_MILESTONES = [0.10, 0.25, 0.50, 0.75, 0.90, 1.00]

# Defaults for command-line direct output mode. These mirror the BOSS script.
DIRECT_OUTPUT_DIR = "capture_sim_output"
DIRECT_N_CELLS = 500
DIRECT_N_RUNS = 3
DIRECT_F_DIRECTED = 0.0

# Random seeds are now automatic.  They are still recorded in exports for QC,
# but the GUI no longer exposes a fixed seed field.  This prevents two separate
# export runs from producing identical data unless the user intentionally edits
# the code.
def new_random_seed() -> int:
    return random.SystemRandom().randint(1, 2_000_000_000)



# -----------------------------------------------------------------------------
# Vector helpers
# -----------------------------------------------------------------------------
def v_add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def v_sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def v_mul(a: Vec3, s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def v_dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def v_cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def v_len(a: Vec3) -> float:
    return math.sqrt(max(0.0, v_dot(a, a)))


def v_norm(a: Vec3) -> Vec3:
    L = v_len(a)
    if L <= 1e-12:
        return (1.0, 0.0, 0.0)
    return (a[0] / L, a[1] / L, a[2] / L)


def random_unit_vector(rng: random.Random) -> Vec3:
    z = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * math.pi)
    r = math.sqrt(max(0.0, 1.0 - z * z))
    return (r * math.cos(phi), r * math.sin(phi), z)


def orthonormal_basis(axis: Vec3) -> Tuple[Vec3, Vec3, Vec3]:
    w = v_norm(axis)
    helper = (0.0, 0.0, 1.0) if abs(w[2]) < 0.9 else (0.0, 1.0, 0.0)
    u = v_norm(v_cross(helper, w))
    v = v_norm(v_cross(w, u))
    return u, v, w


def sample_cone_direction(axis: Vec3, half_angle_rad: float, rng: random.Random) -> Vec3:
    """Uniform sample in solid angle within a cone."""
    half_angle_rad = max(0.0, min(math.pi, half_angle_rad))
    u, v, w = orthonormal_basis(axis)
    cos_min = math.cos(half_angle_rad)
    cos_t = rng.uniform(cos_min, 1.0)
    sin_t = math.sqrt(max(0.0, 1.0 - cos_t * cos_t))
    phi = rng.uniform(0.0, 2.0 * math.pi)
    return v_norm(v_add(v_add(v_mul(u, sin_t * math.cos(phi)), v_mul(v, sin_t * math.sin(phi))), v_mul(w, cos_t)))


def angular_distance(d1: Vec3, d2: Vec3) -> float:
    return math.acos(max(-1.0, min(1.0, v_dot(v_norm(d1), v_norm(d2)))))


def inside_ellipsoid(p: Vec3, radii: Vec3) -> bool:
    x, y, z = p
    a, b, c = radii
    if a <= 0 or b <= 0 or c <= 0:
        return False
    return (x * x) / (a * a) + (y * y) / (b * b) + (z * z) / (c * c) <= 1.0 + 1e-12


def random_point_in_ellipsoid(rng: random.Random, radii: Vec3) -> Vec3:
    d = random_unit_vector(rng)
    r = rng.random() ** (1.0 / 3.0)
    return (d[0] * radii[0] * r, d[1] * radii[1] * r, d[2] * radii[2] * r)


def rotate_vec(p: Vec3, yaw: float, pitch: float) -> Vec3:
    x, y, z = p
    cy, sy = math.cos(yaw), math.sin(yaw)
    x1 = cy * x - sy * y
    y1 = sy * x + cy * y
    z1 = z
    cp, sp = math.cos(pitch), math.sin(pitch)
    y2 = cp * y1 - sp * z1
    z2 = sp * y1 + cp * z1
    return (x1, y2, z2)

# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------
@dataclass
class KTUnit:
    index: int
    pos: Vec3
    is_distal: bool
    alive: bool = True
    capture_time: Optional[float] = None

    @property
    def group(self) -> str:
        return "distal" if self.is_distal else "proximal"


@dataclass
class MTRecord:
    mt_id: int
    cs_index: int
    emit_time: float
    direction: Vec3
    is_directed: bool
    target_index: Optional[int]
    max_len: float
    hits: List[Tuple[int, float]] = field(default_factory=list)  # (kt_index, contact_distance)
    hit_capture_modes: List[str] = field(default_factory=list)  # plus_end or lateral_shaft
    no_hit_deflection: Optional[float] = None
    nearest_kt_index: Optional[int] = None
    final_direction: Optional[Vec3] = None
    pivot_step_count: int = 0
    # For pivoting, this stores the per-MT trajectory mean (not a frame angle).
    mean_deflection_angle: Optional[float] = None


@dataclass
class SimulationParams:
    kt_generation_mode: str = "Uniform random"
    # Effective diameter used for collision, non-overlap placement, and drawing.
    kt_diameter_um: float = 0.80
    # User-facing base diameter and volume scale used to compute kt_diameter_um.
    kt_base_diameter_um: float = 0.80
    kt_volume_multiplier: float = 1.0
    mt_diameter_um: float = 0.025
    mt_velocity_um_s: float = 0.116
    emission_rate_per_centrosome: float = 0.15
    mt_extra_length_um: float = 1.0
    # All conditions use the same engine; "dmt" adds angular diffusion while
    # "none" keeps the initial direction fixed.
    pivoting_mode: str = "none"
    d_mt_rad2_s: float = 0.0
    pivot_dt_s: float = PIVOT_DT_S
    # Internal auto-generated seed used only within one export/animation reset.
    seed: int = 1

    @property
    def kt_radius_um(self) -> float:
        return self.kt_diameter_um / 2.0

    @property
    def kt_base_volume_um3(self) -> float:
        return kt_sphere_volume_um3(self.kt_base_diameter_um)

    @property
    def kt_effective_volume_um3(self) -> float:
        return kt_sphere_volume_um3(self.kt_diameter_um)

    @property
    def pivoting_enabled(self) -> bool:
        return self.pivoting_mode == "dmt"


def params_for_condition(params: SimulationParams, condition: ConditionSpec) -> SimulationParams:
    """Return condition-specific parameters without mutating shared base params."""
    if condition.pivoting_mode == "dmt":
        return replace(
            params,
            pivoting_mode="dmt",
            d_mt_rad2_s=max(0.0, float(condition.d_mt_rad2_s)),
            pivot_dt_s=PIVOT_DT_S,
        )
    return replace(
        params,
        pivoting_mode="none",
        d_mt_rad2_s=0.0,
        pivot_dt_s=PIVOT_DT_S,
    )


def condition_seed_offset(condition: ConditionSpec) -> int:
    """Conditions share initial RNG streams for paired comparisons."""
    return 0


def condition_metadata(condition: ConditionSpec) -> Dict[str, object]:
    """Common export metadata, including the relevant pivoting parameter."""
    if condition.pivoting_mode == "dmt":
        parameter = f"D_MT={condition.d_mt_rad2_s:.6g} rad^2/s"
    else:
        parameter = "none"
    return {
        "condition": condition.name,
        "condition_family": condition.family,
        "directed_percent": condition.directed_prob * 100.0,
        "pivoting_mode": condition.pivoting_mode,
        "pivoting_parameter": parameter,
        "d_mt_rad2_s": condition.d_mt_rad2_s,
        "pivot_dt_s": PIVOT_DT_S,
        "theta_sample_dt_s": THETA_SAMPLE_DT_S,
    }

# -----------------------------------------------------------------------------
# Uniform-random KT layout generation
# -----------------------------------------------------------------------------


def script_dir() -> str:
    try:
        return os.path.dirname(os.path.abspath(__file__))
    except Exception:
        return os.getcwd()

def is_distal_point(p: Vec3) -> bool:
    return v_len(v_sub(p, CS_A)) > D_HALF and v_len(v_sub(p, CS_B)) > D_HALF


def centers_nonoverlap(p: Vec3, existing: Sequence[KTUnit], min_dist: float) -> bool:
    return all(v_len(v_sub(p, kt.pos)) >= min_dist - 1e-12 for kt in existing)


def candidate_point(params: SimulationParams, rng: random.Random, effective_radii: Vec3) -> Optional[Vec3]:
    """Generate one candidate KT center uniformly within the accessible nuclear ellipsoid."""
    return random_point_in_ellipsoid(rng, effective_radii)


def generate_kt_layout(params: SimulationParams, layout_seed: int, max_attempts: int = 1000000) -> List[KTUnit]:
    """Generate exactly 20 distal + 26 proximal non-overlapping KT unit centers."""
    rng = random.Random(layout_seed)
    r = params.kt_radius_um
    effective_radii = (NUC_A - r, NUC_B - r, NUC_C - r)
    if min(effective_radii) <= 0:
        raise ValueError("KT unit diameter is too large for the nucleus. Please reduce KT unit diameter.")

    distal: List[KTUnit] = []
    proximal: List[KTUnit] = []

    def all_units() -> List[KTUnit]:
        return distal + proximal

    attempts = 0
    while (len(distal) < DISTAL_KT_UNITS or len(proximal) < PROXIMAL_KT_UNITS) and attempts < max_attempts:
        attempts += 1
        p = candidate_point(params, rng, effective_radii)
        if p is None:
            continue
        if not centers_nonoverlap(p, all_units(), params.kt_diameter_um):
            continue
        isd = is_distal_point(p)
        idx = len(all_units())
        if isd and len(distal) < DISTAL_KT_UNITS:
            distal.append(KTUnit(index=idx, pos=p, is_distal=True))
        elif (not isd) and len(proximal) < PROXIMAL_KT_UNITS:
            proximal.append(KTUnit(index=idx, pos=p, is_distal=False))

    if len(distal) < DISTAL_KT_UNITS or len(proximal) < PROXIMAL_KT_UNITS:
        raise ValueError(
            f"Could not place {DISTAL_KT_UNITS} distal and {PROXIMAL_KT_UNITS} proximal non-overlapping KT units. "
            f"Try a smaller KT unit diameter."
        )

    units = distal + proximal
    for idx, kt in enumerate(units):
        kt.index = idx
    return units


def clone_kt_units(units: Sequence[KTUnit]) -> List[KTUnit]:
    return [KTUnit(index=u.index, pos=u.pos, is_distal=u.is_distal, alive=True) for u in units]

# -----------------------------------------------------------------------------
# BOSS-style event-driven core
# -----------------------------------------------------------------------------
def tangent_cone_half_angle(cs_pos: Vec3) -> float:
    r_bound = max(NUC_A, NUC_B, NUC_C)
    d = v_len(cs_pos)
    if d <= r_bound:
        return math.radians(89.0)
    return math.asin(min(1.0, r_bound / max(d, 1e-12)))


def max_layout_distance(units: Sequence[KTUnit]) -> float:
    max_dist = 0.0
    for kt in units:
        max_dist = max(max_dist, v_len(v_sub(kt.pos, CS_A)), v_len(v_sub(kt.pos, CS_B)))
    return max_dist


def deflection_to_kt_boundary(cs_pos: Vec3, mt_dir: Vec3, kt_pos: Vec3, kt_radius: float) -> float:
    """BOSS theta: zero inside the KT cone, otherwise angle to its nearest boundary."""
    v = v_sub(kt_pos, cs_pos)
    d = v_len(v)
    if d <= 1e-12:
        return 0.0
    axis = v_mul(v, 1.0 / d)
    theta = math.acos(max(-1.0, min(1.0, v_dot(v_norm(mt_dir), axis))))
    alpha = math.asin(min(1.0, kt_radius / d)) if d > kt_radius else math.pi / 2.0
    return math.degrees(max(0.0, theta - alpha))


def captures_along_path(cs_pos: Vec3, mt_dir: Vec3, units: Sequence[KTUnit], kt_radius: float, max_len: float) -> List[Tuple[int, float]]:
    """Return all alive KT units intersected by the MT path, sorted by contact distance."""
    hits: List[Tuple[int, float]] = []
    d = v_norm(mt_dir)
    for kt in units:
        if not kt.alive:
            continue
        v = v_sub(kt.pos, cs_pos)
        proj = v_dot(v, d)
        if proj < 0:
            continue
        closest = v_add(cs_pos, v_mul(d, proj))
        perp = v_len(v_sub(kt.pos, closest))
        if perp <= kt_radius + 1e-12:
            back = math.sqrt(max(0.0, kt_radius * kt_radius - perp * perp))
            contact = max(0.0, proj - back)
            if contact <= max_len + 1e-12:
                hits.append((kt.index, contact))
    hits.sort(key=lambda item: item[1])
    return hits


def _great_circle_arc_basis(old_dir: Vec3, new_dir: Vec3) -> Tuple[Vec3, Vec3, float]:
    """Return a deterministic basis for the shortest great-circle arc."""
    start = v_norm(old_dir)
    end = v_norm(new_dir)
    dot = max(-1.0, min(1.0, v_dot(start, end)))
    angle = math.acos(dot)
    if angle <= 1e-12:
        return start, orthonormal_basis(start)[0], 0.0
    if math.pi - angle <= 1e-10:
        return start, orthonormal_basis(start)[0], math.pi
    tangent = v_norm(v_sub(end, v_mul(start, dot)))
    return start, tangent, angle


def _direction_on_great_circle_arc(old_dir: Vec3, new_dir: Vec3, fraction: float) -> Vec3:
    start, tangent, angle = _great_circle_arc_basis(old_dir, new_dir)
    if angle <= 1e-12:
        return start
    phi = max(0.0, min(1.0, float(fraction))) * angle
    return v_norm(v_add(v_mul(start, math.cos(phi)), v_mul(tangent, math.sin(phi))))


def _closest_direction_on_great_circle_arc(
    direction: Vec3, old_dir: Vec3, new_dir: Vec3,
) -> Tuple[float, float, Vec3]:
    """Return minimum angle, arc fraction, and orientation closest to direction."""
    q = v_norm(direction)
    start, tangent, angle = _great_circle_arc_basis(old_dir, new_dir)
    start_angle = angular_distance(q, start)
    if angle <= 1e-12:
        return start_angle, 0.0, start
    end = _direction_on_great_circle_arc(old_dir, new_dir, 1.0)
    best_angle, best_fraction, best_direction = start_angle, 0.0, start
    end_angle = angular_distance(q, end)
    if end_angle < best_angle:
        best_angle, best_fraction, best_direction = end_angle, 1.0, end

    plane_normal = v_norm(v_cross(start, tangent))
    projected = v_sub(q, v_mul(plane_normal, v_dot(q, plane_normal)))
    projected_len = v_len(projected)
    if projected_len > 1e-12:
        projected = v_mul(projected, 1.0 / projected_len)
        phi = math.atan2(v_dot(projected, tangent), v_dot(projected, start))
        if phi < 0.0:
            phi += 2.0 * math.pi
        if phi <= angle + 1e-12:
            fraction = max(0.0, min(1.0, phi / angle))
            candidate = _direction_on_great_circle_arc(old_dir, new_dir, fraction)
            candidate_angle = angular_distance(q, candidate)
            if candidate_angle < best_angle:
                best_angle, best_fraction, best_direction = candidate_angle, fraction, candidate
    return best_angle, best_fraction, best_direction


def _minimum_angle_to_great_circle_arc(direction: Vec3, old_dir: Vec3, new_dir: Vec3) -> float:
    return _closest_direction_on_great_circle_arc(direction, old_dir, new_dir)[0]


def _swept_shaft_contact_geometry(
    cs_pos: Vec3,
    old_dir: Vec3,
    new_dir: Vec3,
    kt_pos: Vec3,
    kt_radius: float,
    max_len: float,
) -> Optional[Tuple[float, Vec3, float]]:
    """Return radial contact, first contact orientation, and arc fraction.

    Contact orientation is the first orientation on old->new whose ray touches
    the KT sphere. This is distinct from the post-pivot direction when contact
    occurs partway through a swept sector.
    """
    max_len = max(0.0, float(max_len))
    radius = max(0.0, float(kt_radius))
    offset = v_sub(kt_pos, cs_pos)
    distance = v_len(offset)
    if distance <= radius + 1e-12:
        return (0.0, v_norm(old_dir), 0.0) if max_len >= 0.0 else None
    if distance <= 1e-12 or radius <= 0.0:
        return None

    q_dir = v_mul(offset, 1.0 / distance)
    alpha = math.asin(min(1.0, radius / distance))
    min_angle, min_fraction, _ = _closest_direction_on_great_circle_arc(q_dir, old_dir, new_dir)
    if min_angle > alpha + 1e-12:
        return None

    if angular_distance(q_dir, v_norm(old_dir)) <= alpha + 1e-12:
        fraction = 0.0
        contact_direction = v_norm(old_dir)
    else:
        # Angle to the KT center decreases monotonically from old_dir to the
        # closest point on the arc, so bisection locates first sphere contact.
        lo, hi = 0.0, min_fraction
        for _ in range(64):
            mid = (lo + hi) / 2.0
            candidate = _direction_on_great_circle_arc(old_dir, new_dir, mid)
            if angular_distance(q_dir, candidate) <= alpha:
                hi = mid
            else:
                lo = mid
        fraction = hi
        contact_direction = _direction_on_great_circle_arc(old_dir, new_dir, fraction)

    delta = angular_distance(q_dir, contact_direction)
    cos_delta = math.cos(delta)
    if cos_delta <= 0.0:
        return None
    perpendicular_sq = distance * distance * max(0.0, 1.0 - cos_delta * cos_delta)
    discriminant = radius * radius - perpendicular_sq
    if discriminant < -1e-9:
        return None
    contact = distance * cos_delta - math.sqrt(max(0.0, discriminant))
    if contact < -1e-10 or contact > max_len + 1e-10:
        return None
    return max(0.0, contact), contact_direction, fraction


def _swept_shaft_contact_distance(
    cs_pos: Vec3,
    old_dir: Vec3,
    new_dir: Vec3,
    kt_pos: Vec3,
    kt_radius: float,
    max_len: float,
) -> Optional[float]:
    """Compatibility wrapper returning the swept-sector radial contact."""
    result = _swept_shaft_contact_geometry(cs_pos, old_dir, new_dir, kt_pos, kt_radius, max_len)
    return None if result is None else result[0]


def captures_along_swept_shaft(
    cs_pos: Vec3,
    old_dir: Vec3,
    new_dir: Vec3,
    units: Sequence[KTUnit],
    kt_radius: float,
    max_len: float,
) -> List[Tuple[int, float, Vec3]]:
    """Return alive KTs and true first-contact orientations in a pivot sweep."""
    hits: List[Tuple[int, float, Vec3]] = []
    for kt in units:
        if not kt.alive:
            continue
        result = _swept_shaft_contact_geometry(
            cs_pos, old_dir, new_dir, kt.pos, kt_radius, max_len,
        )
        if result is not None:
            contact, orientation, _fraction = result
            hits.append((kt.index, contact, orientation))
    hits.sort(key=lambda item: (item[1], item[0]))
    return hits


def nearest_alive_deflection(cs_pos: Vec3, mt_dir: Vec3, units: Sequence[KTUnit], kt_radius: float) -> Tuple[Optional[float], Optional[int]]:
    best_angle: Optional[float] = None
    best_idx: Optional[int] = None
    for kt in units:
        if not kt.alive:
            continue
        a = deflection_to_kt_boundary(cs_pos, mt_dir, kt.pos, kt_radius)
        if best_angle is None or a < best_angle:
            best_angle = a
            best_idx = kt.index
    return best_angle, best_idx


def occupied_direction(new_dir: Vec3, dirs_list: Sequence[Vec3], angular_width_rad: float) -> bool:
    return any(angular_distance(new_dir, old) < angular_width_rad for old in dirs_list)


def sample_random_direction_with_occupancy(
    cs_pos: Vec3,
    axis: Vec3,
    half_angle: float,
    dirs_list: List[Vec3],
    angular_width_rad: float,
    rng: random.Random,
) -> Vec3:
    last = sample_cone_direction(axis, half_angle, rng)
    for _ in range(50):
        cand = sample_cone_direction(axis, half_angle, rng)
        last = cand
        if not occupied_direction(cand, dirs_list, angular_width_rad):
            return cand
    return last


def diffuse_pivot_direction(current: Vec3, d_mt_rad2_s: float, dt_s: float, rng: random.Random) -> Vec3:
    """One isotropic angular-diffusion update in the tangent plane.

    Each of the two local tangent components has standard deviation
    sqrt(2 * D_MT * dt), which is the small-angle form of MSAD = 2 D_MT t.
    """
    sigma = math.sqrt(max(0.0, 2.0 * float(d_mt_rad2_s) * max(0.0, float(dt_s))))
    if sigma <= 0.0:
        return v_norm(current)
    u, v, w = orthonormal_basis(current)
    du = rng.gauss(0.0, sigma)
    dv = rng.gauss(0.0, sigma)
    return v_norm(v_add(w, v_add(v_mul(u, du), v_mul(v, dv))))


@dataclass
class PivotingMTState:
    """One MT trajectory in the shared chronological growth engine."""
    record: MTRecord
    current_direction: Vec3
    current_length: float
    last_time: float
    next_pivot_time: float
    max_growth_time: float
    next_theta_sample_time: float
    analysis_samples: List[Tuple[float, Vec3, Tuple[int, ...]]] = field(default_factory=list)
    closest_approach_distance_um: float = float("inf")
    closest_approach_time_s: Optional[float] = None
    closest_approach_kt_index: Optional[int] = None
    closest_approach_alive_indices: Tuple[int, ...] = ()
    closest_approach_direction: Optional[Vec3] = None
    deflection_finalized: bool = False


class EventDrivenSimulator:
    """Shared chronological simulator for Random, Directed, and Pivoting."""

    def __init__(self, params: SimulationParams, directed_prob: float, seed: int,
                 units: Sequence[KTUnit], *, pivot_seed: Optional[int] = None):
        self.params = params
        self.directed_prob = max(0.0, min(1.0, directed_prob))
        # Layout is generated by generate_kt_layout(layout_seed); this stream
        # controls only initial directions/emission target choices.
        self.direction_rng = random.Random(seed)
        self.rng = self.direction_rng  # Backward-compatible internal alias.
        self.pivot_rng = random.Random(seed ^ 0x5EED5EED if pivot_seed is None else pivot_seed)
        self.units = clone_kt_units(units)
        self.max_dist = max_layout_distance(self.units)
        self.max_len = self.max_dist + params.mt_extra_length_um
        self.mt_radius = params.mt_diameter_um / 2.0
        self.angular_width = 2.0 * math.atan(self.mt_radius / max(0.5, self.max_dist * 0.7))
        self.dirs_by_cs: List[List[Vec3]] = [[], []]
        self.axes = [v_norm(v_mul(CS_A, -1.0)), v_norm(v_mul(CS_B, -1.0))]
        self.halves = [tangent_cone_half_angle(CS_A), tangent_cone_half_angle(CS_B)]
        self.mt_records: List[MTRecord] = []
        self.capture_log: List[Dict[str, object]] = []
        self.deflection_rows: List[Dict[str, object]] = []
        self.mt_id_counter = 0
        self.active_mts: List[PivotingMTState] = []
        self._states_by_mt: Dict[int, PivotingMTState] = {}

    def alive_indices(self) -> List[int]:
        return [kt.index for kt in self.units if kt.alive]

    def alive_indices_at(self, sample_time: float) -> Tuple[int, ...]:
        """KTs alive at a time, treating an exact capture-time sample as valid."""
        eps = 1e-10
        return tuple(
            kt.index for kt in self.units
            if kt.alive or (kt.capture_time is not None and kt.capture_time >= sample_time - eps)
        )

    def captured_count(self) -> int:
        return sum(1 for kt in self.units if not kt.alive)

    def completed(self) -> bool:
        return self.captured_count() >= len(self.units)

    def choose_directed_target(self) -> Optional[int]:
        alive = self.alive_indices()
        return self.direction_rng.choice(alive) if alive else None

    def _append_analysis_sample(self, state: PivotingMTState, sample_time: float,
                                direction: Optional[Vec3] = None) -> None:
        if state.deflection_finalized or state.record.hits:
            return
        sample_direction = v_norm(direction if direction is not None else state.current_direction)
        sample = (float(sample_time), sample_direction, self.alive_indices_at(sample_time))
        eps = 1e-10
        if state.analysis_samples and abs(state.analysis_samples[-1][0] - sample_time) <= eps:
            state.analysis_samples[-1] = sample
        else:
            state.analysis_samples.append(sample)

    def _track_uncaptured_closest_segment(
        self, state: PivotingMTState, start_time: float, start_length: float,
        end_time: float, end_length: float, direction: Vec3,
    ) -> None:
        """Track true plus-end/KT Euclidean minimum on a growing straight segment."""
        if state.deflection_finalized or state.record.hits:
            return
        cs_pos = CENTROSOMES[state.record.cs_index]
        direction = v_norm(direction)
        velocity = max(1e-12, self.params.mt_velocity_um_s)
        lo, hi = min(start_length, end_length), max(start_length, end_length)
        for kt in self.units:
            kt_hi = hi
            if kt.capture_time is not None:
                if kt.capture_time < start_time - 1e-10:
                    continue
                kt_hi = min(kt_hi, max(lo, (kt.capture_time - state.record.emit_time) * velocity))
            if kt_hi < lo - 1e-10:
                continue
            offset = v_sub(kt.pos, cs_pos)
            closest_length = max(lo, min(kt_hi, v_dot(offset, direction)))
            plus_end = v_add(cs_pos, v_mul(direction, closest_length))
            distance = v_len(v_sub(kt.pos, plus_end))
            if distance < state.closest_approach_distance_um - 1e-12:
                closest_time = state.record.emit_time + closest_length / velocity
                alive = set(self.alive_indices_at(closest_time))
                alive.add(kt.index)  # Preserve a terminal KT captured at this exact time.
                state.closest_approach_distance_um = distance
                state.closest_approach_time_s = closest_time
                state.closest_approach_kt_index = kt.index
                state.closest_approach_alive_indices = tuple(sorted(alive))
                state.closest_approach_direction = direction

    def _track_uncaptured_closest_arc(
        self, state: PivotingMTState, old_direction: Vec3,
        new_direction: Vec3, pivot_time: float,
    ) -> None:
        """Include the continuous plus-end arc made by a direction update."""
        if state.deflection_finalized or state.record.hits:
            return
        cs_pos = CENTROSOMES[state.record.cs_index]
        length = state.current_length
        for kt in self.units:
            if kt.capture_time is not None and kt.capture_time < pivot_time - 1e-10:
                continue
            offset = v_sub(kt.pos, cs_pos)
            center_distance = v_len(offset)
            if center_distance <= 1e-12:
                direction, distance = v_norm(old_direction), length
            else:
                _delta, _fraction, direction = _closest_direction_on_great_circle_arc(
                    v_mul(offset, 1.0 / center_distance), old_direction, new_direction,
                )
                distance_sq = (
                    center_distance * center_distance + length * length
                    - 2.0 * center_distance * length
                    * max(-1.0, min(1.0, v_dot(v_mul(offset, 1.0 / center_distance), direction)))
                )
                distance = math.sqrt(max(0.0, distance_sq))
            if distance < state.closest_approach_distance_um - 1e-12:
                alive = set(self.alive_indices_at(pivot_time))
                alive.add(kt.index)
                state.closest_approach_distance_um = distance
                state.closest_approach_time_s = pivot_time
                state.closest_approach_kt_index = kt.index
                state.closest_approach_alive_indices = tuple(sorted(alive))
                state.closest_approach_direction = direction

    def _advance_state(self, state: PivotingMTState, target_time: float) -> None:
        """Advance growth and fixed 0.25 s analysis samples to a global event."""
        eps = 1e-10
        target_time = max(state.last_time, float(target_time))
        velocity = max(1e-12, self.params.mt_velocity_um_s)
        new_length = min(state.record.max_len, max(0.0, (target_time - state.record.emit_time) * velocity))
        self._track_uncaptured_closest_segment(
            state, state.last_time, state.current_length, target_time,
            new_length, state.current_direction,
        )
        if not state.deflection_finalized and not state.record.hits:
            while state.next_theta_sample_time <= target_time + eps:
                self._append_analysis_sample(state, state.next_theta_sample_time)
                state.next_theta_sample_time += THETA_SAMPLE_DT_S
        state.current_length = new_length
        state.last_time = target_time

    def _next_growth_contacts(
        self, state: PivotingMTState, limit_time: float,
    ) -> List[Tuple[float, float, int]]:
        """Find this MT's earliest alive-KT plus-end contacts before limit_time."""
        velocity = max(1e-12, self.params.mt_velocity_um_s)
        max_length = min(state.record.max_len, max(0.0, (limit_time - state.record.emit_time) * velocity))
        hits = captures_along_path(
            CENTROSOMES[state.record.cs_index], state.current_direction,
            self.units, self.params.kt_radius_um, max_length,
        )
        eps = 1e-10
        candidates = [
            (max(state.last_time, state.record.emit_time + contact / velocity), contact, kt_idx)
            for kt_idx, contact in hits
            if contact + eps >= state.current_length
        ]
        if not candidates:
            return []
        candidates.sort(key=lambda item: (item[0], item[1], item[2]))
        first_time = candidates[0][0]
        return [c for c in candidates if abs(c[0] - first_time) <= eps]

    def _next_global_growth_contacts(
        self, current_time: float, physical_limit: float,
    ) -> Tuple[float, List[Tuple[PivotingMTState, float, float, int]]]:
        best_time = float("inf")
        candidates: List[Tuple[PivotingMTState, float, float, int]] = []
        eps = 1e-10
        for state in self.active_mts:
            for capture_time, contact, kt_idx in self._next_growth_contacts(state, physical_limit):
                if capture_time < current_time - eps:
                    continue
                if capture_time < best_time - eps:
                    best_time = capture_time
                    candidates = [(state, capture_time, contact, kt_idx)]
                elif abs(capture_time - best_time) <= eps:
                    candidates.append((state, capture_time, contact, kt_idx))
        return best_time, candidates

    def _nearest_deflection_from_indices(
        self, cs_pos: Vec3, direction: Vec3, kt_indices: Sequence[int],
    ) -> Tuple[Optional[float], Optional[int]]:
        """Find the angularly nearest KT in a sample's alive set."""
        best_angle: Optional[float] = None
        best_idx: Optional[int] = None
        for kt_idx in kt_indices:
            kt = self.units[kt_idx]
            angle = deflection_to_kt_boundary(cs_pos, direction, kt.pos, self.params.kt_radius_um)
            if best_angle is None or angle < best_angle:
                best_angle, best_idx = angle, kt_idx
        return best_angle, best_idx

    def _finalize_captured_mt(
        self, state: PivotingMTState, kt_idx: int, contact: float,
        capture_time: float, capture_mode: str, terminal_direction: Vec3,
    ) -> None:
        if state.deflection_finalized:
            return
        eps = 1e-10
        samples = [s for s in state.analysis_samples if s[0] < capture_time - eps]
        exact = next((s for s in reversed(state.analysis_samples) if abs(s[0] - capture_time) <= eps), None)
        terminal = (capture_time, v_norm(terminal_direction), self.alive_indices_at(capture_time))
        if exact is not None:
            samples.append(terminal)
        else:
            samples.append(terminal)
        if not samples:
            samples = [(state.record.emit_time, state.record.direction, self.alive_indices_at(state.record.emit_time)), terminal]
        cs_pos = CENTROSOMES[state.record.cs_index]
        kt = self.units[kt_idx]
        angles = [
            deflection_to_kt_boundary(cs_pos, direction, kt.pos, self.params.kt_radius_um)
            for _time, direction, _alive in samples
        ]
        mean_theta = statistics.fmean(angles)
        state.deflection_finalized = True
        state.record.mean_deflection_angle = mean_theta
        state.record.nearest_kt_index = kt_idx
        state.record.final_direction = state.current_direction
        self.deflection_rows.append({
            "mt_id": state.record.mt_id,
            "mt_status": "captured",
            "mean_deflection_angle_deg": mean_theta,
            "deflection_to_boundary_deg": mean_theta,
            "nearest_kt_index": kt_idx,
            "nearest_is_distal": kt.is_distal,
            "nearest_kt_group": kt.group,
            "cs_index": state.record.cs_index,
            "emit_time_s": state.record.emit_time,
            "capture_time_s": capture_time,
            "analysis_start_time_s": state.record.emit_time,
            "analysis_end_time_s": capture_time,
            "analysis_sample_count": len(samples),
            "analysis_sampling_interval_s": THETA_SAMPLE_DT_S,
            "analysis_terminal_rule": "first KT capture/engagement",
            "analysis_nearest_rule": "the first KT captured by this MT, fixed across its trajectory",
            "contact_distance_um": contact,
            "capture_mode": capture_mode,
            "first_capture_mode": capture_mode,
            "closest_approach_time_s": "",
            "closest_approach_distance_um": "",
            "closest_approach_kt_index": "",
            "is_directed": state.record.is_directed,
            "target_kt_index": "" if state.record.target_index is None else state.record.target_index,
            "target_is_distal": "" if state.record.target_index is None else self.units[state.record.target_index].is_distal,
            "target_kt_group": "" if state.record.target_index is None else self.units[state.record.target_index].group,
            "dir_x": terminal_direction[0], "dir_y": terminal_direction[1], "dir_z": terminal_direction[2],
            "pivot_step_count": state.record.pivot_step_count,
        })
        state.analysis_samples.clear()

    def _append_capture(
        self, state: PivotingMTState, kt_idx: int, contact: float,
        capture_time: float, capture_mode: str, contact_direction: Vec3,
    ) -> None:
        kt = self.units[kt_idx]
        if not kt.alive:
            return
        kt.alive = False
        kt.capture_time = capture_time
        capture_index = self.captured_count()
        state.record.hits.append((kt_idx, contact))
        state.record.hit_capture_modes.append(capture_mode)
        self.capture_log.append({
            "capture_index": capture_index,
            "captured_fraction": capture_index / float(len(self.units)),
            "capture_time_s": capture_time,
            "mt_total_emitted_at_capture": self.mt_id_counter,
            "kt_index": kt_idx,
            "is_distal": kt.is_distal,
            "kt_group": kt.group,
            "mt_id": state.record.mt_id,
            "cs_index": state.record.cs_index,
            "emit_time_s": state.record.emit_time,
            "contact_distance_um": contact,
            "capture_mode": capture_mode,
            "pivoting_mode": self.params.pivoting_mode,
            "d_mt_rad2_s": self.params.d_mt_rad2_s,
        })
        self._finalize_captured_mt(state, kt_idx, contact, capture_time, capture_mode, contact_direction)

    def _finalize_uncaptured_mt(self, state: PivotingMTState) -> None:
        if state.deflection_finalized or state.record.hits:
            return
        terminal_time = state.closest_approach_time_s
        closest_idx = state.closest_approach_kt_index
        if terminal_time is None or closest_idx is None:
            # This is a defensive fallback for degenerate layouts with no alive
            # KT at emission; every analyzed MT still receives one result row.
            terminal_time = state.record.emit_time
            closest_idx = self.units[0].index if self.units else None
            state.closest_approach_alive_indices = tuple(self.alive_indices_at(terminal_time))
            state.closest_approach_direction = state.record.direction
            state.closest_approach_distance_um = float("nan")
        if closest_idx is None:
            raise RuntimeError(f"Cannot finalize MT {state.record.mt_id}: no KT units exist.")
        terminal_direction = v_norm(state.closest_approach_direction or state.record.direction)
        eps = 1e-10
        samples = [sample for sample in state.analysis_samples if sample[0] < terminal_time - eps]
        terminal_alive = set(state.closest_approach_alive_indices)
        terminal_alive.add(closest_idx)
        terminal_sample = (terminal_time, terminal_direction, tuple(sorted(terminal_alive)))
        samples.append(terminal_sample)
        if not samples:
            samples = [(terminal_time, terminal_direction, tuple(sorted(terminal_alive)))]

        cs_pos = CENTROSOMES[state.record.cs_index]
        theta_values: List[float] = []
        terminal_nearest_idx: Optional[int] = None
        for sample_time, direction, alive_indices in samples:
            angle, nearest_idx = self._nearest_deflection_from_indices(cs_pos, direction, alive_indices)
            if angle is None or nearest_idx is None:
                angle = deflection_to_kt_boundary(
                    cs_pos, direction, self.units[closest_idx].pos, self.params.kt_radius_um,
                )
                nearest_idx = closest_idx
            theta_values.append(angle)
            if abs(sample_time - terminal_time) <= eps:
                terminal_nearest_idx = nearest_idx
        if not theta_values:
            raise RuntimeError(f"Cannot finalize MT {state.record.mt_id}: no theta samples were available.")
        mean_theta = statistics.fmean(theta_values)
        terminal_nearest_idx = terminal_nearest_idx if terminal_nearest_idx is not None else closest_idx
        nearest = self.units[terminal_nearest_idx]
        state.deflection_finalized = True
        state.record.no_hit_deflection = mean_theta
        state.record.mean_deflection_angle = mean_theta
        state.record.nearest_kt_index = terminal_nearest_idx
        state.record.final_direction = state.current_direction
        self.deflection_rows.append({
            "mt_id": state.record.mt_id,
            "mt_status": "uncaptured",
            "mean_deflection_angle_deg": mean_theta,
            "deflection_to_boundary_deg": mean_theta,
            "nearest_kt_index": terminal_nearest_idx,
            "nearest_is_distal": nearest.is_distal,
            "nearest_kt_group": nearest.group,
            "cs_index": state.record.cs_index,
            "emit_time_s": state.record.emit_time,
            "capture_time_s": "",
            "analysis_start_time_s": state.record.emit_time,
            "analysis_end_time_s": terminal_time,
            "analysis_sample_count": len(theta_values),
            "analysis_sampling_interval_s": THETA_SAMPLE_DT_S,
            "analysis_terminal_rule": "minimum Euclidean plus-end-to-KT distance",
            "analysis_nearest_rule": "angularly nearest KT alive at each sample; terminal classification uses the terminal sample",
            "contact_distance_um": "",
            "capture_mode": "",
            "first_capture_mode": "",
            "closest_approach_time_s": terminal_time,
            "closest_approach_distance_um": state.closest_approach_distance_um,
            "closest_approach_kt_index": closest_idx,
            "is_directed": state.record.is_directed,
            "target_kt_index": "" if state.record.target_index is None else state.record.target_index,
            "target_is_distal": "" if state.record.target_index is None else self.units[state.record.target_index].is_distal,
            "target_kt_group": "" if state.record.target_index is None else self.units[state.record.target_index].group,
            "dir_x": terminal_direction[0], "dir_y": terminal_direction[1], "dir_z": terminal_direction[2],
            "pivot_step_count": state.record.pivot_step_count,
        })
        state.analysis_samples.clear()

    def _emit_one(self, cs_index: int, emit_time: float) -> MTRecord:
        cs_pos = CENTROSOMES[cs_index]
        is_directed = (self.directed_prob > 0.0) and (self.direction_rng.random() < self.directed_prob)
        target_idx: Optional[int] = None
        if is_directed:
            target_idx = self.choose_directed_target()
            if target_idx is not None:
                direction = v_norm(v_sub(self.units[target_idx].pos, cs_pos))
            else:
                direction = sample_random_direction_with_occupancy(
                    cs_pos, self.axes[cs_index], self.halves[cs_index],
                    self.dirs_by_cs[cs_index], self.angular_width, self.direction_rng,
                )
        else:
            direction = sample_random_direction_with_occupancy(
                cs_pos, self.axes[cs_index], self.halves[cs_index],
                self.dirs_by_cs[cs_index], self.angular_width, self.direction_rng,
            )
        self.dirs_by_cs[cs_index].append(direction)
        self.mt_id_counter += 1
        mt = MTRecord(
            mt_id=self.mt_id_counter, cs_index=cs_index, emit_time=emit_time,
            direction=direction, is_directed=is_directed, target_index=target_idx,
            max_len=self.max_len,
        )
        velocity = max(1e-12, self.params.mt_velocity_um_s)
        pivot_enabled = self.params.pivoting_enabled and self.params.d_mt_rad2_s > 0.0
        state = PivotingMTState(
            record=mt, current_direction=direction, current_length=0.0,
            last_time=emit_time,
            next_pivot_time=(emit_time + max(1e-12, self.params.pivot_dt_s)) if pivot_enabled else float("inf"),
            max_growth_time=emit_time + mt.max_len / velocity,
            next_theta_sample_time=emit_time + THETA_SAMPLE_DT_S,
            analysis_samples=[(emit_time, direction, self.alive_indices_at(emit_time))],
        )
        self._track_uncaptured_closest_segment(state, emit_time, 0.0, emit_time, 0.0, direction)
        self.mt_records.append(mt)
        self.active_mts.append(state)
        self._states_by_mt[mt.mt_id] = state
        return mt

    def emit_one(self, cs_index: int, emit_time: float) -> MTRecord:
        """Backward-compatible public name for emitting a chronological MT."""
        return self._emit_one(cs_index, emit_time)

    def _process_growth_candidates(
        self, candidates: Sequence[Tuple[PivotingMTState, float, float, int]],
    ) -> None:
        for state, capture_time, contact, kt_idx in sorted(
            candidates, key=lambda c: (c[0].record.mt_id, c[1], c[2], c[3]),
        ):
            if self.units[kt_idx].alive:
                self._append_capture(
                    state, kt_idx, contact, capture_time, "plus_end", state.current_direction,
                )

    def _pivot_states_at(self, event_time: float) -> List[PivotingMTState]:
        eps = 1e-10
        return [
            state for state in self.active_mts
            if state.next_pivot_time <= event_time + eps
            and state.next_pivot_time <= state.max_growth_time + eps
        ]

    def _process_pivots(self, states: Sequence[PivotingMTState], event_time: float) -> None:
        if not states:
            return
        changed: List[Tuple[PivotingMTState, Vec3, Vec3]] = []
        for state in sorted(states, key=lambda s: s.record.mt_id):
            old_direction = state.current_direction
            new_direction = diffuse_pivot_direction(
                old_direction, self.params.d_mt_rad2_s, self.params.pivot_dt_s,
                self.pivot_rng,
            )
            state.current_direction = new_direction
            state.record.pivot_step_count += 1
            if not state.record.hits and not state.deflection_finalized:
                self._track_uncaptured_closest_arc(state, old_direction, new_direction, event_time)
                self._append_analysis_sample(state, event_time, new_direction)
            changed.append((state, old_direction, new_direction))
            state.next_pivot_time += max(1e-12, self.params.pivot_dt_s)

        swept_candidates: List[Tuple[PivotingMTState, float, float, int, Vec3]] = []
        for state, old_direction, new_direction in changed:
            # D_MT=0 has exactly the fixed-ray collision set; it receives no
            # extra shaft sweep when old and new orientations are identical.
            if angular_distance(old_direction, new_direction) <= 1e-12:
                continue
            hits = captures_along_swept_shaft(
                CENTROSOMES[state.record.cs_index], old_direction, new_direction,
                self.units, self.params.kt_radius_um, state.current_length,
            )
            swept_candidates.extend(
                (state, event_time, contact, kt_idx, contact_direction)
                for kt_idx, contact, contact_direction in hits
            )
        for state, capture_time, contact, kt_idx, contact_direction in sorted(
            swept_candidates, key=lambda c: (c[0].record.mt_id, c[2], c[3]),
        ):
            if self.units[kt_idx].alive:
                self._append_capture(
                    state, kt_idx, contact, capture_time, "lateral_shaft", contact_direction,
                )

    def _run_chronological(self, max_mts: int, max_time_s: float) -> Dict[str, object]:
        """Run every condition through the same chronological event loop."""
        eps = 1e-10
        for _ in range(FIRST_BURST_PER_CENTROSOME):
            if not self.completed():
                self._emit_one(0, 0.0)
            if not self.completed():
                self._emit_one(1, 0.0)

        emission_dt = 1.0 / max(1e-12, self.params.emission_rate_per_centrosome)
        next_emit_time = emission_dt
        current_time = 0.0
        max_time_s = max(0.0, float(max_time_s))

        while not self.completed() and current_time < max_time_s - eps:
            can_emit = self.mt_id_counter < max_mts
            physical_events = [max_time_s]
            if can_emit:
                physical_events.append(next_emit_time)
            for state in self.active_mts:
                physical_events.append(state.max_growth_time)
                if math.isfinite(state.next_pivot_time):
                    physical_events.append(state.next_pivot_time)
            physical_limit = min(physical_events)
            collision_time, collision_candidates = self._next_global_growth_contacts(
                current_time, physical_limit,
            )
            event_time = min(physical_limit, collision_time)
            if not math.isfinite(event_time) or event_time > max_time_s + eps:
                break

            if event_time > current_time + eps:
                for state in list(self.active_mts):
                    self._advance_state(state, event_time)
            current_time = max(current_time, event_time)

            if collision_candidates and abs(collision_time - event_time) <= eps:
                self._process_growth_candidates(collision_candidates)

            self._process_pivots(self._pivot_states_at(event_time), event_time)

            for state in list(self.active_mts):
                if state.current_length >= state.record.max_len - eps:
                    if not state.record.hits:
                        self._finalize_uncaptured_mt(state)
                    else:
                        state.record.final_direction = state.current_direction
            self.active_mts = [
                state for state in self.active_mts
                if state.current_length < state.record.max_len - eps
            ]

            # Existing MT events win ties with new emissions; this preserves
            # chronological capture ordering at a shared timestamp.
            if self.completed():
                break
            if self.mt_id_counter < max_mts and next_emit_time <= current_time + eps:
                self._emit_one(0, next_emit_time)
                if not self.completed() and self.mt_id_counter < max_mts:
                    self._emit_one(1, next_emit_time)
                next_emit_time += emission_dt
            if current_time >= max_time_s - eps:
                break

        # All emitted MTs are in the deflection analysis. At simulation stop,
        # close any still-growing uncaptured trajectory at its observed closest
        # approach; captured rows were finalized at first engagement.
        for state in self.active_mts:
            if not state.record.hits:
                self._finalize_uncaptured_mt(state)
            state.record.final_direction = state.current_direction
        rows_by_mt: Dict[int, int] = {}
        for row in self.deflection_rows:
            rows_by_mt[int(row["mt_id"])] = rows_by_mt.get(int(row["mt_id"]), 0) + 1
        expected_ids = {mt.mt_id for mt in self.mt_records}
        missing = sorted(expected_ids - set(rows_by_mt))
        duplicates = sorted(mt_id for mt_id, count in rows_by_mt.items() if count != 1)
        if missing or duplicates or len(self.deflection_rows) != len(self.mt_records):
            raise RuntimeError(
                "Deflection QC failed: expected exactly one theta_mean row per emitted/analyzed MT; "
                f"MTs={len(self.mt_records)}, rows={len(self.deflection_rows)}, "
                f"missing={missing[:10]}, duplicate_ids={duplicates[:10]}."
            )

        completed = self.completed()
        final_time = max((float(row["capture_time_s"]) for row in self.capture_log), default=0.0)
        return {
            "completed": completed,
            "n_mts": self.mt_id_counter,
            "n_captured": self.captured_count(),
            "n_kt": len(self.units),
            "final_time_s": final_time,
            "capture_log": self.capture_log,
            "deflection_rows": self.deflection_rows,
            "mt_records": self.mt_records,
            "analysis_mt_count": len(self.mt_records),
            "deflection_rows_count": len(self.deflection_rows),
            "deflection_qc_passed": True,
        }

    def run(self, max_mts: int = 500000, max_time_s: float = 200000.0) -> Dict[str, object]:
        return self._run_chronological(max_mts, max_time_s)


def run_random_vs_zero_dmt_qc(
    params: SimulationParams, seed: int, units: Sequence[KTUnit],
    directed_prob: float = 0.0, max_mts: int = 500000,
    max_time_s: float = 200000.0,
) -> Dict[str, object]:
    """Run paired Random and D_MT=0 simulations and assert semantic identity."""
    random_params = replace(params, pivoting_mode="none", d_mt_rad2_s=0.0)
    zero_dmt_params = replace(params, pivoting_mode="dmt", d_mt_rad2_s=0.0)
    random_result = EventDrivenSimulator(random_params, directed_prob, seed, units).run(max_mts, max_time_s)
    zero_result = EventDrivenSimulator(zero_dmt_params, directed_prob, seed, units).run(max_mts, max_time_s)
    capture_fields = (
        "capture_index", "capture_time_s", "mt_total_emitted_at_capture", "kt_index",
        "mt_id", "cs_index", "emit_time_s", "contact_distance_um", "capture_mode",
    )
    deflection_fields = (
        "mt_id", "mt_status", "mean_deflection_angle_deg", "analysis_start_time_s",
        "analysis_end_time_s", "analysis_sample_count", "nearest_kt_index",
    )
    def normalize(rows, fields):
        return [tuple(row.get(field) for field in fields) for row in rows]
    summaries_equal = all(
        random_result[key] == zero_result[key]
        for key in ("completed", "n_mts", "n_captured", "final_time_s")
    )
    captures_equal = normalize(random_result["capture_log"], capture_fields) == normalize(zero_result["capture_log"], capture_fields)
    deflections_equal = normalize(random_result["deflection_rows"], deflection_fields) == normalize(zero_result["deflection_rows"], deflection_fields)
    if not (summaries_equal and captures_equal and deflections_equal):
        raise AssertionError(
            "Random vs D_MT=0 paired QC failed: shared layout/initial RNG did not produce identical outputs."
        )
    theta_values = [float(row["mean_deflection_angle_deg"]) for row in random_result["deflection_rows"]]
    return {
        "passed": True,
        "summary_equal": summaries_equal,
        "capture_progression_equal": captures_equal,
        "per_mt_theta_mean_equal": deflections_equal,
        "completed": bool(random_result["completed"]),
        "completion_time_s": float(random_result["final_time_s"]),
        "n_mts": int(random_result["n_mts"]),
        "capture_event_count": len(random_result["capture_log"]),
        "deflection_row_count": len(random_result["deflection_rows"]),
        "mean_theta_mean_deg": statistics.fmean(theta_values) if theta_values else float("nan"),
    }


def run_pivot_timestep_convergence(
    params: SimulationParams, directed_prob: float, seed: int,
    units: Sequence[KTUnit], d_mt_rad2_s: float,
    pivot_dts: Sequence[float] = (1.0, 0.5, 0.25),
    max_mts: int = 500000, max_time_s: float = 200000.0,
) -> List[Dict[str, object]]:
    """Diagnostic D_MT run across pivot timesteps; theta sampling stays 0.25 s."""
    reports: List[Dict[str, object]] = []
    for pivot_dt in pivot_dts:
        sim_params = replace(
            params, pivoting_mode="dmt", d_mt_rad2_s=max(0.0, float(d_mt_rad2_s)),
            pivot_dt_s=max(1e-12, float(pivot_dt)),
        )
        result = EventDrivenSimulator(sim_params, directed_prob, seed, units).run(max_mts, max_time_s)
        angles = [float(row["mean_deflection_angle_deg"]) for row in result["deflection_rows"]]
        reports.append({
            "pivot_dt_s": float(pivot_dt),
            "theta_sample_dt_s": THETA_SAMPLE_DT_S,
            "completed": bool(result["completed"]),
            "completion_time_s": float(result["final_time_s"]) if result["completed"] else "",
            "mt_total_emitted": int(result["n_mts"]),
            "mean_theta_mean_deg": statistics.fmean(angles) if angles else float("nan"),
            "deflection_rows_count": int(result["deflection_rows_count"]),
            "analysis_mt_count": int(result["analysis_mt_count"]),
            "qc_passed": bool(result["deflection_qc_passed"]),
        })
    return reports


# -----------------------------------------------------------------------------
# Animation simulation: visualizes growth with the same plus-end and swept-
# shaft capture rules as the event-driven core.
# -----------------------------------------------------------------------------
@dataclass
class VisualKT:
    index: int
    pos: Vec3
    is_distal: bool
    reserved: bool = False
    displayed_captured: bool = False
    capture_time: Optional[float] = None
    display_pos: Optional[Vec3] = None
    retracting: bool = False
    retract_start_time: Optional[float] = None
    retract_end_time: Optional[float] = None
    retract_from_pos: Optional[Vec3] = None
    retract_to_pos: Optional[Vec3] = None
    captured_by_cs_index: Optional[int] = None

    def __post_init__(self) -> None:
        if self.display_pos is None:
            self.display_pos = self.pos

    @property
    def group(self) -> str:
        return "distal" if self.is_distal else "proximal"


@dataclass
class VisualMT:
    mt_id: int
    cs_index: int
    emit_time: float
    direction: Vec3
    max_len: float
    capture_contacts: List[Tuple[int, float, float]]  # kt_index, contact_distance, capture_time
    is_directed: bool
    target_index: Optional[int]
    archived: bool = False
    current_length: float = 0.0
    last_update_time: Optional[float] = None
    next_pivot_time: Optional[float] = None
    max_growth_time: Optional[float] = None
    pivot_step_count: int = 0


class AnimationSimulation:
    def __init__(self, params: SimulationParams, directed_prob: float, units: Sequence[KTUnit], seed: int):
        self.params = params
        self.directed_prob = directed_prob
        self.direction_rng = random.Random(seed)
        self.pivot_rng = random.Random(seed ^ 0x5EED5EED)
        self.rng = self.direction_rng
        self.kt_units: List[VisualKT] = [VisualKT(u.index, u.pos, u.is_distal) for u in units]
        self.max_dist = max_layout_distance([KTUnit(k.index, k.pos, k.is_distal) for k in self.kt_units])
        self.max_len = self.max_dist + params.mt_extra_length_um
        self.mt_radius = params.mt_diameter_um / 2.0
        self.angular_width = 2.0 * math.atan(self.mt_radius / max(0.5, self.max_dist * 0.7))
        self.dirs_by_cs: List[List[Vec3]] = [[], []]
        self.axes = [v_norm(v_mul(CS_A, -1.0)), v_norm(v_mul(CS_B, -1.0))]
        self.halves = [tangent_cone_half_angle(CS_A), tangent_cone_half_angle(CS_B)]
        self.time = 0.0
        self.next_emit_time = 1.0 / max(1e-12, params.emission_rate_per_centrosome)
        self.mt_id_counter = 0
        self.mts: List[VisualMT] = []
        self.archived_mts: List[VisualMT] = []
        self.first_burst_done = False
        self.retract_duration_s = 1.0

    def retraction_target_pos(self, kt_index: int, cs_index: int) -> Vec3:
        """Deterministic low-amplitude offset near the capturing centrosome for display only."""
        rng = random.Random(20260522 + kt_index * 37 + cs_index * 1009)
        d = random_unit_vector(rng)
        # Keep captured KT units close to the centrosome, but not exactly on top of it.
        offset = rng.uniform(0.38, 0.82)
        return v_add(CENTROSOMES[cs_index], v_mul(d, offset))

    def reserved_alive_indices(self) -> List[int]:
        return [k.index for k in self.kt_units if not k.reserved]

    def completed(self) -> bool:
        return all(k.reserved for k in self.kt_units)

    def make_core_units_for_collision(self) -> List[KTUnit]:
        # Reuse event collision code, mapping reserved=False to alive=True.
        return [KTUnit(k.index, k.pos, k.is_distal, alive=(not k.reserved)) for k in self.kt_units]

    def _append_visual_pivot_capture(
        self,
        mt: VisualMT,
        kt_idx: int,
        contact: float,
        capture_time: float,
    ) -> None:
        kt = self.kt_units[kt_idx]
        if kt.reserved:
            return
        kt.reserved = True
        kt.capture_time = capture_time
        mt.capture_contacts.append((kt_idx, contact, capture_time))

    def _visual_pivot_direction(self, mt: VisualMT) -> Vec3:
        return diffuse_pivot_direction(
            mt.direction,
            self.params.d_mt_rad2_s,
            self.params.pivot_dt_s,
            self.pivot_rng,
        )

    def _visual_growth_contacts(self, mt: VisualMT, limit_time: float) -> List[Tuple[float, float, int]]:
        """Find the MT's first plus-end contacts before the next global event."""
        velocity = max(1e-12, self.params.mt_velocity_um_s)
        max_length = min(mt.max_len, max(0.0, (limit_time - mt.emit_time) * velocity))
        hits = captures_along_path(
            CENTROSOMES[mt.cs_index], mt.direction, self.make_core_units_for_collision(),
            self.params.kt_radius_um, max_length,
        )
        eps = 1e-10
        candidates = [
            (max(mt.last_update_time or mt.emit_time, mt.emit_time + contact / velocity), contact, kt_idx)
            for kt_idx, contact in hits if contact + eps >= mt.current_length
        ]
        if not candidates:
            return []
        candidates.sort(key=lambda item: (item[0], item[1], item[2]))
        first_time = candidates[0][0]
        return [candidate for candidate in candidates if abs(candidate[0] - first_time) <= eps]

    def _advance_visual_physics_to(self, target_time: float) -> None:
        """Advance all visible MTs with one global chronological event order."""
        eps = 1e-10
        velocity = max(1e-12, self.params.mt_velocity_um_s)
        current_time = self.time
        while current_time < target_time - eps:
            events = [target_time]
            can_emit = not self.completed()
            if can_emit:
                events.append(self.next_emit_time)
            for mt in self.mts:
                if mt.current_length < mt.max_len - eps and mt.max_growth_time is not None:
                    events.append(mt.max_growth_time)
                if mt.next_pivot_time is not None and math.isfinite(mt.next_pivot_time):
                    events.append(mt.next_pivot_time)
            physical_limit = min(events)

            collision_time = float("inf")
            collision_candidates: List[Tuple[VisualMT, float, float, int]] = []
            for mt in self.mts:
                if mt.current_length >= mt.max_len - eps:
                    continue
                for cap_time, contact, kt_idx in self._visual_growth_contacts(mt, physical_limit):
                    if cap_time < current_time - eps:
                        continue
                    if cap_time < collision_time - eps:
                        collision_time = cap_time
                        collision_candidates = [(mt, cap_time, contact, kt_idx)]
                    elif abs(cap_time - collision_time) <= eps:
                        collision_candidates.append((mt, cap_time, contact, kt_idx))

            event_time = min(physical_limit, collision_time)
            for mt in self.mts:
                mt.current_length = min(mt.max_len, max(0.0, (event_time - mt.emit_time) * velocity))
                mt.last_update_time = event_time
            current_time = event_time

            if collision_candidates and abs(collision_time - event_time) <= eps:
                for mt, cap_time, contact, kt_idx in sorted(
                    collision_candidates, key=lambda c: (c[0].mt_id, c[1], c[2], c[3]),
                ):
                    if not self.kt_units[kt_idx].reserved:
                        self._append_visual_pivot_capture(mt, kt_idx, contact, cap_time)

            pivoting_states = [
                mt for mt in self.mts
                if mt.next_pivot_time is not None
                and mt.next_pivot_time <= event_time + eps
                and mt.next_pivot_time <= (mt.max_growth_time or float("inf")) + eps
            ]
            swept_candidates: List[Tuple[VisualMT, float, int]] = []
            for mt in sorted(pivoting_states, key=lambda item: item.mt_id):
                old_direction = mt.direction
                new_direction = self._visual_pivot_direction(mt)
                mt.direction = new_direction
                mt.pivot_step_count += 1
                mt.next_pivot_time += max(1e-12, self.params.pivot_dt_s)
                if angular_distance(old_direction, new_direction) <= 1e-12:
                    continue
                hits = captures_along_swept_shaft(
                    CENTROSOMES[mt.cs_index], old_direction, new_direction,
                    self.make_core_units_for_collision(), self.params.kt_radius_um,
                    mt.current_length,
                )
                swept_candidates.extend(
                    (mt, contact, kt_idx)
                    for kt_idx, contact, _contact_direction in hits
                )
            for mt, contact, kt_idx in sorted(
                swept_candidates, key=lambda c: (c[0].mt_id, c[1], c[2]),
            ):
                if not self.kt_units[kt_idx].reserved:
                    self._append_visual_pivot_capture(mt, kt_idx, contact, event_time)

            if self.completed() is False and self.next_emit_time <= current_time + eps:
                self.emit_visual(0, self.next_emit_time)
                self.emit_visual(1, self.next_emit_time)
                self.next_emit_time += 1.0 / max(1e-12, self.params.emission_rate_per_centrosome)
            if event_time >= target_time - eps:
                break
            if event_time <= current_time - eps:
                break

    def emit_visual(self, cs_index: int, emit_time: float) -> Optional[VisualMT]:
        if self.completed():
            return None
        cs_pos = CENTROSOMES[cs_index]
        is_directed = (self.directed_prob > 0.0) and (self.direction_rng.random() < self.directed_prob)
        target_idx: Optional[int] = None
        if is_directed and self.reserved_alive_indices():
            target_idx = self.direction_rng.choice(self.reserved_alive_indices())
            direction = v_norm(v_sub(self.kt_units[target_idx].pos, cs_pos))
        else:
            direction = sample_random_direction_with_occupancy(
                cs_pos, self.axes[cs_index], self.halves[cs_index], self.dirs_by_cs[cs_index], self.angular_width, self.direction_rng
            )
            is_directed = False if target_idx is None else is_directed
        self.dirs_by_cs[cs_index].append(direction)
        self.mt_id_counter += 1
        pivot_enabled = self.params.pivoting_enabled and self.params.d_mt_rad2_s > 0.0
        mt = VisualMT(
            self.mt_id_counter, cs_index, emit_time, direction, self.max_len, [], is_directed, target_idx,
            current_length=0.0,
            last_update_time=emit_time,
            next_pivot_time=(emit_time + max(1e-12, self.params.pivot_dt_s)) if pivot_enabled else None,
            max_growth_time=emit_time + self.max_len / max(1e-12, self.params.mt_velocity_um_s),
        )
        self.mts.append(mt)
        return mt

    def emit_first_burst(self) -> None:
        if self.first_burst_done:
            return
        for _ in range(FIRST_BURST_PER_CENTROSOME):
            self.emit_visual(0, 0.0)
            self.emit_visual(1, 0.0)
        self.first_burst_done = True

    def advance(self, dt: float) -> None:
        if not self.first_burst_done:
            self.emit_first_burst()
        if self.completed() and not self.mts:
            return
        target_time = self.time + max(0.0, float(dt))
        self._advance_visual_physics_to(target_time)
        self.time = target_time
        # Update displayed capture status, KT retraction display, and archive old MTs.
        still_live: List[VisualMT] = []
        for mt in self.mts:
            age = self.time - mt.emit_time
            length = max(0.0, min(mt.max_len, age * self.params.mt_velocity_um_s))
            for kt_idx, contact, cap_time in mt.capture_contacts:
                if self.time >= cap_time:
                    kt = self.kt_units[kt_idx]
                    if not kt.displayed_captured:
                        kt.displayed_captured = True
                        kt.retracting = True
                        kt.retract_start_time = cap_time
                        kt.retract_end_time = cap_time + self.retract_duration_s
                        kt.retract_from_pos = kt.pos
                        kt.captured_by_cs_index = mt.cs_index
                        kt.retract_to_pos = self.retraction_target_pos(kt_idx, mt.cs_index)
            if age * self.params.mt_velocity_um_s > mt.max_len + 1.2:
                mt.archived = True
                self.archived_mts.append(mt)
            else:
                still_live.append(mt)
        self.mts = still_live

        # Move displayed captured KT units toward the capturing centrosome cluster.
        for kt in self.kt_units:
            if kt.displayed_captured and kt.retract_start_time is not None and kt.retract_end_time is not None:
                start = kt.retract_from_pos if kt.retract_from_pos is not None else kt.pos
                target = kt.retract_to_pos if kt.retract_to_pos is not None else kt.pos
                denom = max(1e-12, kt.retract_end_time - kt.retract_start_time)
                progress = max(0.0, min(1.0, (self.time - kt.retract_start_time) / denom))
                kt.display_pos = v_add(v_mul(start, 1.0 - progress), v_mul(target, progress))
                if progress >= 1.0:
                    kt.retracting = False
                    kt.display_pos = target

        if len(self.archived_mts) > 220:
            self.archived_mts = self.archived_mts[-220:]

# -----------------------------------------------------------------------------
# Painter widget
# -----------------------------------------------------------------------------
class SimulationView(QWidget):
    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.sim: Optional[AnimationSimulation] = None
        self.yaw = math.radians(-10.0)
        self.pitch = math.radians(78.0)
        self.zoom = 47.0
        self.last_mouse: Optional[QPointF] = None
        self.time_acceleration = 20
        self.setMinimumSize(820, 620)
        self.setMouseTracking(True)

    def set_simulation(self, sim: AnimationSimulation) -> None:
        self.sim = sim
        self.update()

    def set_time_acceleration(self, value: int) -> None:
        self.time_acceleration = max(1, min(200, int(value)))
        self.update()

    def project(self, p: Vec3) -> Tuple[QPointF, float, float]:
        x, y, z = rotate_vec(p, self.yaw, self.pitch)
        cam_dist = 46.0
        denom = max(4.0, cam_dist - z)
        s = self.zoom * cam_dist / denom
        return QPointF(self.width() / 2 + x * s, self.height() / 2 - y * s), z, s

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.last_mouse = event.position()

    def mouseReleaseEvent(self, event):
        self.last_mouse = None

    def mouseMoveEvent(self, event):
        if self.last_mouse is None:
            return
        pos = event.position()
        dx = pos.x() - self.last_mouse.x()
        dy = pos.y() - self.last_mouse.y()
        self.last_mouse = pos
        self.yaw += dx * 0.008
        self.pitch += dy * 0.008
        self.pitch = max(math.radians(-85), min(math.radians(85), self.pitch))
        self.update()

    def wheelEvent(self, event):
        self.zoom *= 1.0 + 0.0012 * event.angleDelta().y()
        self.zoom = max(12.0, min(90.0, self.zoom))
        self.update()

    def ellipsoid_point(self, theta: float, phi: float, radii: Vec3) -> Vec3:
        return (radii[0] * math.cos(theta) * math.cos(phi),
                radii[1] * math.sin(theta) * math.cos(phi),
                radii[2] * math.sin(phi))

    def draw_ellipsoid_wireframe(self, painter: QPainter, radii: Vec3) -> None:
        painter.setBrush(Qt.NoBrush)
        for it in range(9):
            theta = 2.0 * math.pi * it / 9
            path = QPainterPath()
            for j in range(100):
                phi = -math.pi / 2.0 + math.pi * j / 99.0
                q, _, _ = self.project(self.ellipsoid_point(theta, phi, radii))
                if j == 0:
                    path.moveTo(q)
                else:
                    path.lineTo(q)
            painter.setPen(QPen(COL_NUC_WIRE, 1.0))
            painter.drawPath(path)
        for ip in range(1, 7):
            phi = -math.pi / 2.0 + math.pi * ip / 7
            path = QPainterPath()
            for j in range(140):
                theta = 2.0 * math.pi * j / 139.0
                q, _, _ = self.project(self.ellipsoid_point(theta, phi, radii))
                if j == 0:
                    path.moveTo(q)
                else:
                    path.lineTo(q)
            painter.setPen(QPen(COL_NUC_WIRE, 1.0))
            painter.drawPath(path)

    def make_line_primitive(self, a: Vec3, b: Vec3, color: QColor, width: float):
        qa, za, _ = self.project(a)
        qb, zb, _ = self.project(b)
        return ((za + zb) * 0.5, "line", qa, qb, color, width)

    def make_sphere_primitive(self, center: Vec3, radius: float, color: QColor, outline: Optional[QColor] = None):
        q, z, s = self.project(center)
        r_px = max(2.0, radius * s)
        return (z, "sphere", q, r_px, color, outline)

    def current_mt_length(self, mt: VisualMT) -> float:
        if self.sim is None:
            return 0.0
        if self.sim.params.pivoting_enabled:
            return max(0.0, min(mt.max_len, mt.current_length))
        return max(0.0, min(mt.max_len, (self.sim.time - mt.emit_time) * self.sim.params.mt_velocity_um_s))

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.fillRect(self.rect(), COL_BG)
        self.draw_ellipsoid_wireframe(painter, NUC_RADII)
        if self.sim is None:
            return

        prims = []
        # Centrosome guide line
        prims.append(self.make_line_primitive(CS_A, CS_B, QColor(105, 125, 145, 80), 0.9))

        # Archived MTs are intentionally not drawn. Once a microtubule disappears,
        # it should vanish completely instead of remaining as a faint residual line.

        # Active MTs
        for mt in self.sim.mts:
            length = self.current_mt_length(mt)
            if length <= 0:
                continue
            cs = CENTROSOMES[mt.cs_index]
            end = v_add(cs, v_mul(mt.direction, length))
            done = any(self.sim.time >= cap_time for _, _, cap_time in mt.capture_contacts)
            prims.append(self.make_line_primitive(cs, end, COL_MT_DONE if done else COL_MT, 1.35 if not done else 1.8))

        # KT units. Captured units retract to the capturing centrosome and are low-alpha gray
        # so they do not obscure the centrosome cluster.
        r = self.sim.params.kt_radius_um
        for kt in self.sim.kt_units:
            col = COL_KT_DONE if kt.displayed_captured else (COL_KT if kt.is_distal else COL_KT_PROX)
            pos = kt.display_pos if kt.display_pos is not None else kt.pos
            outline = QColor(120, 130, 145, 45 if kt.displayed_captured else 90)
            prims.append(self.make_sphere_primitive(pos, r, col, outline))

        # Centrosomes
        for cs in CENTROSOMES:
            prims.append(self.make_sphere_primitive(cs, 0.42, COL_CS, QColor(70, 230, 140, 180)))
            for k in range(18):
                ang = 2 * math.pi * k / 18.0
                d = (0.0, math.cos(ang), math.sin(ang))
                prims.append(self.make_line_primitive(cs, v_add(cs, v_mul(d, 1.1)), QColor(0, 170, 85, 38), 0.55))

        prims.sort(key=lambda p: p[0])
        for prim in prims:
            if prim[1] == "line":
                _, _, qa, qb, col, width = prim
                painter.setPen(QPen(col, width, Qt.SolidLine, Qt.RoundCap))
                painter.drawLine(qa, qb)
            else:
                _, _, q, radius, col, outline = prim
                painter.setPen(QPen(outline if outline is not None else QColor(120, 130, 145, 80), 0.8))
                painter.setBrush(QBrush(col))
                painter.drawEllipse(q, radius, radius)

        # 5 μm scale bar in the lower-left corner of the 3D animation view.
        # The bar is calibrated at the central focal plane and updates with zoom.
        self.draw_scale_bar(painter, 5.0)

        # Time acceleration is a display-only quantity (not an MT kinetic parameter).
        # Its right edge is aligned with the top-right Time annotation, and its
        # baseline is aligned with the 5 μm scale-bar label.

    def draw_scale_bar(self, painter: QPainter, length_um: float = 5.0) -> None:
        if self.width() <= 0 or self.height() <= 0:
            return

        # Construct a world-space vector that maps to the camera's horizontal
        # axis without a depth change. This gives an exact screen-plane scale
        # at the central focal plane even when the scene is rotated.
        cy = math.cos(self.yaw)
        sy = math.sin(self.yaw)
        delta = (length_um * cy, -length_um * sy, 0.0)
        q0, _, _ = self.project((0.0, 0.0, 0.0))
        q1, _, _ = self.project(delta)
        bar_len_px = abs(q1.x() - q0.x())

        left_margin = 34.0
        bottom_margin = 76.0
        bar_y = self.height() - bottom_margin
        bar_x0 = left_margin
        bar_x1 = min(self.width() - left_margin, bar_x0 + bar_len_px)

        painter.setPen(QPen(QColor(250, 250, 250), 6.0, Qt.SolidLine, Qt.FlatCap))
        painter.drawLine(QPointF(bar_x0, bar_y), QPointF(bar_x1, bar_y))

        label = f"{length_um:g} μm"
        painter.setFont(QFont(DISPLAY_FONT_FAMILY, visual_font_size(20)))
        fm = painter.fontMetrics()
        label_x = bar_x0 + (bar_x1 - bar_x0 - fm.horizontalAdvance(label)) / 2.0
        label_y = bar_y - 14.0
        painter.setPen(QColor(250, 250, 250))
        painter.drawText(int(label_x), int(label_y), label)

    def draw_time_acceleration(self, painter: QPainter) -> None:
        if self.width() <= 0 or self.height() <= 0:
            return
        text = f"Time acceleration: {int(self.time_acceleration)}×"
        painter.setFont(QFont(DISPLAY_FONT_FAMILY, visual_font_size(20)))
        fm = painter.fontMetrics()
        right_margin = 24.0
        # Same baseline as the 5 μm label.
        baseline_y = self.height() - 34.0 - 14.0
        x = self.width() - right_margin - fm.horizontalAdvance(text)
        painter.setPen(QColor(250, 250, 250))
        painter.drawText(int(x), int(baseline_y), text)

# -----------------------------------------------------------------------------
# Animation status strip
# -----------------------------------------------------------------------------
class AnimationStatusWidget(QWidget):
    """Publication-style two-line status strip for the 3D animation panel.

    Row 1: MT emitted at left, Time at right.
    Row 2: Captured count/percentage directly below MT emitted.

    The font is allowed to shrink only when the available panel width requires
    it, preventing collisions in the dual-condition movie while keeping the
    enlarged publication typography whenever space permits.
    """
    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.current_time = 0.0
        self.current_mt = 0
        self.current_captured = 0
        # Two rows need more vertical room than the original one-line strip.
        self.setMinimumHeight(112)
        self.setMaximumHeight(124)

    def set_values(self, current_time: float, current_mt: int, current_captured: int) -> None:
        self.current_time = float(current_time)
        self.current_mt = int(current_mt)
        self.current_captured = int(current_captured)
        self.update()

    def _status_font(
        self, painter: QPainter, emitted_text: str, time_text: str,
        available_width: int, available_height: int
    ) -> QFont:
        """Largest enlarged Arial font that fits both width and two-row height."""
        preferred_pt = visual_font_size(23)
        # At high-resolution export this minimum is still visually large because
        # the logical compositor is rasterized at ~1.77x.
        minimum_pt = max(20, visual_font_size(11))
        gap_px = 38
        for pt in range(preferred_pt, minimum_pt - 1, -1):
            font = QFont(DISPLAY_FONT_FAMILY, pt)
            painter.setFont(font)
            fm = painter.fontMetrics()
            width_need = fm.horizontalAdvance(emitted_text) + gap_px + fm.horizontalAdvance(time_text)
            height_need = 2 * (fm.height() + 2)
            if width_need <= available_width and height_need <= available_height:
                return font
        return QFont(DISPLAY_FONT_FAMILY, minimum_pt)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.fillRect(self.rect(), COL_BG)
        painter.setPen(QColor(245, 248, 250))

        pct = self.current_captured / max(1, TOTAL_KT_UNITS) * 100.0
        emitted_text = f"MT emitted: {self.current_mt}"
        captured_text = f"Captured: {self.current_captured}/{TOTAL_KT_UNITS} ({pct:.0f}%)"
        time_text = f"Time: {self.current_time:.1f} s"

        left_margin = 20
        right_margin = 24
        available_width = max(10, self.width() - left_margin - right_margin)
        available_height = max(20, self.height() - 8)
        font = self._status_font(
            painter, emitted_text, time_text, available_width, available_height
        )
        painter.setFont(font)

        # Keep two clean rows with a small publication-style interline gap.
        fm = painter.fontMetrics()
        row_h = max(fm.height() + 2, int(available_height / 2))
        row1 = QRect(left_margin, 2, available_width, row_h)
        row2 = QRect(left_margin, 4 + row_h, available_width, row_h)

        painter.drawText(row1, Qt.AlignLeft | Qt.AlignVCenter, emitted_text)
        painter.drawText(row1, Qt.AlignRight | Qt.AlignVCenter, time_text)
        # Captured percentage intentionally stays below MT emitted.
        painter.drawText(row2, Qt.AlignLeft | Qt.AlignVCenter, captured_text)



class AnimationLegendWidget(QWidget):
    """Two-column, two-row legend shown below the animation panel.

    Left column: Centrosome, Microtubule
    Right column: Uncaptured kinetochore, Captured kinetochore
    The two columns are packed closely together, with an inter-column spacing
    approximately equal to the visual width of the letters "me".
    """
    def __init__(self, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.setMinimumHeight(96)
        self.setMaximumHeight(112)

    def _draw_circle_item(self, painter: QPainter, rect: QRect, color: QColor, label: str, text_color: QColor) -> None:
        row_center_y = rect.center().y()
        icon_r = 8.5
        icon_x = rect.left() + 4.0
        painter.setPen(Qt.NoPen)
        painter.setBrush(QBrush(color))
        painter.drawEllipse(QPointF(icon_x + icon_r, row_center_y), icon_r, icon_r)
        text_rect = QRect(int(icon_x + 26.0), rect.top(), max(10, rect.width() - int(icon_x + 26.0 - rect.left())), rect.height())
        painter.setPen(text_color)
        painter.drawText(text_rect, Qt.AlignLeft | Qt.AlignVCenter, label)

    def _draw_line_item(self, painter: QPainter, rect: QRect, color: QColor, label: str, text_color: QColor) -> None:
        row_center_y = rect.center().y()
        line_x0 = rect.left() + 2.0
        line_x1 = line_x0 + 28.0
        painter.setPen(QPen(color, 4.5, Qt.SolidLine, Qt.RoundCap))
        painter.drawLine(QPointF(line_x0, row_center_y), QPointF(line_x1, row_center_y))
        text_rect = QRect(int(line_x1 + 12.0), rect.top(), max(10, rect.width() - int(line_x1 + 12.0 - rect.left())), rect.height())
        painter.setPen(text_color)
        painter.drawText(text_rect, Qt.AlignLeft | Qt.AlignVCenter, label)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.fillRect(self.rect(), COL_BG)

        painter.setFont(QFont(DISPLAY_FONT_FAMILY, visual_font_size(15)))
        fm = painter.fontMetrics()
        text_color = QColor(245, 248, 250)

        top_margin = 2
        bottom_margin = 8
        row_gap = 4

        # Pack the two legend columns much closer together. The requested gap is
        # approximately the width of the letters "me" in the current font.
        inter_col_gap = max(10, fm.horizontalAdvance('me'))

        left_block_w = int(max(26 + fm.horizontalAdvance('Centrosome'), 42 + fm.horizontalAdvance('Microtubule')))
        right_block_w = int(max(26 + fm.horizontalAdvance('Uncaptured kinetochore'), 26 + fm.horizontalAdvance('Captured kinetochore')))
        total_block_w = left_block_w + inter_col_gap + right_block_w
        centered_start_x = int((self.width() - total_block_w) / 2.0)
        # Keep the compact v15 column spacing, but shift the whole legend group
        # leftward so it sits closer to the left side like the earlier version.
        start_x = max(8, centered_start_x - 95)

        usable_h = max(10, self.height() - top_margin - bottom_margin)
        row_h = int((usable_h - row_gap) / 2.0)
        row1_y = top_margin
        row2_y = top_margin + row_h + row_gap

        left_x = start_x
        right_x = start_x + left_block_w + inter_col_gap

        left_row1 = QRect(left_x, row1_y, left_block_w, row_h)
        left_row2 = QRect(left_x, row2_y, left_block_w, row_h)
        right_row1 = QRect(right_x, row1_y, right_block_w, row_h)
        right_row2 = QRect(right_x, row2_y, right_block_w, row_h)

        self._draw_circle_item(painter, left_row1, COL_CS, 'Centrosome', text_color)
        self._draw_line_item(painter, left_row2, QColor(0, 170, 82, 230), 'Microtubule', text_color)
        self._draw_circle_item(painter, right_row1, COL_KT, 'Uncaptured kinetochore', text_color)
        self._draw_circle_item(painter, right_row2, COL_KT_DONE_LEGEND, 'Captured kinetochore', text_color)


# -----------------------------------------------------------------------------
# Live animation plot widgets
# -----------------------------------------------------------------------------
class LivePlotWidget(QWidget):
    def __init__(self, mode: str, title: str, parent: Optional[QWidget] = None):
        super().__init__(parent)
        self.mode = mode
        self.title = title
        self.points: List[Tuple[float, float]] = []
        self.current_time = 0.0
        self.current_mt = 0
        self.current_captured = 0
        self.setMinimumSize(300, 280)

    def clear(self) -> None:
        # Keep a biologically meaningful initial point so the dynamic curves
        # start from the expected origin/state.
        self.current_time = 0.0
        self.current_mt = 0
        self.current_captured = 0
        if self.mode == "capture_time":
            self.points = [(0.0, 0.0)]
        else:
            # First burst is 23 MTs from each centrosome at t = 0.
            self.points = [(0.0, float(FIRST_BURST_PER_CENTROSOME * 2))]
        self.update()

    def update_from_sim(self, sim: Optional[AnimationSimulation]) -> None:
        if sim is None:
            self.clear()
            return
        captured = sum(1 for kt in sim.kt_units if kt.displayed_captured)
        pct = captured / max(1, len(sim.kt_units)) * 100.0
        self.current_time = sim.time
        self.current_mt = sim.mt_id_counter
        self.current_captured = captured
        if self.mode == "capture_time":
            if not self.points:
                self.points.append((0.0, 0.0))
            point = (sim.time, pct)
        else:
            if not self.points:
                self.points.append((0.0, float(FIRST_BURST_PER_CENTROSOME * 2)))
            point = (pct, float(sim.mt_id_counter))
        if not self.points or abs(point[0] - self.points[-1][0]) > 1e-9 or abs(point[1] - self.points[-1][1]) > 1e-9:
            self.points.append(point)
            if len(self.points) > 3000:
                # Preserve the initial point and trim only the dense live tail.
                self.points = [self.points[0]] + self.points[-2999:]
        self.update()

    def _nice_max(self, value: float, fallback: float) -> float:
        value = max(value, fallback)
        if value <= 0:
            return fallback
        mag = 10 ** math.floor(math.log10(value))
        norm = value / mag
        if norm <= 1:
            nice = 1
        elif norm <= 2:
            nice = 2
        elif norm <= 5:
            nice = 5
        else:
            nice = 10
        return nice * mag

    def _responsive_font_sizes(self) -> Tuple[int, int]:
        """Return tick/axis point sizes matched to the actual plot panel.

        The dual layout gives each graph much less space than the original
        single-condition view, so using a fixed global 2x font creates label
        collisions. These sizes retain the enlarged typography when there is
        room, but scale down smoothly for compact graph panels.
        """
        min_dim = max(1, min(self.width(), self.height()))
        tick_preferred = visual_font_size(15)
        axis_preferred = visual_font_size(17)
        tick_pt = max(11, min(tick_preferred, int(round(min_dim / 20.0))))
        axis_pt = max(13, min(axis_preferred, int(round(min_dim / 17.0))))
        return tick_pt, axis_pt

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing, True)
        bg = QColor(6, 13, 18)
        painter.fillRect(self.rect(), bg)

        # Thin outer panel border, consistent with a compact publication panel.
        painter.setBrush(Qt.NoBrush)
        painter.setPen(QPen(QColor(235, 240, 243), 1.4))
        painter.drawRect(1, 1, max(1, self.width() - 3), max(1, self.height() - 3))

        if self.mode == "capture_time":
            x_max = 9000.0
            y_max = 100.0
            x_label = "Time (s)"
            y_label = "Captured (%)"
        else:
            x_max = 100.0
            y_max = 3000.0
            x_label = "Captured (%)"
            y_label = "MT emitted"

        tick_pt, axis_pt = self._responsive_font_sizes()
        tick_font = QFont(DISPLAY_FONT_FAMILY, tick_pt)
        axis_font = QFont(DISPLAY_FONT_FAMILY, axis_pt)

        # Measure text first, then derive margins from the actual font metrics.
        # This is what prevents tick values and axis titles from touching.
        painter.setFont(tick_font)
        tick_fm = painter.fontMetrics()
        y_tick_labels = [f"{y_max * i / 4.0:.0f}" for i in range(5)]
        x_tick_labels = [f"{x_max * i / 4.0:.0f}" for i in range(5)]
        max_y_tick_w = max(tick_fm.horizontalAdvance(s) for s in y_tick_labels)
        end_x_tick_half = max(
            tick_fm.horizontalAdvance(x_tick_labels[0]),
            tick_fm.horizontalAdvance(x_tick_labels[-1]),
        ) / 2.0

        painter.setFont(axis_font)
        axis_fm = painter.fontMetrics()

        min_dim = max(1, min(self.width(), self.height()))
        outer_pad = max(7, int(round(min_dim * 0.022)))
        tick_len = max(5, int(round(min_dim * 0.017)))
        tick_gap = max(5, int(round(min_dim * 0.016)))
        label_gap = max(6, int(round(min_dim * 0.018)))

        left = int(outer_pad + axis_fm.height() + label_gap + max_y_tick_w + tick_gap + tick_len)
        right = int(max(outer_pad + 6, end_x_tick_half + outer_pad))
        top = int(max(outer_pad + 4, tick_fm.height() / 2.0 + outer_pad))
        bottom = int(
            outer_pad + tick_len + tick_gap + tick_fm.height() +
            label_gap + axis_fm.height()
        )

        # Safety fallback for very small GUI windows: preserve an actual plotting
        # rectangle rather than allowing the labels to consume the whole panel.
        max_left = max(42, int(self.width() * 0.38))
        max_bottom = max(54, int(self.height() * 0.38))
        left = min(left, max_left)
        bottom = min(bottom, max_bottom)
        right = min(right, max(18, int(self.width() * 0.16)))
        top = min(top, max(18, int(self.height() * 0.16)))

        w = max(20, self.width() - left - right)
        h = max(20, self.height() - top - bottom)
        x0, y0 = left, self.height() - bottom

        frame_color = QColor(238, 242, 245)
        frame_pen = QPen(frame_color, 1.8)
        painter.setPen(frame_pen)
        painter.drawLine(x0, top, x0 + w, top)
        painter.drawLine(x0 + w, top, x0 + w, y0)
        painter.drawLine(x0 + w, y0, x0, y0)
        painter.drawLine(x0, y0, x0, top)

        # Ticks and tick numbers. Rect-based placement keeps endpoint labels
        # inside their reserved margins and avoids visual collisions.
        painter.setFont(tick_font)
        tick_fm = painter.fontMetrics()
        tick_color = QColor(225, 232, 237)
        tick_h = tick_fm.height() + 4

        for i, label in enumerate(y_tick_labels):
            yy = y0 - h * i / 4.0
            painter.setPen(frame_pen)
            painter.drawLine(x0 - tick_len, yy, x0, yy)
            painter.setPen(tick_color)
            label_rect = QRect(
                outer_pad + axis_fm.height() + label_gap,
                int(yy - tick_h / 2),
                max(8, x0 - tick_len - tick_gap - (outer_pad + axis_fm.height() + label_gap)),
                tick_h,
            )
            painter.drawText(label_rect, Qt.AlignRight | Qt.AlignVCenter, label)

        x_tick_y = y0 + tick_len + tick_gap
        for i, label in enumerate(x_tick_labels):
            xx = x0 + w * i / 4.0
            painter.setPen(frame_pen)
            painter.drawLine(xx, y0, xx, y0 + tick_len)
            painter.setPen(tick_color)
            label_w = tick_fm.horizontalAdvance(label) + 12
            painter.drawText(
                QRect(int(xx - label_w / 2), int(x_tick_y), int(label_w), tick_h),
                Qt.AlignCenter,
                label,
            )

        # Axis labels use their own slightly larger responsive font and generous
        # separation from tick numbers, matching common publication figure style.
        painter.setFont(axis_font)
        axis_fm = painter.fontMetrics()
        painter.setPen(QColor(235, 240, 243))

        x_axis_y = int(x_tick_y + tick_h + label_gap - 2)
        painter.drawText(
            QRect(x0, x_axis_y, w, max(axis_fm.height() + 6, self.height() - x_axis_y - outer_pad)),
            Qt.AlignHCenter | Qt.AlignTop,
            x_label,
        )

        painter.save()
        y_axis_center_x = outer_pad + axis_fm.height() / 2.0
        y_axis_center_y = top + h / 2.0
        painter.translate(y_axis_center_x, y_axis_center_y)
        painter.rotate(-90)
        y_text_w = axis_fm.horizontalAdvance(y_label) + 14
        y_text_h = axis_fm.height() + 6
        painter.drawText(
            QRect(int(-y_text_w / 2), int(-y_text_h / 2), int(y_text_w), int(y_text_h)),
            Qt.AlignCenter,
            y_label,
        )
        painter.restore()

        # Dynamic curve. Clamp display to fixed axes so out-of-range values stay
        # inside the frame. Line width scales mildly with panel size.
        curve_width = max(1.8, min(3.0, min_dim / 120.0))
        if len(self.points) >= 2:
            painter.setPen(QPen(QColor(0, 185, 90, 235), curve_width, Qt.SolidLine, Qt.RoundCap))
            prev = None
            for x, y in self.points:
                xc = max(0.0, min(float(x), x_max))
                yc = max(0.0, min(float(y), y_max))
                px = x0 + (xc / max(1e-12, x_max)) * w
                py = y0 - (yc / max(1e-12, y_max)) * h
                if prev is not None:
                    painter.drawLine(prev, QPointF(px, py))
                prev = QPointF(px, py)
        elif len(self.points) == 1:
            x, y = self.points[0]
            xc = max(0.0, min(float(x), x_max))
            yc = max(0.0, min(float(y), y_max))
            px = x0 + (xc / max(1e-12, x_max)) * w
            py = y0 - (yc / max(1e-12, y_max)) * h
            painter.setBrush(QBrush(QColor(0, 185, 90, 235)))
            painter.setPen(Qt.NoPen)
            radius = max(2.8, min(4.5, min_dim / 80.0))
            painter.drawEllipse(QPointF(px, py), radius, radius)

class DualComparisonExportWidget(QWidget):
    """Off-screen Movie S10 compositor used for publication-resolution MP4 output.

    Logical geometry mirrors the two original 1584-px-wide source panels placed
    side by side. Because this widget is independent of the visible GUI, export
    quality and aspect ratio do not depend on monitor resolution or window size.
    """
    LOGICAL_WIDTH = 3300
    LOGICAL_HEIGHT = 1280

    def __init__(
        self,
        conditions: Sequence[ConditionSpec],
        sims: Sequence[AnimationSimulation],
        source_views: Optional[Sequence[SimulationView]] = None,
        parent: Optional[QWidget] = None,
    ):
        super().__init__(parent)
        self.conditions = list(conditions)
        self.sims = list(sims)
        self.setFixedSize(self.LOGICAL_WIDTH, self.LOGICAL_HEIGHT)
        self.setStyleSheet("QWidget {background-color: rgb(0, 0, 0); border: none;}")

        self.views: List[SimulationView] = [SimulationView(), SimulationView()]
        self.capture_time_plots: List[LivePlotWidget] = [
            LivePlotWidget("capture_time", "KT captured (%) vs time"),
            LivePlotWidget("capture_time", "KT captured (%) vs time"),
        ]
        self.mt_capture_plots: List[LivePlotWidget] = [
            LivePlotWidget("mt_capture", "MT emitted vs KT captured (%)"),
            LivePlotWidget("mt_capture", "MT emitted vs KT captured (%)"),
        ]
        self.status_labels: List[AnimationStatusWidget] = [AnimationStatusWidget(), AnimationStatusWidget()]
        self.legends: List[AnimationLegendWidget] = [AnimationLegendWidget(), AnimationLegendWidget()]

        outer = QVBoxLayout(self)
        outer.setSpacing(0)
        outer.setContentsMargins(0, 0, 0, 0)

        movie_title = QLabel("Movie S10: Three-dimensional simulation of random search-and-capture")
        movie_title.setAlignment(Qt.AlignCenter)
        movie_title.setFont(QFont(DISPLAY_FONT_FAMILY, visual_font_size(24)))
        movie_title.setStyleSheet("QLabel {color: rgb(248,248,248); background: black; border: none;}")
        movie_title.setFixedHeight(94)
        outer.addWidget(movie_title)
        outer.addSpacing(96)

        pair_row = QHBoxLayout()
        pair_row.setSpacing(72)
        pair_row.setContentsMargins(0, 0, 0, 0)
        outer.addLayout(pair_row, 1)

        for i in range(2):
            self.views[i].setMinimumSize(900, 760)
            self.capture_time_plots[i].setMinimumSize(470, 360)
            self.mt_capture_plots[i].setMinimumSize(470, 360)
            self.views[i].set_simulation(self.sims[i])
            if source_views is not None and i < len(source_views):
                self.views[i].yaw = source_views[i].yaw
                self.views[i].pitch = source_views[i].pitch
                self.views[i].zoom = source_views[i].zoom

            plot_col = QVBoxLayout()
            plot_col.setSpacing(8)
            plot_col.setContentsMargins(0, 0, 0, 0)
            plot_col.addWidget(self.capture_time_plots[i], 1)
            plot_col.addWidget(self.mt_capture_plots[i], 1)

            animation_inner = QVBoxLayout()
            animation_inner.setSpacing(0)
            animation_inner.setContentsMargins(0, 0, 0, 4)
            animation_inner.addWidget(self.status_labels[i])
            animation_inner.addWidget(self.views[i], 1)
            animation_inner.addWidget(self.legends[i])

            animation_panel = QFrame()
            animation_panel.setStyleSheet(
                "QFrame {background-color: rgb(3, 9, 12); border: 2px solid rgb(235, 240, 243);}"
            )
            animation_panel_layout = QVBoxLayout(animation_panel)
            animation_panel_layout.setSpacing(0)
            animation_panel_layout.setContentsMargins(7, 7, 7, 7)
            animation_panel_layout.addLayout(animation_inner, 1)

            module_visual = QFrame()
            module_visual.setStyleSheet("QFrame {background-color: rgb(0,0,0); border: none;}")
            module_visual_layout = QHBoxLayout(module_visual)
            module_visual_layout.setSpacing(8)
            module_visual_layout.setContentsMargins(0, 0, 0, 0)
            module_visual_layout.addWidget(animation_panel, 2)
            module_visual_layout.addLayout(plot_col, 1)

            name = str(self.conditions[i].name).replace(" directed", " Directed")
            condition_title = QLabel(name)
            condition_title.setAlignment(Qt.AlignCenter)
            condition_title.setFont(QFont(DISPLAY_FONT_FAMILY, visual_font_size(20)))
            condition_title.setStyleSheet("QLabel {color: rgb(248,248,248); background: black; border: none;}")
            condition_title.setFixedHeight(76)

            module = QFrame()
            module.setStyleSheet("QFrame {background-color: rgb(0,0,0); border: none;}")
            module_col = QVBoxLayout(module)
            module_col.setSpacing(8)
            module_col.setContentsMargins(0, 0, 0, 0)
            module_col.addWidget(condition_title)
            module_col.addWidget(module_visual, 1)
            pair_row.addWidget(module, 1)

        self.layout().activate()
        self.update_frame(20)

    def update_frame(self, acceleration: int) -> None:
        acceleration = max(1, min(200, int(acceleration)))
        for i, sim in enumerate(self.sims):
            self.views[i].set_time_acceleration(acceleration)
            captured = sum(1 for kt in sim.kt_units if kt.displayed_captured)
            self.status_labels[i].set_values(sim.time, sim.mt_id_counter, captured)
            self.capture_time_plots[i].update_from_sim(sim)
            self.mt_capture_plots[i].update_from_sim(sim)
            self.views[i].update()
        self.update()


# -----------------------------------------------------------------------------
# Export helpers
# -----------------------------------------------------------------------------
def write_csv(path: str, rows: List[Dict[str, object]], fieldnames: Optional[List[str]] = None) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if fieldnames is None:
        keys: List[str] = []
        for r in rows:
            for k in r.keys():
                if k not in keys:
                    keys.append(k)
        fieldnames = keys
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)



def _safe_float(value) -> Optional[float]:
    try:
        if value == "" or value is None:
            return None
        v = float(value)
        if math.isnan(v):
            return None
        return v
    except Exception:
        return None


def _mean_safe(values: Sequence[object]) -> float:
    nums = [v for v in (_safe_float(x) for x in values) if v is not None]
    return float(statistics.fmean(nums)) if nums else float("nan")


def _sanitize_filename(text: str) -> str:
    out = []
    for ch in str(text):
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        elif ch in ("%",):
            out.append("pct")
        else:
            out.append("_")
    name = "".join(out).strip("_")
    while "__" in name:
        name = name.replace("__", "_")
    return name or "condition"


def write_matrix_csv(path: str, rows: Sequence[Sequence[object]]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerows(rows)


def _wide_columns_by_condition(condition_values: Dict[str, List[object]], condition_order: Sequence[str]) -> List[List[object]]:
    ordered = [c for c in condition_order if c in condition_values]
    max_len = max((len(condition_values.get(c, [])) for c in ordered), default=0)
    rows: List[List[object]] = [list(ordered)]
    for i in range(max_len):
        row = []
        for cond in ordered:
            vals = condition_values.get(cond, [])
            row.append(vals[i] if i < len(vals) else "")
        rows.append(row)
    return rows


def _format_fraction(index: int, total: int = TOTAL_KT_UNITS) -> str:
    return f"{index}/{total}"


def generate_superplot_ready_tables(
    output_dir: str,
    condition_order: Sequence[str],
    summary_rows: Sequence[Dict[str, object]],
    progression_rows: Sequence[Dict[str, object]],
    deflection_rows: Sequence[Dict[str, object]],
) -> List[str]:
    """Create Prism/superplot-ready tables from raw export rows.

    Files are written under Superplot_ready_tables/.
    - Raw metric tables keep run/cell records vertically and conditions horizontally.
    - Capture progression keeps condition blocks horizontally, with all cells stacked vertically.
    - Average tables summarize across all cells/encounters for quick review.
    Returns a lightweight QC log as a list of text lines.
    """
    super_dir = os.path.join(output_dir, "Superplot_ready_tables")
    os.makedirs(super_dir, exist_ok=True)
    conditions = [c for c in condition_order if any(str(r.get("condition")) == c for r in summary_rows)]
    qc_lines: List[str] = []
    qc_lines.append("Superplot-ready table QC")
    qc_lines.append(f"Conditions: {', '.join(conditions)}")
    run_ids = sorted({str(r.get("run", "1")) for r in summary_rows})
    qc_lines.append(f"Independent runs represented: {', '.join(run_ids)}")

    # 1. Completion time and emitted MT count, one value per repeat/cell.
    completion_by_cond: Dict[str, List[object]] = {c: [] for c in conditions}
    mt_by_cond: Dict[str, List[object]] = {c: [] for c in conditions}
    for cond in conditions:
        rows = sorted([r for r in summary_rows if str(r.get("condition")) == cond], key=lambda r: int(r.get("repeat", 0)))
        for r in rows:
            completion_by_cond[cond].append(r.get("completion_time_s", ""))
            mt_by_cond[cond].append(r.get("mt_total_emitted", ""))
        qc_lines.append(f"{cond}: summary repeats = {len(rows)}")
    write_matrix_csv(os.path.join(super_dir, "Superplot_completion_time_s.csv"), _wide_columns_by_condition(completion_by_cond, conditions))
    write_matrix_csv(os.path.join(super_dir, "Superplot_MT_total_emitted.csv"), _wide_columns_by_condition(mt_by_cond, conditions))

    # 2. Deflection angles, all encounters stacked vertically per condition.
    defl_by_cond: Dict[str, List[object]] = {c: [] for c in conditions}
    for cond in conditions:
        rows = [r for r in deflection_rows if str(r.get("condition")) == cond]
        for r in rows:
            defl_by_cond[cond].append(r.get("deflection_to_boundary_deg", ""))
        qc_lines.append(f"{cond}: deflection encounters = {len(rows)}")
    write_matrix_csv(os.path.join(super_dir, "Superplot_deflection_angle_deg.csv"), _wide_columns_by_condition(defl_by_cond, conditions))

    # 3. Capture progression by cell: each condition occupies 3 columns; cells stack downward.
    #    Every completed cell should have exactly 46 capture records.
    # Normalize progression chronologically. During the first burst, multiple MTs
    # have emit_time=0 and the event loop processes them sequentially; the raw
    # capture_index is therefore a processing rank, not always a chronological
    # rank. For "time required to reach n/46", we sort by capture_time_s and
    # re-rank captures from 1/46 to 46/46.
    raw_prog_by_cond_rep: Dict[Tuple[str, int], List[Dict[str, object]]] = {}
    repeat_set = set()
    for r in progression_rows:
        cond = str(r.get("condition"))
        if cond not in conditions:
            continue
        try:
            rep = int(r.get("repeat", 0))
        except Exception:
            continue
        repeat_set.add(rep)
        raw_prog_by_cond_rep.setdefault((cond, rep), []).append(r)
    repeats = sorted(repeat_set)

    prog_by_cond_rep: Dict[Tuple[str, int], Dict[int, Dict[str, object]]] = {}
    for key, rows_for_rep in raw_prog_by_cond_rep.items():
        sorted_rows = sorted(
            rows_for_rep,
            key=lambda r: (
                _safe_float(r.get("capture_time_s")) if _safe_float(r.get("capture_time_s")) is not None else float("inf"),
                _safe_float(r.get("mt_total_emitted_at_capture")) if _safe_float(r.get("mt_total_emitted_at_capture")) is not None else float("inf"),
                _safe_float(r.get("capture_index")) if _safe_float(r.get("capture_index")) is not None else float("inf"),
            ),
        )
        ranked: Dict[int, Dict[str, object]] = {}
        for rank, rec in enumerate(sorted_rows[:TOTAL_KT_UNITS], start=1):
            rec2 = dict(rec)
            rec2["chronological_capture_index"] = rank
            ranked[rank] = rec2
        prog_by_cond_rep[key] = ranked

    header1: List[object] = []
    header2: List[object] = []
    for cond in conditions:
        header1.extend([cond, "", ""])
        header2.extend(["Captured fraction", "Time required (s)", "MT emitted"])
    progression_matrix: List[List[object]] = [header1, header2]

    missing_progression_warnings = 0
    time_monotonic_warnings = 0
    for rep in repeats:
        last_times: Dict[str, Optional[float]] = {c: None for c in conditions}
        last_mts: Dict[str, Optional[float]] = {c: None for c in conditions}
        for idx in range(1, TOTAL_KT_UNITS + 1):
            row_out: List[object] = []
            for cond in conditions:
                rec = prog_by_cond_rep.get((cond, rep), {}).get(idx)
                if rec is None:
                    missing_progression_warnings += 1
                    row_out.extend([_format_fraction(idx), "", ""])
                    continue
                t = rec.get("capture_time_s", "")
                mt = rec.get("mt_total_emitted_at_capture", "")
                t_num = _safe_float(t)
                mt_num = _safe_float(mt)
                if last_times[cond] is not None and t_num is not None and t_num + 1e-9 < float(last_times[cond]):
                    time_monotonic_warnings += 1
                if t_num is not None:
                    last_times[cond] = t_num
                if mt_num is not None:
                    last_mts[cond] = mt_num
                row_out.extend([_format_fraction(idx), t, mt])
            progression_matrix.append(row_out)
        # Blank line between repeats for readability.
        progression_matrix.append(["" for _ in header2])
    write_matrix_csv(os.path.join(super_dir, "Superplot_capture_progression_by_repeat.csv"), progression_matrix)
    qc_lines.append(f"Capture progression repeats = {len(repeats)}")
    qc_lines.append(f"Capture progression missing cells = {missing_progression_warnings}")
    qc_lines.append(f"Capture progression time-order warnings = {time_monotonic_warnings}")

    # 4. Average summary metrics: one compact table for the three core outputs.
    avg_summary: List[List[object]] = [["Metric"] + conditions]
    avg_summary.append(["Mean completion time (s)"] + [f"{_mean_safe(completion_by_cond.get(c, [])):.6g}" for c in conditions])
    avg_summary.append(["Mean MT emitted"] + [f"{_mean_safe(mt_by_cond.get(c, [])):.6g}" for c in conditions])
    avg_summary.append(["Mean deflection angle (deg)"] + [f"{_mean_safe(defl_by_cond.get(c, [])):.6g}" for c in conditions])
    write_matrix_csv(os.path.join(super_dir, "Average_summary_metrics.csv"), avg_summary)

    # 5. Average capture progression: each condition occupies 3 columns, 46 rows total.
    avg_header1: List[object] = []
    avg_header2: List[object] = []
    for cond in conditions:
        avg_header1.extend([cond, "", ""])
        avg_header2.extend(["Captured fraction", "Mean time required (s)", "Mean MT emitted"])
    avg_prog: List[List[object]] = [avg_header1, avg_header2]
    for idx in range(1, TOTAL_KT_UNITS + 1):
        row_out: List[object] = []
        for cond in conditions:
            time_values = []
            mt_values = []
            for rep in repeats:
                rec = prog_by_cond_rep.get((cond, rep), {}).get(idx)
                if rec is None:
                    continue
                time_values.append(rec.get("capture_time_s", ""))
                mt_values.append(rec.get("mt_total_emitted_at_capture", ""))
            row_out.extend([
                _format_fraction(idx),
                f"{_mean_safe(time_values):.6g}",
                f"{_mean_safe(mt_values):.6g}",
            ])
        avg_prog.append(row_out)
    write_matrix_csv(os.path.join(super_dir, "Average_capture_progression.csv"), avg_prog)

    with open(os.path.join(super_dir, "Superplot_QC_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(qc_lines) + "\n")
    return qc_lines



def generate_time_superplot_ready_tables(
    output_dir: str,
    condition_order: Sequence[str],
    summary_rows: Sequence[Dict[str, object]],
) -> List[str]:
    """Create Prism-ready completion-time tables only.

    This is used when the user checks only the dedicated Time output.
    It avoids writing MT-number, deflection, or capture-progression
    superplot tables, so the output folder stays clean.
    """
    super_dir = os.path.join(output_dir, "Superplot_ready_tables")
    os.makedirs(super_dir, exist_ok=True)
    conditions = [c for c in condition_order if any(str(r.get("condition")) == c for r in summary_rows)]

    qc_lines: List[str] = []
    qc_lines.append("Time-only superplot table QC")
    qc_lines.append(f"Conditions: {', '.join(conditions)}")

    completion_by_cond: Dict[str, List[object]] = {c: [] for c in conditions}
    for cond in conditions:
        rows = sorted([r for r in summary_rows if str(r.get("condition")) == cond], key=lambda r: int(r.get("repeat", 0)))
        for r in rows:
            completion_by_cond[cond].append(r.get("completion_time_s", ""))
        completed_values = [v for v in completion_by_cond[cond] if _safe_float(v) is not None]
        qc_lines.append(f"{cond}: repeats = {len(rows)}, completed time values = {len(completed_values)}")

    write_matrix_csv(os.path.join(super_dir, "Superplot_completion_time_s.csv"), _wide_columns_by_condition(completion_by_cond, conditions))

    avg_rows: List[List[object]] = [["Metric"] + conditions]
    avg_rows.append(["Mean completion time (s)"] + [f"{_mean_safe(completion_by_cond.get(c, [])):.6g}" for c in conditions])
    write_matrix_csv(os.path.join(super_dir, "Average_completion_time_s.csv"), avg_rows)

    with open(os.path.join(super_dir, "Time_only_QC_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(qc_lines) + "\n")
    return qc_lines

def generate_metric_only_superplot_ready_tables(
    output_dir: str,
    condition_order: Sequence[str],
    summary_rows: Sequence[Dict[str, object]],
    include_time: bool,
    include_mt_total: bool,
) -> List[str]:
    """Create Prism-ready metric tables only.

    This is used when the user checks only dedicated scalar outputs such as
    Completion time and/or MT total emitted. It avoids writing deflection,
    capture-progression, or per-MT superplot tables, so the folder stays clean.
    """
    super_dir = os.path.join(output_dir, "Superplot_ready_tables")
    os.makedirs(super_dir, exist_ok=True)
    conditions = [c for c in condition_order if any(str(r.get("condition")) == c for r in summary_rows)]

    qc_lines: List[str] = []
    active_metrics = []
    if include_time:
        active_metrics.append("completion_time_s")
    if include_mt_total:
        active_metrics.append("mt_total_emitted")
    qc_lines.append("Metric-only superplot table QC")
    qc_lines.append(f"Metrics: {', '.join(active_metrics)}")
    qc_lines.append(f"Conditions: {', '.join(conditions)}")

    completion_by_cond: Dict[str, List[object]] = {c: [] for c in conditions}
    mt_by_cond: Dict[str, List[object]] = {c: [] for c in conditions}
    for cond in conditions:
        rows = sorted([r for r in summary_rows if str(r.get("condition")) == cond], key=lambda r: int(r.get("repeat", 0)))
        for r in rows:
            completion_by_cond[cond].append(r.get("completion_time_s", ""))
            mt_by_cond[cond].append(r.get("mt_total_emitted", ""))
        completed_times = [v for v in completion_by_cond[cond] if _safe_float(v) is not None]
        completed_mts = [v for v in mt_by_cond[cond] if _safe_float(v) is not None]
        qc_lines.append(f"{cond}: repeats = {len(rows)}, time values = {len(completed_times)}, MT-total values = {len(completed_mts)}")

    avg_rows: List[List[object]] = [["Metric"] + conditions]
    if include_time:
        write_matrix_csv(os.path.join(super_dir, "Superplot_completion_time_s.csv"), _wide_columns_by_condition(completion_by_cond, conditions))
        avg_rows.append(["Mean completion time (s)"] + [f"{_mean_safe(completion_by_cond.get(c, [])):.6g}" for c in conditions])
    if include_mt_total:
        write_matrix_csv(os.path.join(super_dir, "Superplot_MT_total_emitted.csv"), _wide_columns_by_condition(mt_by_cond, conditions))
        avg_rows.append(["Mean MT total emitted"] + [f"{_mean_safe(mt_by_cond.get(c, [])):.6g}" for c in conditions])

    if len(avg_rows) > 1:
        write_matrix_csv(os.path.join(super_dir, "Average_selected_metrics.csv"), avg_rows)
    if include_time and not include_mt_total:
        write_matrix_csv(os.path.join(super_dir, "Average_completion_time_s.csv"), [["Metric"] + conditions, avg_rows[1]])
    if include_mt_total and not include_time:
        write_matrix_csv(os.path.join(super_dir, "Average_MT_total_emitted.csv"), [["Metric"] + conditions, avg_rows[1]])

    with open(os.path.join(super_dir, "Metric_only_QC_report.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(qc_lines) + "\n")
    return qc_lines


# -----------------------------------------------------------------------------
# BOSS-style direct output helpers
# -----------------------------------------------------------------------------
def _mean(values: Sequence[float]) -> float:
    return float(statistics.fmean(values)) if values else float("nan")


def _median(values: Sequence[float]) -> float:
    return float(statistics.median(values)) if values else float("nan")


def _pstdev(values: Sequence[float]) -> float:
    return float(statistics.pstdev(values)) if len(values) >= 2 else 0.0


def _condition_label_from_prob(prob: float) -> str:
    if abs(prob) <= 1e-12:
        return "Random"
    return f"{int(round(prob * 100))}% directed"


def _print_angle_summary(title: str, angles: List[float]) -> None:
    print(title)
    if not angles:
        print("  n = 0")
        print("  no angle data")
        return
    zeros = [a for a in angles if abs(a) <= 1e-12]
    lt5 = [a for a in angles if a < 5.0]
    lt12 = [a for a in angles if a < 12.0]
    gt12 = [a for a in angles if a > 12.0]
    print(f"  n = {len(angles)}")
    print(f"  平均 {_mean(angles):.2f}, 中位数 {_median(angles):.2f} (度)")
    print(f"  theta=0 比例: {len(zeros) / len(angles) * 100:.2f}%")
    print(f"  <5deg: {len(lt5) / len(angles) * 100:.1f}% | <12deg: {len(lt12) / len(angles) * 100:.1f}% | >12deg: {len(gt12) / len(angles) * 100:.1f}%")


def run_boss_style_direct_output(
    params: SimulationParams,
    directed_prob: float = DIRECT_F_DIRECTED,
    n_cells: int = DIRECT_N_CELLS,
    output_dir: str = DIRECT_OUTPUT_DIR,
    progress_milestones: Sequence[float] = PROGRESS_MILESTONES,
    condition: Optional[ConditionSpec] = None,
    n_runs: int = 1,
) -> Dict[str, object]:
    """Run without GUI and export the same core files as the BOSS script.

    Output files:
      - cell_XXXX_progression.csv
      - aggregate_progression.csv
      - per_cell_summary.csv
      - run_summary.csv
      - distal_angles.csv

    Console output mirrors the BOSS summary: completion count, MT count/time,
    distal/all deflection statistics, and milestone time/MT summaries.
    """
    os.makedirs(output_dir, exist_ok=True)
    n_cells = max(1, int(n_cells))
    n_runs = max(1, int(n_runs))
    total_cells = n_cells * n_runs
    if condition is None:
        condition = ConditionSpec(_condition_label_from_prob(directed_prob), directed_prob)
    sim_params = params_for_condition(params, condition)
    condition_label = condition.name

    print("=" * 70)
    print("3D Random Search-and-Capture Simulation (Hybrid KT-unit direct mode)")
    print("=" * 70)
    print(f"细胞核: {NUC_LONG_AXIS:.2f} x {NUC_INTERMEDIATE_AXIS:.2f} x {NUC_SHORT_AXIS:.2f} um")
    print(f"中心体间距: {CENTROSOME_DISTANCE:.2f} um")
    print(f"KT units: {TOTAL_KT_UNITS} ({DISTAL_KT_UNITS} distal, {PROXIMAL_KT_UNITS} proximal)")
    print(f"KT generation: {sim_params.kt_generation_mode}")
    print(
        f"KT base diameter: {sim_params.kt_base_diameter_um:.3f} um | "
        f"volume multiplier: {sim_params.kt_volume_multiplier:.2f}x | "
        f"effective KT diameter: {sim_params.kt_diameter_um:.3f} um | "
        f"effective KT volume: {sim_params.kt_effective_volume_um3:.3f} um^3"
    )
    print(f"MT 直径: {sim_params.mt_diameter_um * 1000:.0f} nm | MT 速度: {sim_params.mt_velocity_um_s:.3f} um/s")
    print(f"发出: burst {FIRST_BURST_PER_CENTROSOME}/中心体, 之后 {sim_params.emission_rate_per_centrosome:.4f} 根/s/中心体")
    print(f"定向程度 f: {condition.directed_prob * 100:.0f}% ({'纯random' if condition.directed_prob == 0 else '全部定向' if condition.directed_prob == 1 else '混合'})")
    print(f"条件组: {condition.family} | pivoting: {condition_metadata(condition)['pivoting_parameter']}")
    print(f"独立运行次数: {n_runs} | 每次运行细胞数: {n_cells} | 总细胞数: {total_cells}")
    print("=" * 70)
    print()

    all_dist_angles: List[float] = []
    all_angles: List[float] = []
    mt_counts: List[int] = []
    times: List[float] = []
    completed_count = 0
    milestone_times: Dict[float, List[float]] = {m: [] for m in progress_milestones}
    milestone_mts: Dict[float, List[int]] = {m: [] for m in progress_milestones}
    per_cell_summary: List[Dict[str, object]] = []

    for run_index in range(1, n_runs + 1):
        for c in range(n_cells):
            cell_index = c + 1
            global_cell_index = (run_index - 1) * n_cells + cell_index
            layout_seed = params.seed + global_cell_index * 100003
            sim_seed = params.seed + global_cell_index * 100003 + 17

            try:
                kt_layout = generate_kt_layout(sim_params, layout_seed)
            except Exception as e:
                print(f"  Run {run_index}, cell {cell_index}: KT 放置失败,跳过 ({e})")
                continue

            sim = EventDrivenSimulator(sim_params, condition.directed_prob, sim_seed, kt_layout)
            result = sim.run()
            completed = bool(result["completed"])
            n_mts = int(result["n_mts"])
            n_captured = int(result["n_captured"])
            final_time = float(result["final_time_s"])

            if completed:
                completed_count += 1
                mt_counts.append(n_mts)
                times.append(final_time)

            for row in result["deflection_rows"]:
                try:
                    ang = float(row["deflection_to_boundary_deg"])
                except Exception:
                    continue
                all_angles.append(ang)
                if bool(row.get("nearest_is_distal", False)):
                    all_dist_angles.append(ang)

            capture_log = list(result["capture_log"])
            cell_file = os.path.join(output_dir, f"cell_{global_cell_index:04d}_progression.csv")
            with open(cell_file, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow([
                    "run", "cell", "global_cell", "capture_index",
                    "fraction_captured", "time_s", "n_mts_emitted",
                    "condition", "condition_family", "pivoting_mode",
                    "d_mt_rad2_s",
                ])
                for row in capture_log:
                    capture_index = int(row["capture_index"])
                    writer.writerow([
                        run_index, cell_index, global_cell_index, capture_index,
                        f"{capture_index / TOTAL_KT_UNITS:.4f}",
                        f"{float(row['capture_time_s']):.2f}",
                        int(row["mt_total_emitted_at_capture"]),
                        condition.name,
                        condition.family,
                        condition.pivoting_mode,
                        condition.d_mt_rad2_s,
                    ])

            for m in progress_milestones:
                target = int(math.ceil(m * TOTAL_KT_UNITS))
                rec = next((row for row in capture_log if int(row["capture_index"]) >= target), None)
                if rec is not None:
                    milestone_times[m].append(float(rec["capture_time_s"]))
                    milestone_mts[m].append(int(rec["mt_total_emitted_at_capture"]))

            per_cell_summary.append({
                "run": run_index,
                "cell": cell_index,
                "global_cell": global_cell_index,
                "completed": completed,
                "n_mts": n_mts,
                "n_captured": n_captured,
                "final_time": final_time,
                "layout_seed": layout_seed,
                "simulation_seed": sim_seed,
                **condition_metadata(condition),
            })

            if cell_index % 10 == 0:
                print(f"  Run {run_index}/{n_runs}: 已完成 {cell_index}/{n_cells} 个细胞...")

    print()
    print("=" * 70)
    print("统计结果 SUMMARY")
    print("=" * 70)
    print(f"完成捕捉的细胞: {completed_count}/{total_cells}")
    print()
    if mt_counts:
        print(f"全捕捉所需微管数: 平均 {_mean(mt_counts):.0f}, 中位数 {_median(mt_counts):.0f}, 标准差 {_pstdev(mt_counts):.0f}")
        print(f"全捕捉所需时间:   平均 {_mean(times):.0f} s ({_mean(times)/60:.1f} min), 中位数 {_median(times)/60:.1f} min")
    print()
    _print_angle_summary("--- 偏转角 (Distal KT encounters) ---", all_dist_angles)
    print()
    _print_angle_summary("--- 偏转角 (All KT encounters) ---", all_angles)
    print()
    print("--- 捕捉进度 / 时间 ---")
    for m in progress_milestones:
        if milestone_times[m]:
            mt_t = _mean(milestone_times[m])
            print(f"  {int(m*100):>3}% ({int(math.ceil(m*TOTAL_KT_UNITS)):>2}/{TOTAL_KT_UNITS}): 平均 {mt_t:>8.1f} s = {mt_t/60:>6.2f} min")
    print()
    print("--- 捕捉进度 / 微管数 ---")
    for m in progress_milestones:
        if milestone_mts[m]:
            mt_n = _mean(milestone_mts[m])
            print(f"  {int(m*100):>3}% ({int(math.ceil(m*TOTAL_KT_UNITS)):>2}/{TOTAL_KT_UNITS}): 平均 {mt_n:>8.0f} 根微管")
    print()

    run_summary_rows: List[Dict[str, object]] = []
    for run_index in range(1, n_runs + 1):
        run_rows = [row for row in per_cell_summary if int(row["run"]) == run_index]
        run_completed = [row for row in run_rows if bool(row["completed"])]
        run_summary_rows.append({
            "run": run_index,
            "cells_requested": n_cells,
            "cells_simulated": len(run_rows),
            "completed_cells": sum(1 for row in run_completed if bool(row["completed"])),
            "mean_completion_time_s": _mean([float(row["final_time"]) for row in run_completed]),
            "mean_mt_total_emitted": _mean([int(row["n_mts"]) for row in run_completed]),
            **condition_metadata(condition),
        })
    write_csv(os.path.join(output_dir, "run_summary.csv"), run_summary_rows)

    # Boss-style aggregate progression output.
    with open(os.path.join(output_dir, "aggregate_progression.csv"), "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow([
            "milestone_fraction", "kt_captured", "mean_time_s", "median_time_s",
            "mean_n_mts", "median_n_mts", "condition_family", "pivoting_mode",
            "d_mt_rad2_s",
        ])
        for m in progress_milestones:
            if milestone_times[m]:
                writer.writerow([
                    m,
                    int(math.ceil(m * TOTAL_KT_UNITS)),
                    f"{_mean(milestone_times[m]):.2f}",
                    f"{_median(milestone_times[m]):.2f}",
                    f"{_mean(milestone_mts[m]):.1f}",
                    f"{_median(milestone_mts[m]):.1f}",
                    condition.family,
                    condition.pivoting_mode,
                    condition.d_mt_rad2_s,
                ])

    with open(os.path.join(output_dir, "per_cell_summary.csv"), "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow([
            "run", "cell", "global_cell", "completed", "n_mts_total",
            "final_time_s", "final_time_min",
            "condition", "condition_family", "pivoting_mode", "d_mt_rad2_s",
        ])
        for s in per_cell_summary:
            writer.writerow([
                s["run"],
                s["cell"],
                s["global_cell"],
                s["completed"],
                s["n_mts"],
                f"{float(s['final_time']):.2f}",
                f"{float(s['final_time'])/60:.2f}",
                s["condition"],
                s["condition_family"],
                s["pivoting_mode"],
                s["d_mt_rad2_s"],
            ])

    with open(os.path.join(output_dir, "distal_angles.csv"), "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow([
            "deflection_angle_deg", "condition_family", "pivoting_mode",
            "d_mt_rad2_s",
        ])
        for a in all_dist_angles:
            writer.writerow([
                f"{a:.3f}", condition.family, condition.pivoting_mode,
                condition.d_mt_rad2_s,
            ])

    # Extra reproducibility file: records the exact GUI/CLI parameters used for
    # this direct BOSS-style calculation. The four core CSVs above keep the same
    # structure as the original BOSS script; this file makes parameter tracking explicit.
    with open(os.path.join(output_dir, "parameters_used.csv"), "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["parameter", "value"])
        writer.writerow(["condition", condition_label])
        writer.writerow(["condition_family", condition.family])
        writer.writerow(["directed_fraction", f"{condition.directed_prob:.6g}"])
        writer.writerow(["directed_percent", f"{condition.directed_prob * 100:.3f}"])
        writer.writerow(["pivoting_mode", condition.pivoting_mode])
        writer.writerow(["pivoting_parameter", condition_metadata(condition)["pivoting_parameter"]])
        writer.writerow(["d_mt_rad2_s", f"{condition.d_mt_rad2_s:.9g}"])
        writer.writerow(["pivot_dt_s", f"{PIVOT_DT_S:.6g}"])
        writer.writerow(["theta_sample_dt_s", f"{THETA_SAMPLE_DT_S:.6g}"])
        writer.writerow(["kt_generation_mode", sim_params.kt_generation_mode])
        writer.writerow(["kt_base_diameter_um", f"{sim_params.kt_base_diameter_um:.3f}"])
        writer.writerow(["kt_volume_multiplier", f"{sim_params.kt_volume_multiplier:.2f}"])
        writer.writerow(["kt_base_volume_um3", f"{sim_params.kt_base_volume_um3:.6f}"])
        writer.writerow(["kt_unit_diameter_um", f"{sim_params.kt_diameter_um:.3f}"])
        writer.writerow(["kt_effective_volume_um3", f"{sim_params.kt_effective_volume_um3:.6f}"])
        writer.writerow(["mt_diameter_um", f"{sim_params.mt_diameter_um:.6f}"])
        writer.writerow(["mt_velocity_um_s", f"{sim_params.mt_velocity_um_s:.6f}"])
        writer.writerow(["post_burst_emission_rate_per_centrosome_MT_per_s", f"{sim_params.emission_rate_per_centrosome:.6f}"])
        writer.writerow(["mt_extra_length_um", f"{sim_params.mt_extra_length_um:.6f}"])
        writer.writerow(["number_of_runs", n_runs])
        writer.writerow(["cells_per_run", n_cells])
        writer.writerow(["total_cells", total_cells])
        writer.writerow(["random_seed", params.seed])
        writer.writerow(["total_kt_units", TOTAL_KT_UNITS])
        writer.writerow(["distal_kt_units", DISTAL_KT_UNITS])
        writer.writerow(["proximal_kt_units", PROXIMAL_KT_UNITS])
        writer.writerow(["first_burst_per_centrosome", FIRST_BURST_PER_CENTROSOME])
        writer.writerow(["first_burst_total", FIRST_BURST_PER_CENTROSOME * 2])
        writer.writerow(["nucleus_long_axis_um", f"{NUC_LONG_AXIS:.6f}"])
        writer.writerow(["nucleus_intermediate_axis_um", f"{NUC_INTERMEDIATE_AXIS:.6f}"])
        writer.writerow(["nucleus_short_axis_um", f"{NUC_SHORT_AXIS:.6f}"])
        writer.writerow(["centrosome_distance_um", f"{CENTROSOME_DISTANCE:.6f}"])

    print("=" * 70)
    print(f"数据已导出到 {output_dir}/ 目录:")
    print(f"  - cell_XXX_progression.csv : 每个细胞的捕捉进度(时间+微管数 vs 捕捉数)")
    print(f"  - aggregate_progression.csv: 各里程碑聚合(时间+微管数)")
    print(f"  - per_cell_summary.csv     : 每个细胞汇总")
    print(f"  - run_summary.csv          : {n_runs} 次独立运行的汇总")
    print(f"  - distal_angles.csv        : distal 偏转角分布(用于作图)")
    print(f"  - parameters_used.csv      : 本次计算使用的参数记录")
    print("=" * 70)

    return {
        "condition": condition_label,
        "condition_family": condition.family,
        "directed_percent": condition.directed_prob * 100.0,
        "pivoting_mode": condition.pivoting_mode,
        "pivoting_parameter": condition_metadata(condition)["pivoting_parameter"],
        "d_mt_rad2_s": condition.d_mt_rad2_s,
        "completed_cells": completed_count,
        "n_runs": n_runs,
        "cells_per_run": n_cells,
        "total_cells": total_cells,
        "mean_completion_time_s": _mean(times),
        "mean_deflection_angle_deg": _mean(all_angles),
        "mean_mt_total_emitted": _mean(mt_counts),
        "output_dir": output_dir,
    }



# -----------------------------------------------------------------------------
# Parallel export worker helpers
# -----------------------------------------------------------------------------
def _run_export_job(job: Dict[str, object]) -> Dict[str, object]:
    """Run one repeat-condition simulation in a worker process.

    This function is intentionally top-level so it can be pickled by
    ProcessPoolExecutor on Windows. It performs CPU-heavy geometry work only;
    Qt widgets and file writing stay in the main process.
    """
    params: SimulationParams = job["params"]  # type: ignore[assignment]
    rep = int(job["repeat"])
    run_index = int(job.get("run", 1))
    cell_index = int(job.get("cell", rep))
    global_cell_index = int(job.get("global_cell", rep))
    layout_seed = int(job["layout_seed"])
    base_units: List[KTUnit] = job["base_units"]  # type: ignore[assignment]
    condition = job.get("condition_spec")
    if not isinstance(condition, ConditionSpec):
        # Compatibility fallback for callers using the original worker shape.
        cond_name = str(job["condition_name"])
        prob = float(job["prob"])
        condition = ConditionSpec(cond_name, prob)
    sim_params = params_for_condition(params, condition)
    include_per_mt = bool(job.get("include_per_mt", False))
    collect_progression = bool(job.get("collect_progression", True))
    collect_deflection = bool(job.get("collect_deflection", True))

    sim_seed = params.seed + rep * 100003 + 17
    sim = EventDrivenSimulator(sim_params, condition.directed_prob, sim_seed, base_units)
    result = sim.run()
    completed = bool(result["completed"])
    n_captured = int(result["n_captured"])
    n_mts = int(result["n_mts"])
    final_time = float(result["final_time_s"])

    common = {
        "run": run_index,
        "cell": cell_index,
        "global_cell": global_cell_index,
        "repeat": rep,
        "layout_seed": layout_seed,
        "simulation_seed": sim_seed,
        **condition_metadata(condition),
        "kt_generation_mode": params.kt_generation_mode,
        "kt_base_diameter_um": params.kt_base_diameter_um,
        "kt_volume_multiplier": params.kt_volume_multiplier,
        "kt_base_volume_um3": params.kt_base_volume_um3,
        "kt_unit_diameter_um": params.kt_diameter_um,
        "kt_effective_volume_um3": params.kt_effective_volume_um3,
        "mt_diameter_um": params.mt_diameter_um,
        "emission_rate_per_centrosome": params.emission_rate_per_centrosome,
        "mt_velocity_um_s": params.mt_velocity_um_s,
        "mt_extra_length_um": params.mt_extra_length_um,
    }

    summary_row = {
        **common,
        "completed": completed,
        "kt_unit_count": TOTAL_KT_UNITS,
        "distal_kt_units": DISTAL_KT_UNITS,
        "proximal_kt_units": PROXIMAL_KT_UNITS,
        "captured_kt_units": n_captured,
        "completion_time_s": final_time if completed else "",
        "mt_total_emitted": n_mts,
        "max_mt_length_um": sim.max_len,
        "angular_occupancy_width_deg": math.degrees(sim.angular_width),
    }
    time_row = {
        "run": run_index,
        "cell": cell_index,
        "global_cell": global_cell_index,
        "repeat": rep,
        "layout_seed": layout_seed,
        "simulation_seed": sim_seed,
        **condition_metadata(condition),
        "completed": completed,
        "completion_time_s": final_time if completed else "",
    }
    per_cell_row = {
        **common,
        "completed": completed,
        "n_captured": n_captured,
        "n_kt": TOTAL_KT_UNITS,
        "n_mts": n_mts,
        "final_time_s": final_time,
    }

    progression_rows = [{**common, **row} for row in result["capture_log"]] if collect_progression else []
    deflection_rows = [{**common, **row} for row in result["deflection_rows"]] if collect_deflection else []
    per_mt_rows: List[Dict[str, object]] = []
    if include_per_mt:
        for mt in result["mt_records"]:
            per_mt_rows.append({
                **common,
                "mt_id": mt.mt_id,
                "cs_index": mt.cs_index,
                "emit_time_s": mt.emit_time,
                "is_directed": mt.is_directed,
                "target_kt_index": "" if mt.target_index is None else mt.target_index,
                "hit_count": len(mt.hits),
                "hit_kt_indices": ";".join(str(h[0]) for h in mt.hits),
                "hit_contact_distances_um": ";".join(f"{h[1]:.6g}" for h in mt.hits),
                "hit_capture_modes": ";".join(mt.hit_capture_modes),
                "no_hit_deflection_to_boundary_deg": "" if mt.no_hit_deflection is None else mt.no_hit_deflection,
                "nearest_kt_index": "" if mt.nearest_kt_index is None else mt.nearest_kt_index,
                "pivot_step_count": mt.pivot_step_count,
                "final_dir_x": "" if mt.final_direction is None else mt.final_direction[0],
                "final_dir_y": "" if mt.final_direction is None else mt.final_direction[1],
                "final_dir_z": "" if mt.final_direction is None else mt.final_direction[2],
                "dir_x": mt.direction[0], "dir_y": mt.direction[1], "dir_z": mt.direction[2],
            })

    return {
        "repeat": rep,
        "condition": condition.name,
        "summary_row": summary_row,
        "time_row": time_row,
        "per_cell_row": per_cell_row,
        "progression_rows": progression_rows,
        "deflection_rows": deflection_rows,
        "per_mt_rows": per_mt_rows,
    }

class ResultsWindow(QWidget):
    """Independent modeless results table showing only the requested metrics.

    Important UI rule: this widget is created as a top-level window, not as a
    child overlay inside the main simulator window. Therefore it will not cover
    or block the main controls; the user can move it aside and continue using
    the main GUI.
    """
    def __init__(self, rows: List[Dict[str, object]], parent: Optional[QWidget] = None):
        # Do NOT pass parent to QWidget, otherwise Qt embeds the result table
        # inside the main window and it overlays/blocks the animation area.
        super().__init__(None)
        self._owner = parent  # keep a weak-style reference for lifetime/debug only
        self.setWindowTitle("Simulation Results")
        self.setWindowFlags(self.windowFlags() | Qt.Window)
        self.setAttribute(Qt.WA_DeleteOnClose, False)
        self.resize(860, 360)
        layout = QVBoxLayout(self)

        title = QLabel("Direct output results")
        title.setStyleSheet("font-size: 16px; font-weight: bold;")
        layout.addWidget(title)

        table = QTableWidget(len(rows), 7)
        table.setHorizontalHeaderLabels([
            "Condition",
            "Scheme",
            "Pivoting parameter",
            "Completed cells",
            "Mean completion time (s)",
            "Mean deflection angle (°)",
            "Mean MT emitted",
        ])

        def fmt(value, digits=2):
            try:
                v = float(value)
                if math.isnan(v):
                    return "NA"
                return f"{v:.{digits}f}"
            except Exception:
                return str(value)

        for r, row in enumerate(rows):
            values = [
                str(row.get("condition", "")),
                str(row.get("condition_family", "")),
                str(row.get("pivoting_parameter", "none")),
                f"{row.get('completed_cells', '')}/{row.get('total_cells', '')}",
                fmt(row.get("mean_completion_time_s", float('nan')), 2),
                fmt(row.get("mean_deflection_angle_deg", float('nan')), 2),
                fmt(row.get("mean_mt_total_emitted", float('nan')), 2),
            ]
            for c, val in enumerate(values):
                table.setItem(r, c, QTableWidgetItem(val))

        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        layout.addWidget(table)

        folders = sorted({str(row.get("output_dir", "")) for row in rows if row.get("output_dir")})
        if folders:
            folder_label = QLabel("CSV source data saved to:\n" + "\n".join(folders))
            folder_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
            layout.addWidget(folder_label)


def run_direct_from_cli(argv: Sequence[str]) -> None:
    parser = argparse.ArgumentParser(description="Run KT-unit MT capture simulation without GUI and export BOSS-style CSV files.")
    parser.add_argument("--direct", action="store_true", help="Run direct output mode instead of opening the GUI.")
    parser.add_argument("--f-directed", type=float, default=DIRECT_F_DIRECTED, help="Directed fraction, e.g. 0, 0.2, 0.4, 0.6, 0.8, 1.0.")
    parser.add_argument("--n-cells", type=int, default=DIRECT_N_CELLS)
    parser.add_argument("--n-runs", type=int, default=DIRECT_N_RUNS, help="Independent runs; each run contains --n-cells cells.")
    parser.add_argument("--output-dir", default=DIRECT_OUTPUT_DIR)
    parser.add_argument("--kt-diameter", type=float, default=0.80, help="Base KT unit diameter in um before volume scaling.")
    parser.add_argument("--kt-volume-multiplier", type=float, default=1.0, help="KT capture volume multiplier. Effective diameter = base diameter * multiplier^(1/3).")
    parser.add_argument("--mt-diameter", type=float, default=0.025, help="MT diameter in um.")
    parser.add_argument("--mt-velocity", type=float, default=0.116, help="MT growth velocity in um/s.")
    parser.add_argument("--emission-rate", type=float, default=0.15, help="Post-burst emission rate per centrosome, MT/s.")
    parser.add_argument("--mt-extra-length", type=float, default=1.0, help="Max MT length = max CS-to-KT distance + this value.")
    parser.add_argument(
        "--pivoting-sweep", choices=("none", "dmt"), default="none",
        help="Run the exact six-condition D_MT pivoting sweep.",
    )
    parser.add_argument(
        "--pivoting-mode", choices=("none", "dmt"), default="none",
        help="Run one condition when --pivoting-sweep=none.",
    )
    parser.add_argument("--d-mt", type=float, default=1e-4, help="D_MT in rad^2/s for one dmt condition.")
    args = parser.parse_args(list(argv))

    kt_base_d = max(0.10, min(2.50, float(args.kt_diameter)))
    kt_vol_mult = max(0.1, min(50.0, float(args.kt_volume_multiplier)))
    params = SimulationParams(
        kt_diameter_um=effective_kt_diameter_um(kt_base_d, kt_vol_mult),
        kt_base_diameter_um=kt_base_d,
        kt_volume_multiplier=kt_vol_mult,
        mt_diameter_um=max(0.025, min(0.500, float(args.mt_diameter))),
        mt_velocity_um_s=max(1e-6, float(args.mt_velocity)),
        emission_rate_per_centrosome=max(0.10, min(2.00, float(args.emission_rate))),
        mt_extra_length_um=max(0.0, float(args.mt_extra_length)),
        seed=new_random_seed(),
    )
    directed_prob = max(0.0, min(1.0, float(args.f_directed)))
    if args.pivoting_sweep == "dmt":
        conditions = [
            ConditionSpec(name, 0.0, "dmt", value, "pivoting_dmt")
            for name, value in PIVOTING_DMT_CONDITIONS
        ]
    elif args.pivoting_mode == "dmt":
        value = max(0.0, float(args.d_mt))
        conditions = [ConditionSpec(f"Pivoting D_MT = {value:.6g} rad^2/s", 0.0, "dmt", value, "pivoting_dmt")]
    else:
        conditions = [ConditionSpec(_condition_label_from_prob(directed_prob), directed_prob)]

    for condition in conditions:
        output_dir = args.output_dir
        if len(conditions) > 1:
            output_dir = os.path.join(
                args.output_dir,
                condition.family,
                _sanitize_filename(condition.name),
            )
        run_boss_style_direct_output(
            params=params,
            directed_prob=condition.directed_prob,
            n_cells=max(1, int(args.n_cells)),
            n_runs=max(1, int(args.n_runs)),
            output_dir=output_dir,
            condition=condition,
        )


# -----------------------------------------------------------------------------
# Main window
# -----------------------------------------------------------------------------
class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("KT-unit MT capture simulator | dual-condition comparison | publication spacing")
        self.current_units: Optional[List[KTUnit]] = None
        self.current_layout_seed: int = new_random_seed()
        self.current_run_seed: int = new_random_seed()
        self.current_animation_seeds: List[int] = [new_random_seed(), new_random_seed()]
        # Backward-compatible alias used by the first-panel SVG export.
        self.current_animation_seed: int = self.current_animation_seeds[0]
        self.running = False
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.tick)

        # Two independently simulated conditions share the same KT layout and GUI parameters.
        # Default comparison matches Movie S10: Random versus 40% directed.
        self.animation_condition_combos: List[QComboBox] = []
        anim_box = QGroupBox("Animation conditions")
        anim_layout = QHBoxLayout()
        default_indices = [0, next((i for i, spec in enumerate(ANIMATION_CONDITIONS) if abs(spec.directed_prob - 0.40) < 1e-12 and spec.pivoting_mode == "none"), 0)]
        for side, default_idx in zip(("Left", "Right"), default_indices):
            anim_layout.addWidget(QLabel(side))
            combo = QComboBox()
            combo.addItems([spec.name for spec in ANIMATION_CONDITIONS])
            combo.setCurrentIndex(default_idx)
            self.animation_condition_combos.append(combo)
            anim_layout.addWidget(combo, 1)
        anim_box.setLayout(anim_layout)
        for combo in self.animation_condition_combos:
            combo.currentIndexChanged.connect(lambda _=0: self.reset_animation(regenerate_layout=False))

        left = QVBoxLayout()
        left.setSpacing(8)
        left.addWidget(anim_box)

        # Speed controls span both complete comparison panels.
        speed_layout = QHBoxLayout()
        speed_layout.addWidget(QLabel("Time acceleration"))
        self.speed_slider = QSlider(Qt.Horizontal)
        self.speed_slider.setRange(1, 200)
        self.speed_slider.setValue(20)
        self.speed_slider.valueChanged.connect(self.update_speed_label)
        speed_layout.addWidget(self.speed_slider, 1)
        self.speed_label = QLabel("20×")
        self.speed_label.setMinimumWidth(56)
        speed_layout.addWidget(self.speed_label)
        left.addLayout(speed_layout)

        self.views: List[SimulationView] = [SimulationView(), SimulationView()]
        self.capture_time_plots: List[LivePlotWidget] = [
            LivePlotWidget("capture_time", "KT captured (%) vs time"),
            LivePlotWidget("capture_time", "KT captured (%) vs time"),
        ]
        self.mt_capture_plots: List[LivePlotWidget] = [
            LivePlotWidget("mt_capture", "MT emitted vs KT captured (%)"),
            LivePlotWidget("mt_capture", "MT emitted vs KT captured (%)"),
        ]
        self.animation_status_labels: List[AnimationStatusWidget] = [AnimationStatusWidget(), AnimationStatusWidget()]
        self.animation_legends: List[AnimationLegendWidget] = [AnimationLegendWidget(), AnimationLegendWidget()]
        self.condition_title_labels: List[QLabel] = [QLabel(), QLabel()]
        self.comparison_modules: List[QFrame] = []

        comparison_row = QHBoxLayout()
        comparison_row.setSpacing(42)
        comparison_row.setContentsMargins(0, 0, 0, 0)

        for i in range(2):
            # Smaller GUI minima keep the two complete panels usable on a desktop screen;
            # MP4 export is independently rendered at ~5200 px width.
            self.views[i].setMinimumSize(420, 320)
            self.capture_time_plots[i].setMinimumSize(210, 150)
            self.mt_capture_plots[i].setMinimumSize(210, 150)
            self.animation_status_labels[i].set_values(0.0, 0, 0)

            plot_col = QVBoxLayout()
            plot_col.setSpacing(8)
            plot_col.setContentsMargins(0, 0, 0, 0)
            plot_col.addWidget(self.capture_time_plots[i], 1)
            plot_col.addWidget(self.mt_capture_plots[i], 1)

            animation_inner_col = QVBoxLayout()
            animation_inner_col.setSpacing(0)
            animation_inner_col.setContentsMargins(0, 0, 0, 0)
            animation_inner_col.addWidget(self.animation_status_labels[i])
            animation_inner_col.addWidget(self.views[i], 1)
            animation_inner_col.addWidget(self.animation_legends[i])

            animation_panel = QFrame()
            animation_panel.setStyleSheet(
                "QFrame {background-color: rgb(3, 9, 12); border: 2px solid rgb(235, 240, 243);}"
            )
            animation_panel_layout = QVBoxLayout(animation_panel)
            animation_panel_layout.setSpacing(0)
            animation_panel_layout.setContentsMargins(7, 7, 7, 7)
            animation_panel_layout.addLayout(animation_inner_col, 1)

            module_visual = QFrame()
            module_visual.setStyleSheet("QFrame {background-color: rgb(0, 0, 0); border: none;}")
            module_visual_layout = QHBoxLayout(module_visual)
            module_visual_layout.setSpacing(8)
            module_visual_layout.setContentsMargins(0, 0, 0, 0)
            module_visual_layout.addWidget(animation_panel, 2)
            module_visual_layout.addLayout(plot_col, 1)

            title = self.condition_title_labels[i]
            title.setAlignment(Qt.AlignCenter)
            title.setFont(QFont(DISPLAY_FONT_FAMILY, visual_font_size(20)))
            title.setStyleSheet("QLabel {color: rgb(248, 248, 248); background-color: rgb(0, 0, 0); border: none;}")
            title.setMinimumHeight(56)

            module = QFrame()
            module.setStyleSheet("QFrame {background-color: rgb(0, 0, 0); border: none;}")
            module_col = QVBoxLayout(module)
            module_col.setSpacing(2)
            module_col.setContentsMargins(0, 0, 0, 0)
            module_col.addWidget(title)
            module_col.addWidget(module_visual, 1)
            self.comparison_modules.append(module)
            comparison_row.addWidget(module, 1)

        self.movie_title_label = QLabel("Movie S10: Three-dimensional simulation of random search-and-capture")
        self.movie_title_label.setAlignment(Qt.AlignCenter)
        self.movie_title_label.setFont(QFont(DISPLAY_FONT_FAMILY, visual_font_size(24)))
        self.movie_title_label.setStyleSheet("QLabel {color: rgb(248, 248, 248); background-color: rgb(0, 0, 0); border: none;}")
        self.movie_title_label.setMinimumHeight(70)

        # This is the exact black publication/video region exported to MP4.
        self.visual_container = QFrame()
        self.visual_container.setStyleSheet("QFrame {background-color: rgb(0, 0, 0); border: none;}")
        visual_outer = QVBoxLayout(self.visual_container)
        visual_outer.setSpacing(0)
        visual_outer.setContentsMargins(0, 0, 0, 0)
        visual_outer.addWidget(self.movie_title_label)
        visual_outer.addSpacing(42)
        visual_outer.addLayout(comparison_row, 1)
        left.addWidget(self.visual_container, 1)

        # First-panel aliases keep legacy frame-export code and other non-comparison
        # utilities working exactly as before.
        self.view = self.views[0]
        self.capture_time_plot = self.capture_time_plots[0]
        self.mt_capture_plot = self.mt_capture_plots[0]
        self.animation_status_label = self.animation_status_labels[0]
        self.animation_legend = self.animation_legends[0]
        self.animation_panel = self.comparison_modules[0]

        right_panel = self.build_right_panel()
        main = QHBoxLayout(self)
        main.addLayout(left, 1)
        main.addWidget(right_panel)

        self.update_condition_titles()
        self.reset_animation(regenerate_layout=True)

    def build_right_panel(self) -> QWidget:
        panel = QFrame()
        panel.setFrameShape(QFrame.StyledPanel)
        layout = QVBoxLayout(panel)
        layout.setSpacing(8)

        params_box = QGroupBox("Parameters")
        form = QFormLayout()
        self.kt_generation_label = QLabel("Uniform random")
        self.kt_generation_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        form.addRow("KT generation mode", self.kt_generation_label)

        self.kt_diam_spin = QDoubleSpinBox()
        self.kt_diam_spin.setRange(0.10, 2.50)
        self.kt_diam_spin.setSingleStep(0.01)
        self.kt_diam_spin.setDecimals(2)
        self.kt_diam_spin.setValue(0.80)
        self.kt_diam_spin.valueChanged.connect(lambda _: self.on_kt_volume_parameter_changed())
        form.addRow("KT base diameter (μm)", self.kt_diam_spin)

        self.kt_volume_mult_spin = QDoubleSpinBox()
        self.kt_volume_mult_spin.setRange(0.1, 50.0)
        self.kt_volume_mult_spin.setSingleStep(0.01)
        self.kt_volume_mult_spin.setDecimals(2)
        self.kt_volume_mult_spin.setValue(1.00)
        self.kt_volume_mult_spin.valueChanged.connect(lambda _: self.on_kt_volume_parameter_changed())
        form.addRow("KT volume multiplier", self.kt_volume_mult_spin)

        self.kt_effective_diameter_label = QLabel("")
        self.kt_effective_diameter_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        form.addRow("Effective diameter (μm)", self.kt_effective_diameter_label)

        self.kt_effective_volume_label = QLabel("")
        self.kt_effective_volume_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        form.addRow("Effective volume (μm³)", self.kt_effective_volume_label)

        self.update_kt_effective_display()

        self.mt_diam_spin = QDoubleSpinBox()
        self.mt_diam_spin.setRange(0.025, 0.500)
        self.mt_diam_spin.setSingleStep(0.005)
        self.mt_diam_spin.setDecimals(3)
        self.mt_diam_spin.setValue(0.025)
        self.mt_diam_spin.valueChanged.connect(lambda _: self.reset_animation(regenerate_layout=False))
        form.addRow("MT diameter (μm)", self.mt_diam_spin)

        self.emission_spin = QDoubleSpinBox()
        self.emission_spin.setRange(0.1000, 2.0000)
        self.emission_spin.setSingleStep(0.0001)
        self.emission_spin.setDecimals(4)
        self.emission_spin.setValue(0.1500)
        self.emission_spin.valueChanged.connect(lambda _: self.reset_animation(regenerate_layout=False))
        form.addRow("Post-burst rate / CS (MT/s)", self.emission_spin)

        self.velocity_spin = QDoubleSpinBox()
        self.velocity_spin.setRange(0.001, 5.000)
        self.velocity_spin.setSingleStep(0.01)
        self.velocity_spin.setDecimals(3)
        self.velocity_spin.setValue(0.116)
        self.velocity_spin.valueChanged.connect(lambda _: self.reset_animation(regenerate_layout=False))
        form.addRow("MT velocity (μm/s)", self.velocity_spin)

        self.extra_len_spin = QDoubleSpinBox()
        self.extra_len_spin.setRange(0.0, 20.0)
        self.extra_len_spin.setSingleStep(0.1)
        self.extra_len_spin.setDecimals(2)
        self.extra_len_spin.setValue(1.0)
        self.extra_len_spin.valueChanged.connect(lambda _: self.reset_animation(regenerate_layout=False))
        form.addRow("MT extra length (μm)", self.extra_len_spin)

        self.runs_spin = QSpinBox()
        self.runs_spin.setRange(1, 100)
        self.runs_spin.setValue(DIRECT_N_RUNS)
        self.runs_spin.setToolTip("Number of independent batch runs for every selected condition.")
        form.addRow("Independent runs", self.runs_spin)

        self.repeats_spin = QSpinBox()
        self.repeats_spin.setRange(1, 10000)
        self.repeats_spin.setValue(DIRECT_N_CELLS)
        self.repeats_spin.setToolTip("Number of cells simulated in each independent run.")
        form.addRow("Cells per run", self.repeats_spin)

        self.workers_spin = QSpinBox()
        self.workers_spin.setRange(1, 12)
        self.workers_spin.setValue(min(12, max(1, os.cpu_count() or 1)))
        self.workers_spin.setToolTip("Number of parallel CPU worker processes for batch export. Use 1 for the old sequential mode.")
        form.addRow("Export workers", self.workers_spin)

        params_box.setLayout(form)
        layout.addWidget(params_box)

        export_cond_box = QGroupBox("Export conditions")
        cond_layout = QVBoxLayout()
        self.export_condition_checks: List[QCheckBox] = []
        for name, _ in EXPORT_CONDITIONS:
            cb = QCheckBox(name)
            cb.setChecked(True)
            self.export_condition_checks.append(cb)
            cond_layout.addWidget(cb)
        export_cond_box.setLayout(cond_layout)
        layout.addWidget(export_cond_box)

        dmt_box = QGroupBox("Pivoting sweep: D_MT (rad²/s)")
        dmt_layout = QVBoxLayout()
        self.pivot_dmt_checks: List[QCheckBox] = []
        for name, _ in PIVOTING_DMT_CONDITIONS:
            cb = QCheckBox(name)
            cb.setChecked(True)
            self.pivot_dmt_checks.append(cb)
            dmt_layout.addWidget(cb)
        dmt_box.setLayout(dmt_layout)
        layout.addWidget(dmt_box)

        outputs_box = QGroupBox("Outputs")
        outputs_layout = QVBoxLayout()
        self.out_time = QCheckBox("Completion time")
        self.out_mt_total = QCheckBox("MT total emitted")
        self.out_summary = QCheckBox("Summary")
        self.out_per_cell = QCheckBox("Per-cell summary")
        self.out_progress = QCheckBox("Capture progression")
        self.out_deflection = QCheckBox("Deflection angles")
        self.out_parameters = QCheckBox("Parameters used")
        self.out_per_mt = QCheckBox("Per-MT details")
        self.out_kt_layout = QCheckBox("KT layout")
        for cb in [self.out_summary, self.out_per_cell, self.out_progress, self.out_deflection, self.out_parameters]:
            cb.setChecked(True)
        for cb in [self.out_time, self.out_mt_total, self.out_summary, self.out_per_cell, self.out_progress, self.out_deflection, self.out_parameters, self.out_per_mt, self.out_kt_layout]:
            outputs_layout.addWidget(cb)
        outputs_box.setLayout(outputs_layout)
        layout.addWidget(outputs_box)

        btn_layout = QVBoxLayout()
        self.regen_btn = QPushButton("Regenerate KT layout")
        self.regen_btn.clicked.connect(lambda: self.reset_animation(regenerate_layout=True, force_new_seed=True))
        btn_layout.addWidget(self.regen_btn)
        self.start_btn = QPushButton("Start")
        self.start_btn.clicked.connect(self.toggle_start)
        btn_layout.addWidget(self.start_btn)
        self.reset_btn = QPushButton("Reset animation")
        self.reset_btn.clicked.connect(lambda: self.reset_animation(regenerate_layout=False))
        btn_layout.addWidget(self.reset_btn)

        self.export_frame_btn = QPushButton("Export current frame (SVG)")
        self.export_frame_btn.clicked.connect(self.export_current_animation_frame_svg)
        btn_layout.addWidget(self.export_frame_btn)

        self.export_movie_btn = QPushButton("Export full capture video (MP4)")
        self.export_movie_btn.setToolTip(
            "Export the complete two-condition comparison at 30 FPS using the current animation speed. "
            "The Movie S10 title, both condition titles, both full simulation panels, and plots are included."
        )
        self.export_movie_btn.clicked.connect(self.export_full_capture_video)
        btn_layout.addWidget(self.export_movie_btn)

        self.export_btn = QPushButton("Export selected results")
        self.export_btn.clicked.connect(self.export_results)
        btn_layout.addWidget(self.export_btn)

        self.direct_output_btn = QPushButton("Direct output export conditions")
        self.direct_output_btn.clicked.connect(self.direct_output_current_condition)
        btn_layout.addWidget(self.direct_output_btn)

        layout.addLayout(btn_layout)
        layout.addStretch(1)
        panel.setFixedWidth(360)
        return panel

    def update_kt_effective_display(self) -> None:
        base_d = float(self.kt_diam_spin.value())
        mult = float(self.kt_volume_mult_spin.value())
        eff_d = effective_kt_diameter_um(base_d, mult)
        base_v = kt_sphere_volume_um3(base_d)
        eff_v = kt_sphere_volume_um3(eff_d)
        self.kt_effective_diameter_label.setText(f"{eff_d:.3f}")
        self.kt_effective_diameter_label.setToolTip(
            f"Effective diameter = {eff_d:.6f} μm; base diameter = {base_d:.6f} μm; multiplier = {mult:.2f}×"
        )
        self.kt_effective_volume_label.setText(f"{eff_v:.3f}  (base {base_v:.3f}; {mult:.2f}×)")
        self.kt_effective_volume_label.setToolTip(
            f"Effective KT volume = {eff_v:.6f} μm³; base KT volume = {base_v:.6f} μm³; multiplier = {mult:.2f}×"
        )

    def on_kt_volume_parameter_changed(self) -> None:
        self.update_kt_effective_display()
        self.reset_animation(regenerate_layout=True)

    def current_params(self) -> SimulationParams:
        base_d = float(self.kt_diam_spin.value())
        mult = float(self.kt_volume_mult_spin.value())
        eff_d = effective_kt_diameter_um(base_d, mult)
        return SimulationParams(
            kt_generation_mode="Uniform random",
            kt_diameter_um=eff_d,
            kt_base_diameter_um=base_d,
            kt_volume_multiplier=mult,
            mt_diameter_um=float(self.mt_diam_spin.value()),
            mt_velocity_um_s=float(self.velocity_spin.value()),
            emission_rate_per_centrosome=float(self.emission_spin.value()),
            mt_extra_length_um=float(self.extra_len_spin.value()),
            seed=int(self.current_run_seed),
        )

    def _display_condition_name(self, name: str) -> str:
        return str(name).replace(" directed", " Directed")

    def selected_animation_conditions(self) -> List[ConditionSpec]:
        selected: List[ConditionSpec] = []
        for combo in self.animation_condition_combos:
            idx = max(0, min(len(ANIMATION_CONDITIONS) - 1, int(combo.currentIndex())))
            selected.append(ANIMATION_CONDITIONS[idx])
        return selected

    def selected_animation_condition(self) -> ConditionSpec:
        # Backward-compatible first-panel selection used by the SVG export.
        return self.selected_animation_conditions()[0]

    def update_condition_titles(self) -> None:
        for label, condition in zip(self.condition_title_labels, self.selected_animation_conditions()):
            label.setText(self._display_condition_name(condition.name))

    def selected_export_conditions(self) -> List[ConditionSpec]:
        selected = [
            ConditionSpec(name, prob)
            for (name, prob), cb in zip(EXPORT_CONDITIONS, self.export_condition_checks)
            if cb.isChecked()
        ]
        selected.extend(
            ConditionSpec(name, 0.0, "dmt", value, "pivoting_dmt")
            for (name, value), cb in zip(PIVOTING_DMT_CONDITIONS, self.pivot_dmt_checks)
            if cb.isChecked()
        )
        return selected

    def reset_animation(self, regenerate_layout: bool, force_new_seed: bool = False) -> None:
        # Both comparison conditions share one KT layout so only the directed
        # fraction differs spatially; each animation receives its own fresh RNG seed.
        self.current_run_seed = new_random_seed()
        params = self.current_params()
        if regenerate_layout or force_new_seed or self.current_units is None:
            self.current_layout_seed = new_random_seed()
            try:
                self.current_units = generate_kt_layout(params, self.current_layout_seed)
            except Exception as e:
                QMessageBox.warning(self, "KT layout generation failed", str(e))
                return

        self.current_animation_seeds = [new_random_seed(), new_random_seed()]
        self.current_animation_seed = self.current_animation_seeds[0]
        conditions = self.selected_animation_conditions()
        self.update_condition_titles()
        for i, condition in enumerate(conditions):
            self.views[i].set_simulation(
                AnimationSimulation(
                    params_for_condition(params, condition),
                    condition.directed_prob,
                    self.current_units,
                    self.current_animation_seeds[i],
                )
            )
            self.capture_time_plots[i].clear()
            self.mt_capture_plots[i].clear()

        self.update_animation_status()
        self.running = False
        self.timer.stop()
        self.start_btn.setText("Start")

    def update_speed_label(self) -> None:
        value = max(1, min(200, int(self.speed_slider.value())))
        self.speed_label.setText(f"{value}×")
        for view in self.views:
            view.set_time_acceleration(value)

    def update_animation_status(self) -> None:
        # Keep both on-screen/exported Time acceleration annotations synchronized.
        acceleration = max(1, min(200, int(self.speed_slider.value())))
        for i, view in enumerate(self.views):
            view.set_time_acceleration(acceleration)
            sim = view.sim
            if sim is None:
                self.animation_status_labels[i].set_values(0.0, 0, 0)
                continue
            captured = sum(1 for kt in sim.kt_units if kt.displayed_captured)
            self.animation_status_labels[i].set_values(sim.time, sim.mt_id_counter, captured)

    def toggle_start(self) -> None:
        self.running = not self.running
        if self.running:
            self.timer.start(16)
            self.start_btn.setText("Pause")
        else:
            self.timer.stop()
            self.start_btn.setText("Start")

    def tick(self) -> None:
        # Both conditions advance by the same simulated dt, so time is directly comparable.
        dt = 0.016 * max(1, int(self.speed_slider.value()))
        for i, view in enumerate(self.views):
            if view.sim is None:
                continue
            view.sim.advance(dt)
            self.capture_time_plots[i].update_from_sim(view.sim)
            self.mt_capture_plots[i].update_from_sim(view.sim)
            view.update()
        self.update_animation_status()

    def _condition_output_folder_name(self, condition_name: str) -> str:
        """Folder-safe condition name for BOSS-style direct outputs."""
        return condition_name.replace("%", "pct").replace(" ", "_").replace("/", "_")

    def direct_output_current_condition(self) -> None:
        """Run BOSS-style direct output for the selected Export conditions.

        Animation condition is only for visualization and is intentionally not
        used here. This button uses the current GUI parameters plus the checked
        Export conditions, then writes BOSS-style files directly to
        ./capture_sim_output.

        If one export condition is selected, files are written directly under
        ./capture_sim_output, matching the original BOSS script. If multiple
        conditions are selected, each condition is written into its own
        subfolder to avoid overwriting files.
        """
        params = self.current_params()
        conditions = self.selected_export_conditions()
        if not conditions:
            QMessageBox.information(self, "No export condition", "Please select at least one Export condition.")
            return

        cells_per_run = int(self.repeats_spin.value())
        n_runs = int(self.runs_spin.value())
        base_output_dir = os.path.join(os.getcwd(), DIRECT_OUTPUT_DIR)

        was_running = self.running
        if self.running:
            self.running = False
            self.timer.stop()
            self.start_btn.setText("Start")

        try:
            self.direct_output_btn.setEnabled(False)
            QApplication.processEvents()

            result_rows: List[Dict[str, object]] = []
            for condition in conditions:
                if len(conditions) == 1:
                    output_dir = base_output_dir
                else:
                    output_dir = os.path.join(base_output_dir, self._condition_output_folder_name(condition.name))
                result_rows.append(run_boss_style_direct_output(
                    params=params,
                    directed_prob=condition.directed_prob,
                    n_cells=cells_per_run,
                    n_runs=n_runs,
                    output_dir=output_dir,
                    progress_milestones=PROGRESS_MILESTONES,
                    condition=condition,
                ))
        except Exception as e:
            QMessageBox.warning(self, "Direct output failed", str(e))
            return
        finally:
            self.direct_output_btn.setEnabled(True)
            if was_running:
                # Do not automatically resume; direct output may have taken time
                # and the user may want to inspect the result first.
                pass

        # Independent top-level table window. It is intentionally not embedded
        # in the main simulator window, so it does not block later operations.
        self.results_window = ResultsWindow(result_rows, None)
        self.results_window.show()
        self.results_window.raise_()
        self.results_window.activateWindow()

    def _find_ffmpeg_executable(self) -> Optional[str]:
        """Return a usable FFmpeg executable if available.

        Preference order:
        1. ffmpeg / ffmpeg.exe on PATH
        2. ffmpeg executable placed next to this script
        3. imageio-ffmpeg's bundled executable, if imageio-ffmpeg is installed
        """
        exe = shutil.which("ffmpeg") or shutil.which("ffmpeg.exe")
        if exe:
            return exe

        local_candidates = [
            os.path.join(script_dir(), "ffmpeg.exe"),
            os.path.join(script_dir(), "ffmpeg"),
        ]
        for candidate in local_candidates:
            if os.path.isfile(candidate):
                return candidate

        try:
            import imageio_ffmpeg  # type: ignore
            candidate = imageio_ffmpeg.get_ffmpeg_exe()
            if candidate and os.path.isfile(candidate):
                return candidate
        except Exception:
            pass
        return None

    def _render_widget_rgba_bytes(self, widget: QWidget, physical_w: int, physical_h: int) -> bytes:
        """Render a QWidget to high-resolution RGBA bytes for FFmpeg.

        The painter is scaled from the current logical widget geometry to the
        requested output pixel dimensions. Text, lines, and plot strokes are
        therefore re-rasterized at export resolution instead of enlarging an
        already-rendered low-resolution frame.
        """
        logical_w = max(1, int(widget.width()))
        logical_h = max(1, int(widget.height()))
        image = QImage(physical_w, physical_h, QImage.Format_RGBA8888)
        image.fill(QColor(0, 0, 0))
        painter = QPainter(image)
        try:
            painter.setRenderHint(QPainter.Antialiasing, True)
            painter.scale(physical_w / float(logical_w), physical_h / float(logical_h))
            widget.render(painter, QPoint(0, 0))
        finally:
            painter.end()
        return bytes(image.bits())

    def export_full_capture_video(self) -> None:
        """Offline-render the Movie S10 two-condition comparison to high-resolution MP4."""
        if self.current_units is None or any(view.sim is None for view in self.views):
            QMessageBox.information(self, "No animation", "Please generate the comparison animation before exporting a video.")
            return

        ffmpeg_exe = self._find_ffmpeg_executable()
        if not ffmpeg_exe:
            QMessageBox.warning(
                self,
                "FFmpeg not found",
                "MP4 export needs FFmpeg.\n\n"
                "Install it with:\n"
                "    pip install imageio-ffmpeg\n\n"
                "or place ffmpeg.exe next to this Python script, then run the export again.",
            )
            return

        conditions = self.selected_animation_conditions()
        safe_pair = "_vs_".join(_sanitize_filename(condition.name) for condition in conditions)
        default_name = f"capture_movie_{safe_pair}.mp4"
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export dual-condition full capture video",
            os.path.join(os.getcwd(), default_name),
            "MP4 video (*.mp4)",
        )
        if not path:
            return
        if not path.lower().endswith(".mp4"):
            path += ".mp4"

        fps = 30
        initial_acceleration = max(1, min(200, int(self.speed_slider.value())))
        min_acceleration_used = initial_acceleration
        max_acceleration_used = initial_acceleration

        # Pause only live playback. The visible comparison remains untouched and
        # the speed slider stays accessible while the hidden high-res compositor records.
        was_running = self.running
        self.running = False
        self.timer.stop()
        self.start_btn.setText("Start")

        params = self.current_params()
        export_sims = [
            AnimationSimulation(
                params_for_condition(params, condition),
                condition.directed_prob,
                self.current_units,
                self.current_animation_seeds[i],
            )
            for i, condition in enumerate(conditions)
        ]
        export_widget = DualComparisonExportWidget(conditions, export_sims, self.views)
        export_widget.update_frame(initial_acceleration)
        export_widget.layout().activate()
        QApplication.processEvents()

        # Fixed logical compositor = 3168 x 1120. Re-rasterize it at ~5600 px
        # width, giving an output near 5600 x 1980 px with no AE downsampling.
        logical_w = int(export_widget.width())
        logical_h = int(export_widget.height())
        render_scale = 5600.0 / float(logical_w)
        physical_w = max(2, int(round(logical_w * render_scale)))
        physical_h = max(2, int(round(logical_h * render_scale)))
        out_w = physical_w if physical_w % 2 == 0 else physical_w + 1
        out_h = physical_h if physical_h % 2 == 0 else physical_h + 1

        cmd = [
            ffmpeg_exe, "-y", "-loglevel", "error",
            "-f", "rawvideo", "-pix_fmt", "rgba",
            "-s", f"{physical_w}x{physical_h}", "-r", str(fps),
            "-i", "-", "-an",
            "-c:v", "libx264", "-preset", "medium", "-crf", "16",
            "-vf", f"pad={out_w}:{out_h}:0:0:black",
            "-pix_fmt", "yuv420p", "-movflags", "+faststart", path,
        ]

        creationflags = getattr(subprocess, "CREATE_NO_WINDOW", 0) if os.name == "nt" else 0
        progress = QProgressDialog(
            "Recording Movie S10 comparison...\nTime acceleration remains adjustable from 1× to 200×.",
            "Cancel", 0, 0, self
        )
        progress.setWindowTitle("Export full capture comparison")
        progress.setWindowModality(Qt.NonModal)
        progress.setMinimumDuration(0)
        progress.show()

        proc = None
        frame_count = 0
        cancelled = False
        completed = False
        error_text = ""
        max_frames = fps * 20 * 60

        try:
            proc = subprocess.Popen(
                cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                creationflags=creationflags,
            )
            if proc.stdin is None:
                raise RuntimeError("Could not open the FFmpeg input pipe.")

            next_frame_wall_time = time.perf_counter()
            while frame_count < max_frames:
                QApplication.processEvents()
                if progress.wasCanceled():
                    cancelled = True
                    break

                current_acceleration = max(1, min(200, int(self.speed_slider.value())))
                min_acceleration_used = min(min_acceleration_used, current_acceleration)
                max_acceleration_used = max(max_acceleration_used, current_acceleration)
                export_widget.update_frame(current_acceleration)
                export_widget.layout().activate()
                QApplication.processEvents()

                proc.stdin.write(
                    self._render_widget_rgba_bytes(export_widget, physical_w, physical_h)
                )
                frame_count += 1

                completed_flags = []
                for sim in export_sims:
                    all_displayed = all(k.displayed_captured for k in sim.kt_units)
                    all_retracted = all((not k.retracting) for k in sim.kt_units if k.displayed_captured)
                    completed_flags.append(all_displayed and all_retracted)
                if all(completed_flags):
                    completed = True
                    break

                for sim in export_sims:
                    sim.advance(current_acceleration / float(fps))

                if frame_count % 10 == 0:
                    cap_text = "   |   ".join(
                        f"{self._display_condition_name(condition.name)}: "
                        f"{sum(1 for k in sim.kt_units if k.displayed_captured)}/{TOTAL_KT_UNITS}"
                        for condition, sim in zip(conditions, export_sims)
                    )
                    progress.setLabelText(
                        f"{cap_text}\n"
                        f"Simulation time {max(sim.time for sim in export_sims):.1f} s   "
                        f"Video {frame_count / fps:.1f} s   "
                        f"Time acceleration {current_acceleration}×"
                    )

                next_frame_wall_time += 1.0 / float(fps)
                while True:
                    remaining = next_frame_wall_time - time.perf_counter()
                    if remaining <= 0:
                        break
                    QApplication.processEvents()
                    if progress.wasCanceled():
                        cancelled = True
                        break
                    time.sleep(min(0.005, remaining))
                if cancelled:
                    break

            if completed and not cancelled:
                final_raw = self._render_widget_rgba_bytes(export_widget, physical_w, physical_h)
                for _ in range(fps):
                    proc.stdin.write(final_raw)
                    frame_count += 1

            try:
                proc.stdin.close()
            except Exception:
                pass

            if cancelled:
                proc.terminate()
                proc.wait(timeout=5)
            else:
                return_code = proc.wait()
                if proc.stderr is not None:
                    try:
                        error_text = proc.stderr.read().decode("utf-8", errors="replace").strip()
                    except Exception:
                        error_text = ""
                if return_code != 0:
                    raise RuntimeError(error_text or f"FFmpeg exited with code {return_code}.")

        except Exception as e:
            error_text = str(e)
            if proc is not None:
                try:
                    proc.terminate()
                except Exception:
                    pass
        finally:
            progress.close()
            export_widget.deleteLater()
            if was_running:
                self.running = True
                self.timer.start(16)
                self.start_btn.setText("Pause")

        if cancelled:
            try:
                if os.path.exists(path):
                    os.remove(path)
            except Exception:
                pass
            QMessageBox.information(self, "Video export cancelled", "The video export was cancelled.")
            return

        if error_text:
            try:
                if os.path.exists(path) and os.path.getsize(path) == 0:
                    os.remove(path)
            except Exception:
                pass
            QMessageBox.warning(self, "Video export failed", error_text)
            return

        if not completed:
            QMessageBox.warning(
                self, "Video export stopped",
                "The safety limit was reached before both conditions captured all kinetochores. "
                "Increase Time acceleration and export again.",
            )
            return

        if not os.path.exists(path) or os.path.getsize(path) <= 0:
            QMessageBox.warning(self, "Video export failed", "The MP4 file was not created correctly.")
            return

        QMessageBox.information(
            self, "Video export complete",
            f"Exported dual-condition Movie S10 comparison to:\n{path}\n\n"
            f"Conditions: {self._display_condition_name(conditions[0].name)} vs "
            f"{self._display_condition_name(conditions[1].name)}\n"
            f"FPS: {fps}\n"
            f"Time acceleration used: {min_acceleration_used}×–{max_acceleration_used}×\n"
            f"Direct high-resolution render: {physical_w} x {physical_h} px\n"
            f"Encoded resolution: {out_w} x {out_h} px\n"
            f"Frames: {frame_count}\n"
            f"Duration: {frame_count / fps:.1f} s",
        )

    def export_current_animation_frame_svg(self) -> None:
        """Export the current 3D animation frame as a scalable vector SVG.

        The exported frame contains the black status strip plus the 3D
        animation view only. Live plots and control panels are intentionally
        excluded so the file can be used as a clean figure panel.
        """
        sim = self.view.sim
        if sim is None:
            QMessageBox.information(self, "No animation frame", "Please generate an animation before exporting a frame.")
            return

        # Freeze the current visual state before exporting, so the saved SVG
        # corresponds exactly to the frame the user is inspecting.
        if self.running:
            self.running = False
            self.timer.stop()
            self.start_btn.setText("Start")

        condition = self.selected_animation_condition()
        safe_condition = _sanitize_filename(condition.name)
        default_name = f"animation_frame_{safe_condition}_{sim.time:.1f}s.svg"
        default_path = os.path.join(os.getcwd(), default_name)
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export current animation frame as SVG",
            default_path,
            "SVG vector image (*.svg)",
        )
        if not path:
            return
        if not path.lower().endswith(".svg"):
            path += ".svg"

        self.update_animation_status()
        self.animation_status_label.update()
        self.view.update()
        QApplication.processEvents()

        label_w = max(1, self.animation_status_label.width())
        label_h = max(1, self.animation_status_label.height())
        view_w = max(1, self.view.width())
        view_h = max(1, self.view.height())
        source_w = max(label_w, view_w)
        source_h = label_h + view_h

        # SVG remains vector-scalable, but a large nominal viewBox makes the
        # file open crisply in Illustrator/PowerPoint and keeps exported PNGs
        # high-resolution if the SVG is rasterized later.
        export_w = 3000
        scale = export_w / float(source_w)
        export_h = max(1, int(math.ceil(source_h * scale)))

        generator = QSvgGenerator()
        generator.setFileName(path)
        generator.setSize(QSize(export_w, export_h))
        generator.setViewBox(QRect(0, 0, export_w, export_h))
        generator.setTitle("Current KT-unit MT capture animation frame")
        generator.setDescription("Black status strip plus current 3D animation view exported from the simulator.")

        painter = QPainter(generator)
        try:
            painter.setRenderHint(QPainter.Antialiasing, True)
            painter.scale(scale, scale)
            painter.fillRect(0, 0, source_w, source_h, COL_BG)

            # Render the status bar and the 3D view as one clean panel.
            self.animation_status_label.render(painter, QPoint(0, 0))
            painter.translate(0, label_h)
            self.view.render(painter, QPoint(0, 0))
        finally:
            painter.end()

        file_ok = os.path.exists(path) and os.path.getsize(path) > 0
        if not file_ok:
            QMessageBox.warning(self, "Frame export failed", "The SVG file was not created correctly.")
            return

        QMessageBox.information(
            self,
            "Frame export complete",
            f"Exported current animation frame to:\n{path}\n\n"
            f"QC summary:\n"
            f"- Format: SVG vector image\n"
            f"- Exported area: status strip + 3D animation view\n"
            f"- Nominal canvas: {export_w} x {export_h} px\n"
            f"- File size: {os.path.getsize(path)} bytes"
        )

    def export_results(self) -> None:
        conditions = self.selected_export_conditions()
        if not conditions:
            QMessageBox.information(self, "No export condition", "Please select at least one export condition.")
            return
        out_dir = QFileDialog.getExistingDirectory(self, "Choose output folder", os.getcwd())
        if not out_dir:
            return

        params = self.current_params()
        # Fresh base seed for this export, so separate export runs no longer reproduce identical tables.
        params.seed = new_random_seed()
        n_runs = int(self.runs_spin.value())
        cells_per_run = int(self.repeats_spin.value())
        workers = min(12, max(1, int(self.workers_spin.value())))
        total_jobs = n_runs * cells_per_run * len(conditions)
        condition_order = [condition.name for condition in conditions]
        condition_rank = {name: i for i, name in enumerate(condition_order)}

        selected_output_checks = [
            self.out_summary.isChecked(),
            self.out_per_cell.isChecked(),
            self.out_progress.isChecked(),
            self.out_deflection.isChecked(),
            self.out_parameters.isChecked(),
            self.out_per_mt.isChecked(),
            self.out_kt_layout.isChecked(),
        ]
        metric_only = (self.out_time.isChecked() or self.out_mt_total.isChecked()) and not any(selected_output_checks)
        include_per_mt = self.out_per_mt.isChecked()
        # Full superplot export still requires capture progression and deflection rows.
        collect_progression = (not metric_only) or self.out_progress.isChecked()
        collect_deflection = (not metric_only) or self.out_deflection.isChecked()

        progress = QProgressDialog(
            f"Running event-driven export with {workers} worker(s)...",
            "Cancel", 0, total_jobs, self
        )
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)

        summary_rows: List[Dict[str, object]] = []
        run_summary_rows: List[Dict[str, object]] = []
        time_rows: List[Dict[str, object]] = []
        mt_total_rows: List[Dict[str, object]] = []
        per_cell_rows: List[Dict[str, object]] = []
        progression_rows: List[Dict[str, object]] = []
        deflection_rows: List[Dict[str, object]] = []
        per_mt_rows: List[Dict[str, object]] = []
        kt_layout_rows: List[Dict[str, object]] = []
        parameters_rows: List[Dict[str, object]] = []

        # Generate each repeat's KT layout once in the main process, then reuse
        # that same layout across all selected export conditions. This preserves
        # the original comparison design while allowing the heavy simulations to
        # run in parallel.
        jobs: List[Dict[str, object]] = []
        try:
            for run_index in range(1, n_runs + 1):
                for cell_index in range(1, cells_per_run + 1):
                    global_cell_index = (run_index - 1) * cells_per_run + cell_index
                    layout_seed = params.seed + global_cell_index * 100003
                    base_units = generate_kt_layout(params, layout_seed)

                    if self.out_kt_layout.isChecked():
                        for kt in base_units:
                            kt_layout_rows.append({
                                "run": run_index,
                                "cell": cell_index,
                                "global_cell": global_cell_index,
                                "repeat": global_cell_index,
                                "layout_seed": layout_seed,
                                "kt_index": kt.index,
                                "x": kt.pos[0], "y": kt.pos[1], "z": kt.pos[2],
                                "is_distal": kt.is_distal,
                                "kt_group": kt.group,
                                "kt_generation_mode": params.kt_generation_mode,
                                "kt_base_diameter_um": params.kt_base_diameter_um,
                                "kt_volume_multiplier": params.kt_volume_multiplier,
                                "kt_base_volume_um3": params.kt_base_volume_um3,
                                "kt_unit_diameter_um": params.kt_diameter_um,
                                "kt_effective_volume_um3": params.kt_effective_volume_um3,
                            })

                    for condition in conditions:
                        jobs.append({
                            "params": params,
                            "repeat": global_cell_index,
                            "run": run_index,
                            "cell": cell_index,
                            "global_cell": global_cell_index,
                            "layout_seed": layout_seed,
                            "base_units": base_units,
                            "condition_spec": condition,
                            "include_per_mt": include_per_mt,
                            "collect_progression": collect_progression,
                            "collect_deflection": collect_deflection,
                        })
        except Exception as e:
            QMessageBox.warning(self, "KT layout generation failed", str(e))
            return

        start_time = time.perf_counter()
        completed_jobs = 0

        def merge_worker_result(result: Dict[str, object]) -> None:
            summary_rows.append(result["summary_row"])  # type: ignore[arg-type]
            time_rows.append(result["time_row"])  # type: ignore[arg-type]
            per_cell_rows.append(result["per_cell_row"])  # type: ignore[arg-type]
            progression_rows.extend(result["progression_rows"])  # type: ignore[arg-type]
            deflection_rows.extend(result["deflection_rows"])  # type: ignore[arg-type]
            per_mt_rows.extend(result["per_mt_rows"])  # type: ignore[arg-type]

        try:
            if workers <= 1:
                # Old deterministic sequential mode, kept for debugging/QC.
                for job in jobs:
                    if progress.wasCanceled():
                        return
                    merge_worker_result(_run_export_job(job))
                    completed_jobs += 1
                    progress.setValue(completed_jobs)
                    QApplication.processEvents()
            else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
                    future_to_key = {executor.submit(_run_export_job, job): (job["repeat"], job["condition_spec"].name) for job in jobs}
                    for future in concurrent.futures.as_completed(future_to_key):
                        if progress.wasCanceled():
                            for f in future_to_key:
                                f.cancel()
                            executor.shutdown(wait=False, cancel_futures=True)
                            return
                        try:
                            merge_worker_result(future.result())
                        except Exception as e:
                            rep, cond_name = future_to_key[future]
                            QMessageBox.warning(self, "Export worker failed", f"Repeat {rep}, {cond_name}: {e}")
                            return
                        completed_jobs += 1
                        progress.setValue(completed_jobs)
                        QApplication.processEvents()
        except Exception as e:
            QMessageBox.warning(self, "Parallel export failed", str(e))
            return

        # Sort rows after parallel completion so output tables remain stable and
        # easy to compare with old sequential exports.
        def cond_idx(row: Dict[str, object]) -> int:
            return condition_rank.get(str(row.get("condition")), 999)
        def safe_int(value: object) -> int:
            try:
                return int(value)
            except Exception:
                return 0
        summary_rows.sort(key=lambda r: (safe_int(r.get("repeat", 0)), cond_idx(r)))
        time_rows.sort(key=lambda r: (safe_int(r.get("repeat", 0)), cond_idx(r)))
        mt_total_rows = [{
            "run": r.get("run", ""),
            "cell": r.get("cell", ""),
            "global_cell": r.get("global_cell", ""),
            "repeat": r.get("repeat", ""),
            "layout_seed": r.get("layout_seed", ""),
            "simulation_seed": r.get("simulation_seed", ""),
            "condition": r.get("condition", ""),
            "condition_family": r.get("condition_family", ""),
            "directed_percent": r.get("directed_percent", ""),
            "pivoting_mode": r.get("pivoting_mode", ""),
            "pivoting_parameter": r.get("pivoting_parameter", ""),
            "d_mt_rad2_s": r.get("d_mt_rad2_s", ""),
            "completed": r.get("completed", ""),
            "mt_total_emitted": r.get("mt_total_emitted", ""),
        } for r in summary_rows]
        mt_total_rows.sort(key=lambda r: (safe_int(r.get("repeat", 0)), cond_idx(r)))
        per_cell_rows.sort(key=lambda r: (safe_int(r.get("repeat", 0)), cond_idx(r)))
        progression_rows.sort(key=lambda r: (
            safe_int(r.get("repeat", 0)), cond_idx(r),
            _safe_float(r.get("capture_time_s")) if _safe_float(r.get("capture_time_s")) is not None else float("inf"),
        ))
        deflection_rows.sort(key=lambda r: (safe_int(r.get("repeat", 0)), cond_idx(r), safe_int(r.get("mt_id", 0))))
        per_mt_rows.sort(key=lambda r: (safe_int(r.get("repeat", 0)), cond_idx(r), safe_int(r.get("mt_id", 0))))

        for run_index in range(1, n_runs + 1):
            for condition in conditions:
                condition_rows = [
                    row for row in summary_rows
                    if int(row.get("run", 0)) == run_index
                    and str(row.get("condition", "")) == condition.name
                ]
                completed_condition_rows = [
                    row for row in condition_rows if bool(row.get("completed", False))
                ]
                run_summary_rows.append({
                    "run": run_index,
                    "cells_requested": cells_per_run,
                    "cells_simulated": len(condition_rows),
                    "completed_cells": len(completed_condition_rows),
                    "mean_completion_time_s": _mean([
                        float(row["completion_time_s"])
                        for row in completed_condition_rows
                        if _safe_float(row.get("completion_time_s")) is not None
                    ]),
                    "mean_mt_total_emitted": _mean([
                        float(row["mt_total_emitted"])
                        for row in completed_condition_rows
                        if _safe_float(row.get("mt_total_emitted")) is not None
                    ]),
                    **condition_metadata(condition),
                })

        if self.out_time.isChecked():
            write_csv(os.path.join(out_dir, "completion_time_s.csv"), time_rows)
        if self.out_mt_total.isChecked():
            write_csv(os.path.join(out_dir, "mt_total_emitted.csv"), mt_total_rows)
        if self.out_summary.isChecked():
            write_csv(os.path.join(out_dir, "summary.csv"), summary_rows)
            write_csv(os.path.join(out_dir, "run_summary.csv"), run_summary_rows)
        if self.out_per_cell.isChecked():
            write_csv(os.path.join(out_dir, "per_cell_summary.csv"), per_cell_rows)
        if self.out_progress.isChecked():
            write_csv(os.path.join(out_dir, "capture_progression.csv"), progression_rows)
        if self.out_deflection.isChecked():
            write_csv(os.path.join(out_dir, "deflection_angles.csv"), deflection_rows)
        if self.out_parameters.isChecked():
            parameters_rows.append({
                "kt_generation_mode": params.kt_generation_mode,
                "kt_base_diameter_um": params.kt_base_diameter_um,
                "kt_volume_multiplier": params.kt_volume_multiplier,
                "kt_base_volume_um3": params.kt_base_volume_um3,
                "kt_unit_diameter_um": params.kt_diameter_um,
                "kt_effective_volume_um3": params.kt_effective_volume_um3,
                "mt_diameter_um": params.mt_diameter_um,
                "post_burst_emission_rate_per_centrosome_MT_per_s": params.emission_rate_per_centrosome,
                "mt_velocity_um_s": params.mt_velocity_um_s,
                "mt_extra_length_um": params.mt_extra_length_um,
                "pivot_dt_s": PIVOT_DT_S,
                "theta_sample_dt_s": THETA_SAMPLE_DT_S,
                "pivoting_schemes": "dmt",
                "d_mt_conditions_rad2_s": ";".join(f"{value:.9g}" for _, value in PIVOTING_DMT_CONDITIONS),
                "number_of_runs": n_runs,
                "cells_per_run": cells_per_run,
                "total_cells_per_condition": n_runs * cells_per_run,
                "seed_mode": "automatic_random_each_export",
                "auto_random_seed_for_this_export": params.seed,
                "export_conditions": ";".join(condition.name for condition in conditions),
                "output_completion_time": self.out_time.isChecked(),
                "output_mt_total_emitted": self.out_mt_total.isChecked(),
                "export_workers": workers,
                "total_export_jobs": total_jobs,
                "parallel_backend": "ProcessPoolExecutor" if workers > 1 else "sequential",
                "total_kt_units": TOTAL_KT_UNITS,
                "distal_kt_units": DISTAL_KT_UNITS,
                "proximal_kt_units": PROXIMAL_KT_UNITS,
                "first_burst_per_centrosome": FIRST_BURST_PER_CENTROSOME,
                "first_burst_total": FIRST_BURST_PER_CENTROSOME * 2,
            })
            write_csv(os.path.join(out_dir, "parameters_used.csv"), parameters_rows)
        if self.out_per_mt.isChecked():
            write_csv(os.path.join(out_dir, "per_mt_details.csv"), per_mt_rows)
        if self.out_kt_layout.isChecked():
            write_csv(os.path.join(out_dir, "kt_layout.csv"), kt_layout_rows)

        # Keep the two pivoting sweeps physically separate as well as
        # distinguishable by columns in the combined tables.  Random/Directed
        # output remains in the original top-level files.
        if self.out_parameters.isChecked():
            manifest_rows = [condition_metadata(condition) for condition in conditions]
            write_csv(os.path.join(out_dir, "condition_manifest.csv"), manifest_rows)

        family_dirs = {
            "pivoting_dmt": "Pivoting_DMT",
        }
        for family, family_dir_name in family_dirs.items():
            family_conditions = [
                condition.name for condition in conditions if condition.family == family
            ]
            if not family_conditions:
                continue
            family_dir = os.path.join(out_dir, family_dir_name)
            family_summary = [r for r in summary_rows if r.get("condition_family") == family]
            family_time = [r for r in time_rows if r.get("condition_family") == family]
            family_mt_total = [r for r in mt_total_rows if r.get("condition_family") == family]
            family_per_cell = [r for r in per_cell_rows if r.get("condition_family") == family]
            family_run_summary = [r for r in run_summary_rows if r.get("condition_family") == family]
            family_progression = [r for r in progression_rows if r.get("condition_family") == family]
            family_deflection = [r for r in deflection_rows if r.get("condition_family") == family]
            family_per_mt = [r for r in per_mt_rows if r.get("condition_family") == family]
            if self.out_time.isChecked():
                write_csv(os.path.join(family_dir, "completion_time_s.csv"), family_time)
            if self.out_mt_total.isChecked():
                write_csv(os.path.join(family_dir, "mt_total_emitted.csv"), family_mt_total)
            if self.out_summary.isChecked():
                write_csv(os.path.join(family_dir, "summary.csv"), family_summary)
                write_csv(os.path.join(family_dir, "run_summary.csv"), family_run_summary)
            if self.out_per_cell.isChecked():
                write_csv(os.path.join(family_dir, "per_cell_summary.csv"), family_per_cell)
            if self.out_progress.isChecked():
                write_csv(os.path.join(family_dir, "capture_progression.csv"), family_progression)
            if self.out_deflection.isChecked():
                write_csv(os.path.join(family_dir, "deflection_angles.csv"), family_deflection)
            if self.out_per_mt.isChecked():
                write_csv(os.path.join(family_dir, "per_mt_details.csv"), family_per_mt)
            if self.out_parameters.isChecked():
                family_parameter_row = dict(parameters_rows[0]) if parameters_rows else {}
                family_parameter_row.update({
                    "condition_family": family,
                    "export_conditions": ";".join(family_conditions),
                })
                write_csv(os.path.join(family_dir, "parameters_used.csv"), [family_parameter_row])
            if metric_only:
                generate_metric_only_superplot_ready_tables(
                    output_dir=family_dir,
                    condition_order=family_conditions,
                    summary_rows=family_summary,
                    include_time=self.out_time.isChecked(),
                    include_mt_total=self.out_mt_total.isChecked(),
                )
            else:
                generate_superplot_ready_tables(
                    output_dir=family_dir,
                    condition_order=family_conditions,
                    summary_rows=family_summary,
                    progression_rows=family_progression,
                    deflection_rows=family_deflection,
                )

        try:
            if metric_only:
                qc_lines = generate_metric_only_superplot_ready_tables(
                    output_dir=out_dir,
                    condition_order=condition_order,
                    summary_rows=summary_rows,
                    include_time=self.out_time.isChecked(),
                    include_mt_total=self.out_mt_total.isChecked(),
                )
            else:
                qc_lines = generate_superplot_ready_tables(
                    output_dir=out_dir,
                    condition_order=condition_order,
                    summary_rows=summary_rows,
                    progression_rows=progression_rows,
                    deflection_rows=deflection_rows,
                )
        except Exception as e:
            QMessageBox.warning(self, "Superplot table generation failed", str(e))
            return

        elapsed = time.perf_counter() - start_time
        progress.setValue(total_jobs)
        qc_text = "\n".join(qc_lines[-4:]) if qc_lines else "QC finished."
        QMessageBox.information(
            self,
            "Export complete",
            f"Exported selected results to:\n{out_dir}\n\n"
            f"Added superplot-ready tables in:\n{os.path.join(out_dir, 'Superplot_ready_tables')}\n\n"
            f"Parallel export QC:\n"
            f"- Workers used: {workers}\n"
            f"- Completed jobs: {completed_jobs}/{total_jobs}\n"
            f"- Elapsed compute/write time: {elapsed:.2f} s\n"
            f"- Backend: {'ProcessPoolExecutor' if workers > 1 else 'sequential'}\n\n"
            f"QC summary:\n{qc_text}"
        )


def main() -> None:
    if "--direct" in sys.argv[1:]:
        run_direct_from_cli(sys.argv[1:])
        return

    if not QT_AVAILABLE:
        print("PySide6 is required for GUI mode. Install it with: pip install PySide6")
        print("For direct output without GUI, run: python random_directed_microtubule_kinetochore_capture_3D_simulation_unified_engine.py --direct")
        sys.exit(1)

    app = QApplication(sys.argv)
    app.setFont(QFont(DISPLAY_FONT_FAMILY, 10))
    w = MainWindow()
    w.resize(1850, 900)
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
