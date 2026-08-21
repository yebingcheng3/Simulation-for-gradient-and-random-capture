#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MT Directed Capture Simulation v10 - Python version for Fig.5
Converted from the HTML/JavaScript model.

Core model:
- 2D single-centrosome / single-kinetochore simulation.
- MTs are emitted from the centrosome with directions capped to the KT-centered gradient zone.
- Optional KT-centered directional bias bends growing MT tips toward the KT inside the gradient radius.
- Optional dynamic instability: grow -> pause -> continue or catastrophe/shrink -> grow.
- At first KT capture, new MT emission stops, but already-emitted MTs continue to their own endpoint.
- Capture economy readouts use the true first-contact time and MT count at first contact.
- Theta-vs-frame uses each MT's own endpoint: capture frame for captured MTs, max-length/stuck frame for uncaptured MTs.
- Batch export includes exact plot tables, per-trial capture tables, and a full per-MT theta-vs-time CSV.

Dependencies:
- Required: numpy, matplotlib
- GUI: tkinter, matplotlib TkAgg backend
- XLSX export: openpyxl is optional; if unavailable, CSV files are exported instead.

Author: ChatGPT, generated for user workflow.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
import time
import threading
import queue
from dataclasses import dataclass, asdict, field
from pathlib import Path
from statistics import mean, median, pstdev
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np


# -----------------------------
# Data structures
# -----------------------------

@dataclass
class SimParams:
    # Sliders / UI defaults copied from the updated HTML program.
    mt_rate: float = 0.02                     # MTs emitted per second; confirmed default
    angle_spread_pi: float = 0.30             # half-angle spread, in units of pi; confirmed default
    kt_diam: float = 0.80                     # um; confirmed default
    grad_r: float = 2.50                      # um; confirmed default
    kt_dist: float = 6.70                     # um; confirmed default
    cycle_ms: float = 500.0                   # ms; confirmed default
    grow_dist: float = 0.110                  # um per growth step; confirmed default
    shrink_dist: float = 0.100                # um per catastrophe; confirmed default
    base_cat_percent: float = 20.0            # percent; confirmed default
    cont_bias: float = 0.20                   # 0-1; confirmed default
    max_len: float = 12.00                    # um; confirmed default
    max_time: float = 3000.0                 # s; confirmed default
    mech_bias: bool = True
    mech_inst: bool = True

    # Batch / output settings
    n_trials: int = 3000
    max_back: int = 40                        # frames before endpoint used for theta-vs-frame plot
    hist_bins: int = 25
    seed: Optional[int] = None                # None = true random sequence

    @property
    def kt_r(self) -> float:
        return self.kt_diam / 2.0

    @property
    def dt(self) -> float:
        return self.cycle_ms / 1000.0

    @property
    def base_cat(self) -> float:
        return self.base_cat_percent / 100.0

    @property
    def angle_spread_rad(self) -> float:
        return self.angle_spread_pi * math.pi


@dataclass
class Geometry:
    centro: Tuple[float, float]
    kt: Tuple[float, float]


@dataclass
class Microtubule:
    mt_id: int
    emitted_step: int
    emitted_time_s: float
    x: float
    y: float
    angle: float
    state: str = "grow"
    path_len: float = 0.0
    captured: bool = False
    capture_step: int = -1
    capture_frame: int = -1
    stuck_step: int = -1
    age_cycles: int = 0
    min_dist: float = math.inf
    min_dist_frame: int = -1
    path: List[Tuple[float, float]] = field(default_factory=list)
    theta_trace: List[float] = field(default_factory=list)
    theta_steps: List[int] = field(default_factory=list)                 # simulation step for each theta value
    theta_tip_positions: List[Tuple[float, float]] = field(default_factory=list)
    distance_to_kt_trace: List[float] = field(default_factory=list)

    def endpoint_frame(self) -> int:
        if self.captured:
            return self.capture_frame
        return len(self.theta_trace) - 1

    def final_theta(self) -> Optional[float]:
        if not self.theta_trace:
            return None
        return self.theta_trace[-1]


@dataclass
class TrialResult:
    condition: str
    trial: int
    captured: bool
    capture_time_s: float
    capture_step: int
    n_mt_at_capture: int
    first_captured_mt_id: Optional[int]
    stop_step: int
    stop_time_s: float
    total_mts_emitted: int
    mts: List[Microtubule]
    history: Optional[List[Dict[str, object]]] = None


@dataclass
class DemoMovie:
    title: str
    left_title: str
    right_title: str
    left_trial: TrialResult
    right_trial: TrialResult
    params: SimParams
    demo_seed: Optional[int] = None
    width_px: int = 3000
    height_px: int = 1500
    fps: float = 10.0
    scale_bar_um: float = 5.0

    @property
    def n_frames(self) -> int:
        left_n = len(self.left_trial.history or [])
        right_n = len(self.right_trial.history or [])
        # The comparison movie ends when the guided MT reaches the KT.
        # The unguided panel is shown at the same synchronized time.
        if self.right_trial.captured and right_n > 0:
            return right_n
        return max(left_n, right_n, 1)


@dataclass
class ConditionResult:
    condition: str
    trials: List[TrialResult]
    theta_cap_sum: List[float]
    theta_cap_n: List[int]
    theta_unc_sum: List[float]
    theta_unc_n: List[int]
    theta_tot_sum: List[float]
    theta_tot_n: List[int]

    def capture_times(self) -> List[float]:
        return [t.capture_time_s for t in self.trials if t.captured]

    def n_mts_at_capture(self) -> List[int]:
        return [t.n_mt_at_capture for t in self.trials if t.captured]

    def cap_count(self) -> int:
        return sum(1 for t in self.trials if t.captured)

    def theta_mean(self, kind: str) -> List[Optional[float]]:
        if kind == "captured":
            sums, ns = self.theta_cap_sum, self.theta_cap_n
        elif kind == "uncaptured":
            sums, ns = self.theta_unc_sum, self.theta_unc_n
        elif kind == "total":
            sums, ns = self.theta_tot_sum, self.theta_tot_n
        else:
            raise ValueError(f"Unknown theta kind: {kind}")
        return [s / n if n else None for s, n in zip(sums, ns)]

    def theta_n(self, kind: str) -> List[int]:
        if kind == "captured":
            return self.theta_cap_n
        if kind == "uncaptured":
            return self.theta_unc_n
        if kind == "total":
            return self.theta_tot_n
        raise ValueError(f"Unknown theta kind: {kind}")


@dataclass
class BatchResult:
    params: SimParams
    gradient: ConditionResult
    random_control: ConditionResult
    started_at: str
    elapsed_s: float


# -----------------------------
# Core math and simulation logic
# -----------------------------

def build_geometry(p: SimParams) -> Geometry:
    half = p.kt_dist / 2.0
    return Geometry(centro=(-half, 0.0), kt=(half, 0.0))


def normalize_angle(a: float) -> float:
    while a > math.pi:
        a -= 2 * math.pi
    while a < -math.pi:
        a += 2 * math.pi
    return a


def tangent_half_angle(origin: Tuple[float, float], target: Tuple[float, float], radius: float) -> float:
    d = math.hypot(target[0] - origin[0], target[1] - origin[1])
    if d <= radius:
        return math.pi / 2
    return math.asin(min(1.0, radius / d))


def theta_fig3(tip: Tuple[float, float], centro: Tuple[float, float], kt: Tuple[float, float], kt_r: float) -> float:
    """Fig.3-style deflection angle in degrees.

    Line 1: centrosome -> MT tip.
    Lines 2/3: tangents from centrosome to KT disc.
    theta = 0 if line 1 is inside the tangent cone; otherwise angle to nearest tangent.
    """
    dir_tip = math.atan2(tip[1] - centro[1], tip[0] - centro[0])
    dir_kt = math.atan2(kt[1] - centro[1], kt[0] - centro[0])
    half = tangent_half_angle(centro, kt, kt_r)
    d = normalize_angle(dir_tip - dir_kt)
    absd = abs(d)
    if absd <= half:
        return 0.0
    return (absd - half) * 180.0 / math.pi


def emission_angle(p: SimParams, geom: Geometry, rng: random.Random) -> float:
    base_a = math.atan2(geom.kt[1] - geom.centro[1], geom.kt[0] - geom.centro[0])
    cone_half = tangent_half_angle(geom.centro, geom.kt, p.grad_r)
    spread = min(p.angle_spread_rad, cone_half)
    return base_a + (rng.random() - 0.5) * 2.0 * spread


def make_mt(mt_id: int, emitted_step: int, p: SimParams, geom: Geometry, rng: random.Random) -> Microtubule:
    x, y = geom.centro
    mt = Microtubule(
        mt_id=mt_id,
        emitted_step=emitted_step,
        emitted_time_s=emitted_step * p.dt,
        x=x,
        y=y,
        angle=emission_angle(p, geom, rng),
        path=[(x, y)],
    )
    return mt


def retract_path(mt: Microtubule, amount: float) -> None:
    """Retract the MT tip by a true distance along the existing path."""
    remaining = amount
    while remaining > 1e-12 and len(mt.path) >= 2:
        tip = mt.path[-1]
        prev = mt.path[-2]
        seg_len = math.hypot(tip[0] - prev[0], tip[1] - prev[1])
        if seg_len <= remaining:
            mt.path.pop()
            remaining -= seg_len
        else:
            f = (seg_len - remaining) / seg_len
            nx = prev[0] + (tip[0] - prev[0]) * f
            ny = prev[1] + (tip[1] - prev[1]) * f
            mt.path[-1] = (nx, ny)
            remaining = 0.0

    mt.x, mt.y = mt.path[-1]
    mt.path_len = max(0.0, mt.path_len - amount)
    if len(mt.path) >= 2:
        a = mt.path[-2]
        b = mt.path[-1]
        mt.angle = math.atan2(b[1] - a[1], b[0] - a[0])


def kinetochore_lobes(center: Tuple[float, float], radius: float) -> List[Tuple[Tuple[float, float], float]]:
    """Return the exact two-circle geometry used for BOTH drawing and capture.

    Keeping capture geometry identical to the visible magenta symbol prevents
    premature capture before the MT tip actually touches the kinetochore.
    """
    lob_r = max(radius * 0.58, 0.09)
    offset = lob_r
    return [
        ((center[0], center[1] + offset), lob_r),
        ((center[0], center[1] - offset), lob_r),
    ]


def segment_circle_first_hit(
    start: Tuple[float, float],
    end: Tuple[float, float],
    center: Tuple[float, float],
    radius: float,
) -> Optional[Tuple[float, Tuple[float, float]]]:
    """Return the first segment-circle contact as (fraction, point), or None."""
    sx, sy = start
    ex, ey = end
    cx, cy = center
    dx = ex - sx
    dy = ey - sy
    fx = sx - cx
    fy = sy - cy

    # If numerical drift has already placed the starting point inside/on the KT,
    # count the current tip position as contact rather than advancing further.
    if fx * fx + fy * fy <= radius * radius + 1e-12:
        return 0.0, start

    a = dx * dx + dy * dy
    if a <= 1e-18:
        return None
    b = 2.0 * (fx * dx + fy * dy)
    c = fx * fx + fy * fy - radius * radius
    disc = b * b - 4.0 * a * c
    if disc < 0.0:
        return None

    root = math.sqrt(max(0.0, disc))
    t1 = (-b - root) / (2.0 * a)
    t2 = (-b + root) / (2.0 * a)
    candidates = [t for t in (t1, t2) if -1e-12 <= t <= 1.0 + 1e-12]
    if not candidates:
        return None
    t = max(0.0, min(1.0, min(candidates)))
    return t, (sx + dx * t, sy + dy * t)


def first_kinetochore_contact(
    start: Tuple[float, float],
    end: Tuple[float, float],
    kt_center: Tuple[float, float],
    kt_radius: float,
) -> Optional[Tuple[float, Tuple[float, float]]]:
    """First exact contact of an MT growth segment with either visible KT lobe."""
    hits = []
    for lobe_center, lobe_radius in kinetochore_lobes(kt_center, kt_radius):
        hit = segment_circle_first_hit(start, end, lobe_center, lobe_radius)
        if hit is not None:
            hits.append(hit)
    return min(hits, key=lambda item: item[0]) if hits else None


def step_mt(mt: Microtubule, with_gradient: bool, p: SimParams, geom: Geometry, rng: random.Random, current_step: int) -> bool:
    """Advance one MT by one DI cycle. Returns True if it captures during this step."""
    if mt.captured or mt.state == "stuck":
        return False

    mt.age_cycles += 1
    if mt.age_cycles > 50000:
        mt.state = "stuck"
        mt.stuck_step = current_step
        return False

    kt = geom.kt
    centro = geom.centro
    d_kt = math.hypot(mt.x - kt[0], mt.y - kt[1])
    in_grad = with_gradient and d_kt < p.grad_r

    if mt.state == "grow":
        if in_grad and p.mech_bias and p.cont_bias > 0:
            to_kt = math.atan2(kt[1] - mt.y, kt[0] - mt.x)
            d = normalize_angle(to_kt - mt.angle)
            mt.angle += d * p.cont_bias

        start = (mt.x, mt.y)
        proposed_end = (
            mt.x + math.cos(mt.angle) * p.grow_dist,
            mt.y + math.sin(mt.angle) * p.grow_dist,
        )
        contact = first_kinetochore_contact(start, proposed_end, kt, p.kt_r)

        if contact is not None:
            frac, contact_point = contact
            mt.x, mt.y = contact_point
            mt.path_len += p.grow_dist * frac
            mt.path.append(contact_point)
        else:
            mt.x, mt.y = proposed_end
            mt.path_len += p.grow_dist
            mt.path.append(proposed_end)

        th = theta_fig3((mt.x, mt.y), centro, kt, p.kt_r)
        mt.theta_trace.append(th)

        new_dist = math.hypot(mt.x - kt[0], mt.y - kt[1])
        mt.theta_steps.append(current_step)
        mt.theta_tip_positions.append((mt.x, mt.y))
        mt.distance_to_kt_trace.append(new_dist)
        if new_dist < mt.min_dist:
            mt.min_dist = new_dist
            mt.min_dist_frame = len(mt.theta_trace) - 1

        if contact is not None:
            mt.captured = True
            mt.capture_frame = len(mt.theta_trace) - 1
            mt.capture_step = current_step
            mt.state = "captured"
            return True

        mt.state = "pause"

    elif mt.state == "pause":
        if mt.path_len >= p.max_len:
            mt.state = "stuck"
            mt.stuck_step = current_step
            return False
        if not p.mech_inst:
            mt.state = "grow"
        else:
            mt.state = "shrink" if rng.random() < p.base_cat else "grow"

    elif mt.state == "shrink":
        retract_path(mt, p.shrink_dist)
        mt.state = "grow"

    return False


def snapshot_trial_state(mts: List[Microtubule], step: int, p: SimParams) -> Dict[str, object]:
    """Create a lightweight history snapshot for GUI animation/export."""
    items = []
    for mt in mts:
        items.append({
            "mt_id": mt.mt_id,
            "state": mt.state,
            "captured": mt.captured,
            "path": list(mt.path),
            "capture_step": mt.capture_step,
            "emitted_step": mt.emitted_step,
        })
    return {
        "step": step,
        "time_s": step * p.dt,
        "mts": items,
        "total_mts": len(mts),
        "captured_count": sum(1 for mt in mts if mt.captured),
    }


def run_trial(with_gradient: bool, p: SimParams, trial_index: int, rng: random.Random, record_history: bool = False) -> TrialResult:
    geom = build_geometry(p)
    max_steps = max(1, int(round(p.max_time / p.dt)))
    condition = "Gradient" if with_gradient else "Random"

    mts: List[Microtubule] = []
    step = 0
    emit_carry = 0.0
    captured = False
    capture_step = -1
    n_mt_at_capture = 0
    first_captured_mt_id: Optional[int] = None
    stop_emit = False
    next_mt_id = 1
    history: Optional[List[Dict[str, object]]] = [] if record_history else None

    if record_history:
        history.append(snapshot_trial_state(mts, 0, p))

    while step < max_steps:
        step += 1

        if not stop_emit:
            emit_carry += p.mt_rate * p.dt
            while emit_carry >= 1.0:
                mts.append(make_mt(next_mt_id, step, p, geom, rng))
                next_mt_id += 1
                emit_carry -= 1.0

        for mt in mts:
            if step_mt(mt, with_gradient, p, geom, rng, step):
                if not captured:
                    # True first-contact event.
                    # These values are the ONLY values used for:
                    #   1) time_to_capture_one_KT
                    #   2) MTs_emitted_to_capture_one_KT
                    captured = True
                    capture_step = step
                    n_mt_at_capture = len(mts)
                    first_captured_mt_id = mt.mt_id
                    stop_emit = True

        if record_history:
            history.append(snapshot_trial_state(mts, step, p))

        if captured:
            # Restored legacy theta-plot logic:
            # after first capture, stop emitting new MTs, but let every already-emitted MT
            # continue until it either captures the KT or reaches max length / stuck.
            # This gives each MT its own real endpoint for theta-vs-frame alignment.
            still_active = any((not mt.captured) and mt.state != "stuck" for mt in mts)
            if not still_active:
                break

    return TrialResult(
        condition=condition,
        trial=trial_index,
        captured=captured,
        capture_time_s=capture_step * p.dt if capture_step > 0 else p.max_time,
        capture_step=capture_step,
        n_mt_at_capture=n_mt_at_capture,
        first_captured_mt_id=first_captured_mt_id,
        stop_step=step,
        stop_time_s=step * p.dt,
        total_mts_emitted=len(mts),
        mts=mts,
        history=history,
    )


def collect_condition(with_gradient: bool, p: SimParams, rng: random.Random, progress_cb=None) -> ConditionResult:
    label = "Gradient" if with_gradient else "Random"
    max_back = p.max_back
    cap_sum = [0.0] * (max_back + 1)
    cap_n = [0] * (max_back + 1)
    unc_sum = [0.0] * (max_back + 1)
    unc_n = [0] * (max_back + 1)
    tot_sum = [0.0] * (max_back + 1)
    tot_n = [0] * (max_back + 1)
    trials: List[TrialResult] = []

    for i in range(1, p.n_trials + 1):
        r = run_trial(with_gradient, p, i, rng)
        trials.append(r)

        for mt in r.mts:
            trace = mt.theta_trace
            if not trace:
                continue
            f0 = mt.capture_frame if mt.captured else len(trace) - 1
            if f0 < 0:
                continue
            for k in range(max_back + 1):
                idx = f0 - k
                if idx < 0:
                    break
                value = trace[idx]
                if mt.captured:
                    cap_sum[k] += value
                    cap_n[k] += 1
                else:
                    unc_sum[k] += value
                    unc_n[k] += 1
                tot_sum[k] += value
                tot_n[k] += 1

        if progress_cb is not None:
            progress_cb(label, i, p.n_trials)

    return ConditionResult(label, trials, cap_sum, cap_n, unc_sum, unc_n, tot_sum, tot_n)


def run_batch(p: SimParams, progress_cb=None) -> BatchResult:
    started = time.strftime("%Y-%m-%d %H:%M:%S")
    t0 = time.time()
    rng = random.Random(p.seed) if p.seed is not None else random.Random()

    grad = collect_condition(True, p, rng, progress_cb=progress_cb)
    rand = collect_condition(False, p, rng, progress_cb=progress_cb)

    return BatchResult(
        params=p,
        gradient=grad,
        random_control=rand,
        started_at=started,
        elapsed_s=time.time() - t0,
    )


# -----------------------------
# Tables and summaries
# -----------------------------

def safe_mean(vals: List[float]) -> Optional[float]:
    return mean(vals) if vals else None


def safe_median(vals: List[float]) -> Optional[float]:
    return median(vals) if vals else None


def safe_sd(vals: List[float]) -> Optional[float]:
    return pstdev(vals) if vals else None


def summary_rows(batch: BatchResult) -> List[Dict[str, object]]:
    rows = []
    for cond in [batch.gradient, batch.random_control]:
        times = cond.capture_times()
        n_mts = cond.n_mts_at_capture()
        rows.append({
            "condition": cond.condition,
            "n_trials": len(cond.trials),
            "captured_trials": cond.cap_count(),
            "capture_fraction": cond.cap_count() / len(cond.trials) if cond.trials else None,
            "capture_time_mean_s": safe_mean(times),
            "capture_time_median_s": safe_median(times),
            "capture_time_sd_s": safe_sd(times),
            "mean_time_to_capture_one_KT_s": safe_mean(times),
            "median_time_to_capture_one_KT_s": safe_median(times),
            "sd_time_to_capture_one_KT_s": safe_sd(times),
            "mt_at_capture_mean": safe_mean(n_mts),
            "mt_at_capture_median": safe_median(n_mts),
            "mt_at_capture_sd": safe_sd([float(x) for x in n_mts]),
            "mean_MTs_emitted_to_capture_one_KT": safe_mean(n_mts),
            "median_MTs_emitted_to_capture_one_KT": safe_median(n_mts),
            "sd_MTs_emitted_to_capture_one_KT": safe_sd([float(x) for x in n_mts]),
            "mt_at_capture_min": min(n_mts) if n_mts else None,
            "mt_at_capture_max": max(n_mts) if n_mts else None,
            "theta_captured_MT_count": sum(1 for r in cond.trials for mt in r.mts if mt.theta_trace and mt.captured),
            "theta_uncaptured_MT_count": sum(1 for r in cond.trials for mt in r.mts if mt.theta_trace and not mt.captured),
            "theta_total_MT_count": sum(1 for r in cond.trials for mt in r.mts if mt.theta_trace),
            "theta_captured_n_at_leftmost_frame": cond.theta_n("captured")[batch.params.max_back],
            "theta_uncaptured_n_at_leftmost_frame": cond.theta_n("uncaptured")[batch.params.max_back],
            "theta_total_n_at_leftmost_frame": cond.theta_n("total")[batch.params.max_back],
        })
    return rows


def params_signature(p: SimParams) -> str:
    """Stable string signature for detecting whether GUI params match the last batch."""
    return json.dumps(asdict(p), ensure_ascii=False, sort_keys=True)


def compact_param_caption(p: SimParams) -> str:
    """Short parameter caption shown on plots to avoid mixing old plots with new GUI values."""
    seed_txt = "random" if p.seed is None else str(p.seed)
    return (
        f"rate={p.mt_rate:.2f}/s, spread={p.angle_spread_pi:.2f}π, KT={p.kt_diam:.2f} μm, "
        f"gradR={p.grad_r:.1f}, dist={p.kt_dist:.1f}, grow={p.grow_dist:.3f}, "
        f"shrink={p.shrink_dist:.3f}, cat={p.base_cat_percent:.1f}%, "
        f"bias={p.cont_bias:.2f}, maxL={p.max_len:.1f}, maxT={p.max_time:.0f}, seed={seed_txt}"
    )


def theta_endpoint_counts(batch: BatchResult) -> Dict[str, int]:
    """Count MTs contributing to each theta class across the batch."""
    out: Dict[str, int] = {}
    for label, cond in [("Gradient", batch.gradient), ("Random", batch.random_control)]:
        captured = 0
        uncaptured = 0
        for r in cond.trials:
            for mt in r.mts:
                if not mt.theta_trace:
                    continue
                if mt.captured:
                    captured += 1
                else:
                    uncaptured += 1
        out[f"{label}_captured"] = captured
        out[f"{label}_uncaptured"] = uncaptured
        out[f"{label}_total"] = captured + uncaptured
    return out


def theta_leftmost_counts(batch: BatchResult) -> Dict[str, int]:
    """n at the leftmost frame shown in theta plot, useful for QC of unstable curves."""
    k = batch.params.max_back
    out = {}
    for label, cond in [("Gradient", batch.gradient), ("Random", batch.random_control)]:
        out[f"{label}_captured_left_n"] = cond.theta_n("captured")[k]
        out[f"{label}_uncaptured_left_n"] = cond.theta_n("uncaptured")[k]
        out[f"{label}_total_left_n"] = cond.theta_n("total")[k]
    return out


def trial_summary_rows(batch: BatchResult) -> Iterable[Dict[str, object]]:
    for cond in [batch.gradient, batch.random_control]:
        for r in cond.trials:
            captured_mts = sum(1 for mt in r.mts if mt.captured)
            stuck_mts = sum(1 for mt in r.mts if mt.state == "stuck")
            yield {
                "condition": cond.condition,
                "trial": r.trial,
                "first_capture_occurred": r.captured,
                "first_capture_time_s": r.capture_time_s if r.captured else None,
                "first_capture_step": r.capture_step if r.captured else None,
                "n_mt_at_first_capture": r.n_mt_at_capture if r.captured else None,
                "first_captured_mt_id": r.first_captured_mt_id,
                "total_mts_emitted_at_trial_stop": r.total_mts_emitted,
                "captured_mt_count_at_trial_stop": captured_mts,
                "uncaptured_stuck_mt_count_at_trial_stop": stuck_mts,
                "trial_stop_reason": "all_emitted_MTs_reached_endpoint_after_first_capture" if r.captured else "max_time_no_capture",
                "post_first_capture_MT_continuation_for_theta": bool(r.captured),
                "trial_stop_step": r.stop_step,
                "trial_stop_time_s": r.stop_time_s,
            }



def per_trial_time_to_capture_one_kt_rows(batch: BatchResult) -> Iterable[Dict[str, object]]:
    """Per-trial table: how long it takes to capture one kinetochore.

    Because this Fig.5 model is a single-kinetochore simulation, each successful
    trial captures one kinetochore. Therefore first_capture_time_s is the
    time_to_capture_one_KT_s for that trial.
    """
    for cond in [batch.gradient, batch.random_control]:
        for r in cond.trials:
            yield {
                "condition": cond.condition,
                "trial": r.trial,
                "capture_occurred": r.captured,
                "time_to_capture_one_KT_s": r.capture_time_s if r.captured else None,
                "capture_step": r.capture_step if r.captured else None,
                "max_time_s": batch.params.max_time,
                "note": "First true MT-KT contact; not trial_stop_time" if r.captured else "No capture within max_time",
            }


def per_trial_mts_to_capture_one_kt_rows(batch: BatchResult) -> Iterable[Dict[str, object]]:
    """Per-trial table: how many MTs are emitted before capturing one kinetochore.

    Because this Fig.5 model is a single-kinetochore simulation, n_mt_at_capture
    is the number of MTs emitted to capture one kinetochore in each successful
    trial. total_mts_emitted is also included as a QC column.
    """
    for cond in [batch.gradient, batch.random_control]:
        for r in cond.trials:
            yield {
                "condition": cond.condition,
                "trial": r.trial,
                "capture_occurred": r.captured,
                "MTs_emitted_to_capture_one_KT": r.n_mt_at_capture if r.captured else None,
                "first_captured_mt_id": r.first_captured_mt_id if r.captured else None,
                "total_mts_emitted_at_trial_stop": r.total_mts_emitted,
                "note": "MT count at first true MT-KT contact; not total MT endpoint count" if r.captured else "No capture within max_time",
            }


def mt_endpoint_rows(batch: BatchResult) -> Iterable[Dict[str, object]]:
    for cond in [batch.gradient, batch.random_control]:
        for r in cond.trials:
            for mt in r.mts:
                end_frame = mt.endpoint_frame()
                yield {
                    "condition": cond.condition,
                    "trial": r.trial,
                    "mt_id": mt.mt_id,
                    "emitted_step": mt.emitted_step,
                    "emitted_time_s": mt.emitted_time_s,
                    "captured": mt.captured,
                    "final_state": mt.state,
                    "capture_step": mt.capture_step if mt.captured else None,
                    "capture_time_s": mt.capture_step * batch.params.dt if mt.captured else None,
                    "capture_frame_in_theta_trace": mt.capture_frame if mt.captured else None,
                    "endpoint_frame_in_theta_trace": end_frame,
                    "theta_trace_length": len(mt.theta_trace),
                    "initial_theta_deg": mt.theta_trace[0] if mt.theta_trace else None,
                    "final_theta_deg": mt.final_theta(),
                    "min_tip_to_KT_distance_um": mt.min_dist if math.isfinite(mt.min_dist) else None,
                    "min_distance_frame": mt.min_dist_frame if mt.min_dist_frame >= 0 else None,
                    "final_x_um": mt.x,
                    "final_y_um": mt.y,
                    "final_path_length_um": mt.path_len,
                    "age_cycles": mt.age_cycles,
                }


def theta_vs_frame_rows(batch: BatchResult) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for k in range(batch.params.max_back, -1, -1):
        row: Dict[str, object] = {"relative_frame": -k}
        for cond_key, cond in [("Gradient", batch.gradient), ("Random", batch.random_control)]:
            for kind in ["captured", "uncaptured", "total"]:
                m = cond.theta_mean(kind)[k]
                n = cond.theta_n(kind)[k]
                row[f"{cond_key}_{kind}_mean_theta_deg"] = m
                row[f"{cond_key}_{kind}_n"] = n
        rows.append(row)
    return rows


def theta_trace_aligned_rows(batch: BatchResult) -> Iterable[Dict[str, object]]:
    """Long table with the exact aligned MT-level values used to calculate theta-vs-frame means.

    Each row is one MT at one relative frame. This is not the full unaligned trace; it is the
    plot-specific table for frame -max_back ... 0 aligned to capture or max-length endpoint.
    """
    max_back = batch.params.max_back
    for cond in [batch.gradient, batch.random_control]:
        for r in cond.trials:
            for mt in r.mts:
                trace = mt.theta_trace
                if not trace:
                    continue
                f0 = mt.capture_frame if mt.captured else len(trace) - 1
                if f0 < 0:
                    continue
                endpoint_type = "capture" if mt.captured else "max_length_or_stuck"
                for k in range(max_back, -1, -1):
                    idx = f0 - k
                    if idx < 0:
                        continue
                    yield {
                        "condition": cond.condition,
                        "trial": r.trial,
                        "mt_id": mt.mt_id,
                        "captured": mt.captured,
                        "endpoint_type": endpoint_type,
                        "relative_frame": -k,
                        "theta_trace_index": idx,
                        "theta_deg": trace[idx],
                    }


def mt_theta_time_full_rows(batch: BatchResult) -> Iterable[Dict[str, object]]:
    """Full unaligned MT-level theta-over-time table.

    Each row is one theta measurement for one MT at one simulated time point.
    This table is different from Theta_vs_frame_plot_data:
    - Theta_vs_frame_plot_data is condition-level mean aligned to endpoint.
    - MT_theta_time_full is the raw per-MT time series.
    """
    p = batch.params
    geom = build_geometry(p)
    for cond in [batch.gradient, batch.random_control]:
        for r in cond.trials:
            for mt in r.mts:
                endpoint_type = "capture" if mt.captured else "max_length_or_stuck"
                for idx, theta in enumerate(mt.theta_trace):
                    step = mt.theta_steps[idx] if idx < len(mt.theta_steps) else None
                    tip = mt.theta_tip_positions[idx] if idx < len(mt.theta_tip_positions) else (None, None)
                    dist = mt.distance_to_kt_trace[idx] if idx < len(mt.distance_to_kt_trace) else None
                    abs_time = step * p.dt if step is not None else None
                    rel_time = abs_time - mt.emitted_time_s if abs_time is not None else None
                    yield {
                        "condition": cond.condition,
                        "trial": r.trial,
                        "mt_id": mt.mt_id,
                        "captured": mt.captured,
                        "endpoint_type": endpoint_type,
                        "is_first_captured_mt": bool(r.first_captured_mt_id == mt.mt_id),
                        "trace_index": idx,
                        "simulation_step": step,
                        "absolute_time_s": abs_time,
                        "time_since_MT_emission_s": rel_time,
                        "emitted_step": mt.emitted_step,
                        "emitted_time_s": mt.emitted_time_s,
                        "capture_step": mt.capture_step if mt.captured else None,
                        "capture_time_s": mt.capture_step * p.dt if mt.captured else None,
                        "tip_x_um": tip[0],
                        "tip_y_um": tip[1],
                        "centrosome_x_um": geom.centro[0],
                        "centrosome_y_um": geom.centro[1],
                        "kt_x_um": geom.kt[0],
                        "kt_y_um": geom.kt[1],
                        "tip_to_KT_center_distance_um": dist,
                        "theta_boundary_corrected_deg": theta,
                    }


def capture_histogram_rows(batch: BatchResult) -> List[Dict[str, object]]:
    g_times = batch.gradient.capture_times()
    r_times = batch.random_control.capture_times()
    all_times = g_times + r_times
    if not all_times:
        return []
    max_t = max(all_times) * 1.1
    n_bins = max(1, batch.params.hist_bins)
    bin_w = max_t / n_bins if max_t > 0 else 1.0
    g_counts = [0] * n_bins
    r_counts = [0] * n_bins
    for t in g_times:
        idx = min(int(t // bin_w), n_bins - 1)
        g_counts[idx] += 1
    for t in r_times:
        idx = min(int(t // bin_w), n_bins - 1)
        r_counts[idx] += 1
    rows = []
    for i in range(n_bins):
        left = i * bin_w
        right = (i + 1) * bin_w
        rows.append({
            "bin_index": i + 1,
            "bin_left_s": left,
            "bin_right_s": right,
            "bin_mid_s": (left + right) / 2.0,
            "Gradient_count": g_counts[i],
            "Random_count": r_counts[i],
        })
    return rows


def parameter_rows(p: SimParams) -> List[Dict[str, object]]:
    rows = []
    labels = {
        "mt_rate": "MTs emitted per second",
        "angle_spread_pi": "Angle spread (x pi)",
        "kt_diam": "KT diameter (um)",
        "grad_r": "Gradient radius (um)",
        "kt_dist": "KT-centrosome distance (um)",
        "cycle_ms": "DI cycle interval (ms)",
        "grow_dist": "Grow per cycle (um)",
        "shrink_dist": "Shrink per cycle (um)",
        "base_cat_percent": "Base catastrophe (%)",
        "cont_bias": "Continuous bias (0-1)",
        "max_len": "Max MT length (um)",
        "max_time": "Max time (s)",
        "mech_bias": "Directional bias switch",
        "mech_inst": "Dynamic instability switch",
        "n_trials": "Trials per condition",
        "max_back": "Theta plot frames before endpoint",
        "hist_bins": "Capture-time histogram bins",
        "seed": "Random seed; blank/None means non-fixed random run",
    }
    d = asdict(p)
    for key, value in d.items():
        rows.append({"parameter_id": key, "label": labels.get(key, key), "value": value})
    rows.append({"parameter_id": "dt", "label": "Time per DI cycle (s)", "value": p.dt})
    rows.append({"parameter_id": "kt_r", "label": "KT radius (um)", "value": p.kt_r})
    rows.append({"parameter_id": "base_cat", "label": "Base catastrophe probability", "value": p.base_cat})
    rows.append({"parameter_id": "angle_spread_rad", "label": "Angle spread (radian)", "value": p.angle_spread_rad})
    rows.append({
        "parameter_id": "trial_endpoint_rule",
        "label": "Trial endpoint rule",
        "value": "At first KT capture, stop emitting new MTs; already-emitted MTs continue to capture or max length/stuck. Capture economy uses first_capture_time_s and n_mt_at_first_capture only.",
    })
    return rows


# -----------------------------
# Export helpers
# -----------------------------

def write_csv(path: Path, rows: Iterable[Dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8-sig")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def autosize_openpyxl(ws, max_width: int = 38) -> None:
    for column_cells in ws.columns:
        max_len = 0
        letter = column_cells[0].column_letter
        for cell in column_cells:
            try:
                value = "" if cell.value is None else str(cell.value)
            except Exception:
                value = ""
            max_len = max(max_len, len(value))
        ws.column_dimensions[letter].width = min(max(max_len + 2, 10), max_width)


def write_xlsx(path: Path, sheets: Dict[str, List[Dict[str, object]]]) -> bool:
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
        from openpyxl.utils import get_column_letter
    except Exception:
        return False

    wb = Workbook()
    default_ws = wb.active
    wb.remove(default_ws)
    header_fill = PatternFill("solid", fgColor="D9EAD3")
    header_font = Font(bold=True)
    thin = Side(style="thin", color="DDDDDD")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    for sheet_name, rows in sheets.items():
        ws = wb.create_sheet(sheet_name[:31])
        rows = list(rows)
        if not rows:
            ws.append(["No data"])
            continue
        headers = list(rows[0].keys())
        ws.append(headers)
        for cell in ws[1]:
            cell.fill = header_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            cell.border = border
        for row in rows:
            ws.append([row.get(h) for h in headers])
        for row in ws.iter_rows():
            for cell in row:
                cell.border = border
                cell.alignment = Alignment(vertical="center")
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        autosize_openpyxl(ws)

    wb.save(path)
    return True


def current_run_label(p: SimParams) -> str:
    """Label the manually selected condition based on current checkboxes.

    The user manually controls Directional bias and Dynamic instability.
    Export therefore reports only the current selected run, not an automatic ON/OFF pair.
    """
    bias = "BiasON" if p.mech_bias else "BiasOFF"
    di = "DI_ON" if p.mech_inst else "DI_OFF"
    return f"{bias}_{di}"


def _capture_pairs_for_condition(cond: ConditionResult) -> List[Tuple[int, float]]:
    """Successful capture rows for one condition: (MT count, first-contact time)."""
    pairs: List[Tuple[int, float]] = []
    for r in cond.trials:
        if not r.captured:
            continue
        pairs.append((r.n_mt_at_capture, r.capture_time_s))
    return pairs


def _four_column_capture_rows(
    random_pairs: List[Tuple[int, float]],
    gradient_pairs: List[Tuple[int, float]],
) -> List[Dict[str, object]]:
    """Put Random on the left two columns and Gradient on the right two columns.

    Rows are aligned by simulation result index only. They are not interpreted as
    paired biological repeats. If one condition has fewer successful captures,
    the remaining cells for that condition are left blank.
    """
    rows: List[Dict[str, object]] = []
    max_n = max(len(random_pairs), len(gradient_pairs)) if (random_pairs or gradient_pairs) else 0
    for i in range(max_n):
        random_mt, random_time = random_pairs[i] if i < len(random_pairs) else (None, None)
        gradient_mt, gradient_time = gradient_pairs[i] if i < len(gradient_pairs) else (None, None)
        rows.append({
            "Random_MTs_emitted_to_capture_one_KT": random_mt,
            "Random_Time_to_capture_one_KT_s": random_time,
            "Gradient_MTs_emitted_to_capture_one_KT": gradient_mt,
            "Gradient_Time_to_capture_one_KT_s": gradient_time,
        })
    return rows


def article_capture_full_rows(batch: BatchResult) -> List[Dict[str, object]]:
    """Article-ready capture economy table with explicit Random vs Gradient columns.

    Column order requested by the user:
    1) no-gradient/random MT count
    2) no-gradient/random capture time
    3) gradient MT count
    4) gradient capture time
    """
    random_pairs = _capture_pairs_for_condition(batch.random_control)
    gradient_pairs = _capture_pairs_for_condition(batch.gradient)
    return _four_column_capture_rows(random_pairs, gradient_pairs)


def article_capture_sample_rows(batch: BatchResult, n_sample: int = 500) -> List[Dict[str, object]]:
    """Sample up to n_sample successful capture rows per condition for display points."""
    random_pairs = _capture_pairs_for_condition(batch.random_control)
    gradient_pairs = _capture_pairs_for_condition(batch.gradient)

    # If a seed is set, sampling is reproducible. If seed is blank, it is random.
    rng = random.Random(batch.params.seed) if batch.params.seed is not None else random.Random()

    def sample_pairs(pairs: List[Tuple[int, float]]) -> List[Tuple[int, float]]:
        if len(pairs) <= n_sample:
            return list(pairs)
        idxs = sorted(rng.sample(range(len(pairs)), n_sample))
        return [pairs[i] for i in idxs]

    return _four_column_capture_rows(sample_pairs(random_pairs), sample_pairs(gradient_pairs))


def article_theta_mean_curve_rows(batch: BatchResult) -> List[Dict[str, object]]:
    """Article-ready theta curve: six mean-value series matching the GUI plot.

    This sheet intentionally contains ONLY final mean curves, not per-MT traces.
    Columns correspond to the six plotted lines:
    Gradient/Random × captured/uncaptured/total.
    """
    rows: List[Dict[str, object]] = []
    g_cap = batch.gradient.theta_mean("captured")
    g_unc = batch.gradient.theta_mean("uncaptured")
    g_tot = batch.gradient.theta_mean("total")
    r_cap = batch.random_control.theta_mean("captured")
    r_unc = batch.random_control.theta_mean("uncaptured")
    r_tot = batch.random_control.theta_mean("total")

    for k in range(batch.params.max_back, -1, -1):
        rows.append({
            "relative_frame": -k,
            "Gradient_captured_mean_theta_deg": g_cap[k],
            "Gradient_uncaptured_mean_theta_deg": g_unc[k],
            "Gradient_total_mean_theta_deg": g_tot[k],
            "Random_captured_mean_theta_deg": r_cap[k],
            "Random_uncaptured_mean_theta_deg": r_unc[k],
            "Random_total_mean_theta_deg": r_tot[k],
        })
    return rows


def article_qc_summary_rows(batch: BatchResult, n_sample: int = 500) -> List[Dict[str, object]]:
    """Small QC table for article-ready export, split by Random and Gradient."""
    random_pairs = _capture_pairs_for_condition(batch.random_control)
    gradient_pairs = _capture_pairs_for_condition(batch.gradient)
    full = article_capture_full_rows(batch)
    sample = article_capture_sample_rows(batch, n_sample=n_sample)
    theta_counts = theta_endpoint_counts(batch)
    theta_left = theta_leftmost_counts(batch)

    def cond_values(prefix: str, cond: ConditionResult, pairs: List[Tuple[int, float]]) -> List[Tuple[str, object]]:
        times = [t for _, t in pairs]
        mts = [m for m, _ in pairs]
        return [
            (f"{prefix}_n_trials", len(cond.trials)),
            (f"{prefix}_successful_captures", len(pairs)),
            (f"{prefix}_failed_captures", len(cond.trials) - len(pairs)),
            (f"{prefix}_capture_fraction", len(pairs) / len(cond.trials) if cond.trials else None),
            (f"{prefix}_mean_MTs_emitted_to_capture_one_KT", safe_mean(mts)),
            (f"{prefix}_median_MTs_emitted_to_capture_one_KT", safe_median(mts)),
            (f"{prefix}_sd_MTs_emitted_to_capture_one_KT", safe_sd([float(x) for x in mts])),
            (f"{prefix}_mean_Time_to_capture_one_KT_s", safe_mean(times)),
            (f"{prefix}_median_Time_to_capture_one_KT_s", safe_median(times)),
            (f"{prefix}_sd_Time_to_capture_one_KT_s", safe_sd(times)),
        ]

    values: List[Tuple[str, object]] = [
        ("run_label", current_run_label(batch.params)),
        ("Directional_bias", bool(batch.params.mech_bias)),
        ("Dynamic_instability", bool(batch.params.mech_inst)),
        ("Capture_full_rows", len(full)),
        ("Capture_sample_requested_n_per_condition", n_sample),
        ("Capture_sample_rows", len(sample)),
        ("Capture_full_column_order", "Random_MTs, Random_Time_s, Gradient_MTs, Gradient_Time_s"),
    ]
    values.extend(cond_values("Random", batch.random_control, random_pairs))
    values.extend(cond_values("Gradient", batch.gradient, gradient_pairs))
    values.extend([
        ("Random_theta_captured_MT_count", theta_counts.get("Random_captured")),
        ("Random_theta_uncaptured_MT_count", theta_counts.get("Random_uncaptured")),
        ("Random_theta_total_MT_count", theta_counts.get("Random_total")),
        ("Gradient_theta_captured_MT_count", theta_counts.get("Gradient_captured")),
        ("Gradient_theta_uncaptured_MT_count", theta_counts.get("Gradient_uncaptured")),
        ("Gradient_theta_total_MT_count", theta_counts.get("Gradient_total")),
        ("Random_theta_captured_n_at_leftmost_frame", theta_left.get("Random_captured_left_n")),
        ("Random_theta_uncaptured_n_at_leftmost_frame", theta_left.get("Random_uncaptured_left_n")),
        ("Random_theta_total_n_at_leftmost_frame", theta_left.get("Random_total_left_n")),
        ("Gradient_theta_captured_n_at_leftmost_frame", theta_left.get("Gradient_captured_left_n")),
        ("Gradient_theta_uncaptured_n_at_leftmost_frame", theta_left.get("Gradient_uncaptured_left_n")),
        ("Gradient_theta_total_n_at_leftmost_frame", theta_left.get("Gradient_total_left_n")),
        ("note_capture_metrics", "Capture metrics use true first MT-KT contact only, not trial_stop_time."),
        ("note_capture_row_alignment", "Random and Gradient are placed side-by-side by result index for plotting convenience; rows are not paired biological repeats."),
        ("note_theta_metrics", "Theta curve uses each MT endpoint after already-emitted MTs continue to capture or max length/stuck."),
    ])
    return [{"item": k, "value": v} for k, v in values]


def article_parameter_rows(p: SimParams) -> List[Dict[str, object]]:
    """Simplified parameter table for article export."""
    items = [
        ("run_label", current_run_label(p)),
        ("MTs_emitted_per_second", p.mt_rate),
        ("Angle_spread_x_pi", p.angle_spread_pi),
        ("KT_diameter_um", p.kt_diam),
        ("Gradient_radius_um", p.grad_r),
        ("KT_centrosome_distance_um", p.kt_dist),
        ("DI_cycle_interval_ms", p.cycle_ms),
        ("Grow_per_cycle_um", p.grow_dist),
        ("Shrink_per_cycle_um", p.shrink_dist),
        ("Base_catastrophe_percent", p.base_cat_percent),
        ("Continuous_bias", p.cont_bias),
        ("Max_MT_length_um", p.max_len),
        ("Max_time_s", p.max_time),
        ("Directional_bias", p.mech_bias),
        ("Dynamic_instability", p.mech_inst),
        ("Trials_per_condition", p.n_trials),
        ("Theta_plot_frames_before_endpoint", p.max_back),
        ("Seed", "blank/random" if p.seed is None else p.seed),
    ]
    return [{"parameter": k, "value": v} for k, v in items]


def build_article_ready_export_tables(batch: BatchResult, n_sample: int = 500) -> Dict[str, List[Dict[str, object]]]:
    """Minimal publication-plot-ready workbook tables.

    Sheets:
    1. Capture_full: all successful trials, four columns: Random MT/time and Gradient MT/time.
    2. Capture_sample500: up to 500 randomly sampled successful trials per condition, same four-column format.
    3. Theta_mean_curve: final six-line mean theta curve by relative frame.
    4. QC_summary: compact verification.
    5. Parameters: current manually selected settings.
    """
    return {
        "Capture_full": article_capture_full_rows(batch),
        "Capture_sample500": article_capture_sample_rows(batch, n_sample=n_sample),
        "Theta_mean_curve": article_theta_mean_curve_rows(batch),
        "QC_summary": article_qc_summary_rows(batch, n_sample=n_sample),
        "Parameters": article_parameter_rows(batch.params),
    }


def build_export_tables(batch: BatchResult, include_long_trace: bool = True) -> Dict[str, List[Dict[str, object]]]:
    sheets: Dict[str, List[Dict[str, object]]] = {
        "Parameters": parameter_rows(batch.params),
        "Summary": summary_rows(batch),
        "Trial_summary": list(trial_summary_rows(batch)),
        "Time_to_capture_one_KT": list(per_trial_time_to_capture_one_kt_rows(batch)),
        "MTs_to_capture_one_KT": list(per_trial_mts_to_capture_one_kt_rows(batch)),
        "Capture_time_histogram": capture_histogram_rows(batch),
        "Theta_vs_frame_plot_data": theta_vs_frame_rows(batch),
        "MT_endpoint_summary": list(mt_endpoint_rows(batch)),
        "Batch_QC": [{
            "params_signature": params_signature(batch.params),
            "parameter_caption": compact_param_caption(batch.params),
            **theta_endpoint_counts(batch),
            **theta_leftmost_counts(batch),
        }],
        "MT_theta_time_full_README": [{
            "table": "MT_theta_time_full.csv",
            "description": "Full per-microtubule theta-vs-time data are exported as a separate CSV file to avoid freezing Excel. After first capture, already-emitted MTs continue to their own endpoint for theta statistics.",
            "theta_definition": "theta_boundary_corrected_deg = max(|angle(centrosome->tip) - angle(centrosome->KT)| - asin(KT_radius / centrosome_to_KT_distance), 0), in degrees.",
        }],
    }
    if include_long_trace:
        sheets["Theta_trace_aligned_long"] = list(theta_trace_aligned_rows(batch))
    return sheets


def export_batch(batch: BatchResult, out_dir: Path, prefix: str = "Fig5_article_ready_simulation", include_long_trace: bool = False) -> Dict[str, str]:
    """Export article-ready data only.

    The export is intentionally minimal:
    - no per-MT theta trace CSV
    - no endpoint debugging sheets
    - no automatic ON/OFF pairing
    - only the current manually selected mechanism settings are exported
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    label = current_run_label(batch.params)
    base = f"{prefix}_{label}_{timestamp}"

    n_sample = 500
    tables = build_article_ready_export_tables(batch, n_sample=n_sample)
    xlsx_path = out_dir / f"{base}.xlsx"

    exported: Dict[str, str] = {}
    if write_xlsx(xlsx_path, tables):
        exported["article_ready_xlsx"] = str(xlsx_path)
    else:
        csv_dir = out_dir / f"{base}_csv_tables"
        csv_dir.mkdir(parents=True, exist_ok=True)
        for name, rows in tables.items():
            path = csv_dir / f"{name}.csv"
            write_csv(path, rows)
        exported["article_ready_csv_dir"] = str(csv_dir)

    # Compact metadata for reproducibility, not for plotting.
    meta_path = out_dir / f"{base}_metadata.json"
    meta = {
        "started_at": batch.started_at,
        "elapsed_s": batch.elapsed_s,
        "run_label": label,
        "params": asdict(batch.params),
        "qc_summary": article_qc_summary_rows(batch, n_sample=n_sample),
        "export_note": (
            "Article-ready export. Capture_full and Capture_sample500 contain four columns: "
            "Random MT count/time on the left and Gradient MT count/time on the right. "
            "Theta_mean_curve contains only mean theta values by relative frame."
        ),
    }
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    exported["metadata_json"] = str(meta_path)

    return exported


# -----------------------------
# Plotting
# -----------------------------

def ensure_matplotlib_backend(gui: bool = False) -> None:
    import matplotlib
    if not gui:
        matplotlib.use("Agg", force=True)

    # Keep exported SVG/PDF as real editable vector text whenever possible.
    # Without these settings, matplotlib may convert text to paths or embed
    # Type 3 fonts, which can cause editing/export issues in Illustrator.
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["font.family"] = ["Arial", "DejaVu Sans"]


def plot_capture_histogram(batch: BatchResult, save_path: Optional[Path] = None, ax=None):
    import matplotlib.pyplot as plt

    rows = capture_histogram_rows(batch)
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        created_fig = True
    else:
        fig = ax.figure
    ax.clear()

    if not rows:
        ax.text(0.5, 0.5, "No captures within time limit", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    else:
        mids = [r["bin_mid_s"] for r in rows]
        width = rows[0]["bin_right_s"] - rows[0]["bin_left_s"] if rows else 1.0
        g_counts = [r["Gradient_count"] for r in rows]
        r_counts = [r["Random_count"] for r in rows]
        ax.bar([m - width * 0.18 for m in mids], r_counts, width=width * 0.36, label="Random", alpha=0.55, color="#D85A30")
        ax.bar([m + width * 0.18 for m in mids], g_counts, width=width * 0.36, label="Gradient", alpha=0.65, color="#1D9E75")
        ax.set_xlabel("Capture time (s)")
        ax.set_ylabel("Count")
        ax.set_title("Capture-time distribution")
        ax.legend(frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        if created_fig:
            plt.close(fig)
    return fig, ax


def plot_theta_vs_frame(batch: BatchResult, save_path: Optional[Path] = None, ax=None):
    import matplotlib.pyplot as plt

    rows = theta_vs_frame_rows(batch)
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.8, 4.6))
        created_fig = True
    else:
        fig = ax.figure
    ax.clear()

    rel = [r["relative_frame"] for r in rows]
    series = [
        ("Gradient captured", "Gradient_captured_mean_theta_deg", "#1D9E75", "-"),
        ("Gradient uncaptured", "Gradient_uncaptured_mean_theta_deg", "#1D9E75", "--"),
        ("Gradient total", "Gradient_total_mean_theta_deg", "#1D9E75", ":"),
        ("Random captured", "Random_captured_mean_theta_deg", "#D85A30", "-"),
        ("Random uncaptured", "Random_uncaptured_mean_theta_deg", "#D85A30", "--"),
        ("Random total", "Random_total_mean_theta_deg", "#D85A30", ":"),
    ]
    counts = theta_endpoint_counts(batch)
    left_counts = theta_leftmost_counts(batch)
    label_count_key = {
        "Gradient captured": "Gradient_captured",
        "Gradient uncaptured": "Gradient_uncaptured",
        "Gradient total": "Gradient_total",
        "Random captured": "Random_captured",
        "Random uncaptured": "Random_uncaptured",
        "Random total": "Random_total",
    }
    left_count_key = {
        "Gradient captured": "Gradient_captured_left_n",
        "Gradient uncaptured": "Gradient_uncaptured_left_n",
        "Gradient total": "Gradient_total_left_n",
        "Random captured": "Random_captured_left_n",
        "Random uncaptured": "Random_uncaptured_left_n",
        "Random total": "Random_total_left_n",
    }
    for label, key, color, linestyle in series:
        y = [r[key] for r in rows]
        count_key = label_count_key.get(label)
        left_key = left_count_key.get(label)
        display_label = label
        if count_key is not None:
            display_label = f"{label} (n={counts.get(count_key, 0)}, left={left_counts.get(left_key, 0)})"
        ax.plot(rel, y, label=display_label, color=color, linestyle=linestyle, linewidth=2)
    ax.set_xlabel("Frame relative to endpoint (0 = capture or max length)")
    ax.set_ylabel("Deflection angle θ (deg)")
    ax.set_title("θ vs frame aligned to endpoint\n" + compact_param_caption(batch.params), fontsize=10)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        if created_fig:
            plt.close(fig)
    return fig, ax


SCENE_BG = "#DADADA"
CENTROSOME_GREEN = "#0F6B57"
MT_GREEN = "#1FA33C"
KT_MAGENTA = "#F414DD"
GREY_TRACE = "#9B9B9B"
GRADIENT_STOPS = [
    (0.00, "#ff4a0c", 1.00),
    (0.40, "#ff990c", 1.00),
    (0.70, "#ffd40c", 0.85),
    (1.00, "#fff4ac", 0.00),
]
MOVIE_S12_TITLE = "Movie S12: Simulation of microtubule–kinetochore capture with and without a guidance gradient"


def _hex_to_rgb01(hex_color: str):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


def _lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def _gradient_color_at(position: float):
    p = max(0.0, min(1.0, float(position)))
    for (p0, c0, a0), (p1, c1, a1) in zip(GRADIENT_STOPS[:-1], GRADIENT_STOPS[1:]):
        if p <= p1:
            t = 0.0 if p1 == p0 else (p - p0) / (p1 - p0)
            rgb0 = _hex_to_rgb01(c0)
            rgb1 = _hex_to_rgb01(c1)
            rgb = tuple(_lerp(rgb0[i], rgb1[i], t) for i in range(3))
            a = _lerp(a0, a1, t)
            return (*rgb, a)
    rgb = _hex_to_rgb01(GRADIENT_STOPS[-1][1])
    a = GRADIENT_STOPS[-1][2]
    return (*rgb, a)


def normalize_export_save_path(path: str, ext: str) -> Path:
    """Return a writable path with the requested file extension.

    Tk's asksaveasfilename can return a path without an extension or with an
    unintended extension when the user switches to "All files". This helper
    prevents silent saves such as a PDF button producing a file without .pdf.
    """
    clean_ext = "." + ext.lower().lstrip(".")
    save_path = Path(path)
    if save_path.suffix.lower() != clean_ext:
        save_path = save_path.with_suffix(clean_ext)
    return save_path


def save_figure_robust(fig, save_path: Path) -> None:
    """Save a matplotlib figure with settings suitable for SVG/PDF/PNG export."""
    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)
    ext = save_path.suffix.lower().lstrip(".")
    kwargs = {
        "bbox_inches": "tight",
        "facecolor": fig.get_facecolor(),
    }
    if ext in {"png", "jpg", "jpeg", "tif", "tiff"}:
        kwargs["dpi"] = 600
    fig.savefig(str(save_path), **kwargs)


def _draw_soft_gradient(ax, geom: Geometry, p: SimParams) -> None:
    """Draw the orange-yellow guidance gradient matching the AI reference."""
    from matplotlib.patches import Circle

    n_rings = 36
    for i in range(n_rings, 0, -1):
        frac = i / n_rings
        radius = max(p.grad_r * frac, 1e-6)
        rgba = _gradient_color_at(frac)
        ax.add_patch(
            Circle(
                geom.kt,
                radius,
                facecolor=rgba,
                edgecolor="none",
                zorder=0,
            )
        )


def _draw_kinetochore_symbol(ax, center: Tuple[float, float], radius: float) -> None:
    from matplotlib.patches import Circle
    for lobe_center, lobe_radius in kinetochore_lobes(center, radius):
        ax.add_patch(Circle(lobe_center, lobe_radius, facecolor=KT_MAGENTA, edgecolor="#8A007D", linewidth=0.7, zorder=6))


def draw_scene_frame(trial: TrialResult, p: SimParams, with_gradient: bool, frame_index: Optional[int] = None, save_path: Optional[Path] = None, ax=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    geom = build_geometry(p)
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=160)
        created_fig = True
    else:
        fig = ax.figure
    ax.clear()
    ax.set_facecolor(SCENE_BG)
    fig.patch.set_facecolor(SCENE_BG)

    snapshot = None
    if trial.history:
        idx = len(trial.history) - 1 if frame_index is None else max(0, min(frame_index, len(trial.history) - 1))
        snapshot = trial.history[idx]

    if with_gradient:
        _draw_soft_gradient(ax, geom, p)

    ax.plot([geom.centro[0], geom.kt[0]], [geom.centro[1], geom.kt[1]], linestyle="--", color="#C3C3C3", linewidth=1.0, zorder=1)

    if snapshot is not None:
        mts_to_draw = snapshot["mts"]
        time_s = float(snapshot["time_s"])
        total_mts = int(snapshot["total_mts"])
    else:
        mts_to_draw = [{"path": list(mt.path), "state": mt.state, "captured": mt.captured, "mt_id": mt.mt_id} for mt in trial.mts]
        time_s = trial.stop_time_s
        total_mts = trial.total_mts_emitted

    for mt in mts_to_draw:
        path = mt.get("path", [])
        if len(path) < 2:
            continue
        x = [pt[0] for pt in path]
        y = [pt[1] for pt in path]
        state = mt.get("state", "grow")
        captured = bool(mt.get("captured", False))
        if captured:
            color = MT_GREEN
            lw = 2.0
            alpha = 0.98
        elif state == "stuck":
            color = GREY_TRACE
            lw = 1.0
            alpha = 0.45
        else:
            color = MT_GREEN
            lw = 1.7
            alpha = 0.90
        ax.plot(x, y, color=color, linewidth=lw, alpha=alpha, solid_capstyle="round", zorder=3)

    ax.add_patch(Circle(geom.centro, 0.20, facecolor=CENTROSOME_GREEN, edgecolor="none", zorder=7))
    _draw_kinetochore_symbol(ax, geom.kt, p.kt_r)

    status = "Captured" if trial.captured else "Not captured"
    mode = "Directed" if with_gradient else "Random"
    overlay = f"{mode} | {status}\nTime: {time_s:.2f} s\nMT emitted: {total_mts}"
    ax.text(0.02, 0.98, overlay, transform=ax.transAxes, ha="left", va="top", fontsize=11, color="#202020",
            bbox=dict(boxstyle="round,pad=0.35", facecolor=(1, 1, 1, 0.60), edgecolor="none"), zorder=10)
    ax.text((geom.centro[0] + geom.kt[0]) / 2, 0.22, f"{p.kt_dist:.1f} μm", ha="center", va="bottom", fontsize=10, color="#666666", zorder=9)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-8, 8)
    ax.set_ylim(-4.5, 4.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout(pad=0.4)
    if save_path is not None:
        save_figure_robust(fig, Path(save_path))
        if created_fig:
            plt.close(fig)
    return fig, ax


def plot_scene_from_trial(trial: TrialResult, p: SimParams, with_gradient: bool, save_path: Optional[Path] = None, ax=None):
    return draw_scene_frame(trial, p, with_gradient, frame_index=None, save_path=save_path, ax=ax)


def run_single_mt_trial(with_gradient: bool, p: SimParams, trial_index: int, rng: random.Random, record_history: bool = True) -> TrialResult:
    """Movie/demo mode: emit exactly ONE MT immediately at Time = 0 s.

    This special rule is ONLY for the Movie S12 demonstration/preview. It removes
    the uninformative waiting period before the first MT appears. The quantitative
    batch model is unchanged and still uses mt_rate to determine MT emission timing.
    """
    geom = build_geometry(p)
    max_steps = max(1, int(round(p.max_time / p.dt)))
    condition = "Gradient" if with_gradient else "Random"

    # Movie S12 starts with one MT already emitted at t=0. The MT initially has
    # zero visible length and begins growing on the first simulation step.
    mts: List[Microtubule] = [make_mt(1, 0, p, geom, rng)]
    step = 0
    captured = False
    capture_step = -1
    n_mt_at_capture = 0
    first_captured_mt_id: Optional[int] = None
    history: Optional[List[Dict[str, object]]] = [] if record_history else None

    if record_history:
        history.append(snapshot_trial_state(mts, 0, p))

    while step < max_steps:
        step += 1

        for mt in mts:
            if step_mt(mt, with_gradient, p, geom, rng, step):
                if not captured:
                    captured = True
                    capture_step = step
                    n_mt_at_capture = 1
                    first_captured_mt_id = mt.mt_id

        if record_history:
            history.append(snapshot_trial_state(mts, step, p))

        still_active = any((not mt.captured) and mt.state != "stuck" for mt in mts)
        if not still_active:
            break

    return TrialResult(
        condition=condition,
        trial=trial_index,
        captured=captured,
        capture_time_s=capture_step * p.dt if capture_step > 0 else p.max_time,
        capture_step=capture_step,
        n_mt_at_capture=n_mt_at_capture,
        first_captured_mt_id=first_captured_mt_id,
        stop_step=step,
        stop_time_s=step * p.dt,
        total_mts_emitted=1,
        mts=mts,
        history=history,
    )


def build_movie_s12_demo(base_params: SimParams, fps: float = 10.0) -> DemoMovie:
    left_params = SimParams(**asdict(base_params))
    left_params.mech_bias = False
    left_params.mech_inst = False
    right_params = SimParams(**asdict(base_params))
    right_params.mech_bias = True
    right_params.mech_inst = True

    # Movie S12 uses the SAME seed on both sides so the two panels begin from
    # exactly the same MT emission direction. If no seed was supplied, generate
    # a fresh random seed instead of falling back to one fixed demonstration.
    demo_seed = base_params.seed if base_params.seed is not None else random.SystemRandom().randint(1, 2147483647)
    left_rng = random.Random(demo_seed)
    right_rng = random.Random(demo_seed)

    left_trial = run_single_mt_trial(False, left_params, 1, left_rng, record_history=True)
    right_trial = run_single_mt_trial(True, right_params, 1, right_rng, record_history=True)
    return DemoMovie(
        title=MOVIE_S12_TITLE,
        left_title="Without gradient",
        right_title="With gradient",
        left_trial=left_trial,
        right_trial=right_trial,
        params=base_params,
        demo_seed=demo_seed,
        fps=float(fps),
    )


def _draw_single_demo_panel(ax, trial: TrialResult, p: SimParams, with_gradient: bool, frame_index: int):
    from matplotlib.patches import Circle
    geom = build_geometry(p)
    ax.set_facecolor(SCENE_BG)
    if with_gradient:
        _draw_soft_gradient(ax, geom, p)

    snap_list = trial.history or []
    if snap_list:
        idx = max(0, min(frame_index, len(snap_list) - 1))
        snapshot = snap_list[idx]
        mts_to_draw = snapshot.get("mts", [])
    else:
        mts_to_draw = []

    for mt in mts_to_draw:
        path = mt.get("path", [])
        if len(path) < 2:
            continue
        x = [pt[0] for pt in path]
        y = [pt[1] for pt in path]
        # Movie S12 uses one visual convention only: every MT is a solid green line.
        ax.plot(x, y, color=MT_GREEN, linewidth=2.4, linestyle='-', alpha=1.0,
                solid_capstyle='round', zorder=3)

    ax.add_patch(Circle(geom.centro, 0.20, facecolor=CENTROSOME_GREEN, edgecolor='none', zorder=7))
    _draw_kinetochore_symbol(ax, geom.kt, p.kt_r)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(-5.2, 6.5)
    ax.set_ylim(-3.0, 3.0)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.6)
        spine.set_edgecolor('#202020')


def draw_movie_s12_frame(movie: DemoMovie, frame_index: int, save_path: Optional[Path] = None, fig=None):
    """Draw one Movie S12 frame.

    Important: when no figure is supplied, this function creates an OFF-SCREEN
    Agg figure rather than a TkAgg/pyplot figure. That keeps high-resolution
    video export independent of Tk and avoids Tk_GetPixmap memory failures.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.patches import Rectangle, Circle
    from matplotlib.offsetbox import DrawingArea, AnnotationBbox
    from matplotlib.lines import Line2D

    created_fig = False
    if fig is None:
        export_dpi = 300
        fig = Figure(figsize=(movie.width_px / export_dpi, movie.height_px / export_dpi), dpi=export_dpi)
        FigureCanvasAgg(fig)
        created_fig = True
    fig.patch.set_facecolor('#000000')
    fig.clf()

    # Keep exactly the same point sizes in preview and export. Because both
    # use the same 2:1 composition, down-scaling the export preserves layout.
    title_fs = 13.2       # largest text on the page
    panel_fs = 12.0
    time_fs = 11.5
    legend_fs = 11.0

    # Main title remains outside the large content frame.
    fig.text(0.5, 0.960, movie.title, ha='center', va='top', fontsize=title_fs,
             color='white', family='Arial')

    # Large white frame enclosing everything except the main title.
    frame_left, frame_bottom, frame_right, frame_top = 0.035, 0.065, 0.965, 0.875
    fig.add_artist(Rectangle(
        (frame_left, frame_bottom), frame_right - frame_left, frame_top - frame_bottom,
        transform=fig.transFigure, fill=False, edgecolor='white', linewidth=1.1, zorder=50
    ))

    # Much tighter inner layout: panels are pushed upward and labels sit just above them.
    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        height_ratios=[1.0, 0.19],
        width_ratios=[1.0, 1.0],
        left=0.065,
        right=0.935,
        top=0.785,
        bottom=0.095,
        hspace=0.075,
        wspace=0.075,
    )

    left_history = movie.left_trial.history or []
    right_history = movie.right_trial.history or []
    left_idx = max(0, min(frame_index, len(left_history) - 1)) if left_history else 0
    right_idx = max(0, min(frame_index, len(right_history) - 1)) if right_history else 0
    time_s = 0.0
    if right_history:
        time_s = float(right_history[right_idx].get('time_s', 0.0))
    elif left_history:
        time_s = float(left_history[left_idx].get('time_s', 0.0))

    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1])
    _draw_single_demo_panel(ax_left, movie.left_trial, movie.params, False, left_idx)
    _draw_single_demo_panel(ax_right, movie.right_trial, movie.params, True, right_idx)

    # Condition labels directly above each animation frame with only a small gap.
    left_box = ax_left.get_position()
    right_box = ax_right.get_position()
    label_y = max(left_box.y1, right_box.y1) + 0.014
    fig.text((left_box.x0 + left_box.x1) / 2.0, label_y, movie.left_title,
             ha='center', va='bottom', fontsize=panel_fs, color='white', family='Arial')
    fig.text((right_box.x0 + right_box.x1) / 2.0, label_y, movie.right_title,
             ha='center', va='bottom', fontsize=panel_fs, color='white', family='Arial')

    # Time moved to the upper-right corner of the large frame.
    fig.text(frame_right - 0.025, frame_top - 0.035, f"Time: {time_s:.0f} s",
             ha='right', va='top', fontsize=time_fs, color='white', family='Arial')

    legend_ax = fig.add_subplot(gs[1, :])
    legend_ax.set_facecolor('#000000')
    legend_ax.axis('off')

    # Four-item legend only.
    y = 0.60
    items = [0.075, 0.285, 0.505, 0.715]

    x = items[0]
    legend_ax.scatter([x], [y], s=190, color=CENTROSOME_GREEN,
                      transform=legend_ax.transAxes, clip_on=False)
    legend_ax.text(x + 0.025, y, 'Centrosome', transform=legend_ax.transAxes,
                   va='center', fontsize=legend_fs, color='white', family='Arial')

    x = items[1]
    # Draw the kinetochore legend using a pixel-based mini canvas so the two circles
    # stay perfectly round and vertically tangent, matching the panel symbol.
    da = DrawingArea(30, 40, 0, 0)
    kin_r_px = 7.0
    kin_cx_px = 15.0
    kin_mid_y_px = 20.0
    # Exactly tangent: center-to-center distance = 2 * radius.
    da.add_artist(Circle((kin_cx_px, kin_mid_y_px + kin_r_px), kin_r_px,
                         facecolor=KT_MAGENTA, edgecolor='#8A007D', linewidth=0.9))
    da.add_artist(Circle((kin_cx_px, kin_mid_y_px - kin_r_px), kin_r_px,
                         facecolor=KT_MAGENTA, edgecolor='#8A007D', linewidth=0.9))
    ab = AnnotationBbox(da, (x + 0.004, y), xycoords=legend_ax.transAxes, frameon=False, box_alignment=(0.5, 0.5))
    legend_ax.add_artist(ab)
    legend_ax.text(x + 0.032, y, 'Kinetochore', transform=legend_ax.transAxes,
                   va='center', fontsize=legend_fs, color='white', family='Arial')

    x = items[2]
    legend_ax.plot([x, x + 0.055], [y, y], color=MT_GREEN, linewidth=2.4,
                   linestyle='-', transform=legend_ax.transAxes, solid_capstyle='round')
    legend_ax.text(x + 0.070, y, 'Microtubule', transform=legend_ax.transAxes,
                   va='center', fontsize=legend_fs, color='white', family='Arial')

    x = items[3]
    # Smooth circular radial gradient icon. Scatter markers stay perfectly circular
    # regardless of the legend axes aspect ratio, and many concentric layers mimic
    # the same soft gradient shown in the right animation panel.
    gx = x + 0.022
    max_size = 1050.0
    n_layers = 44
    for i in range(n_layers, 0, -1):
        frac = i / n_layers  # 1 = outer rim, ~0 = center
        marker_size = max_size * (frac ** 2)
        rgba = _gradient_color_at(frac)
        legend_ax.scatter([gx], [y], s=marker_size, color=[rgba],
                          transform=legend_ax.transAxes, clip_on=False,
                          edgecolors='none', zorder=2)
    legend_ax.text(x + 0.060, y, 'Hypothetical guidance\ngradient',
                   transform=legend_ax.transAxes, va='center', fontsize=legend_fs,
                   color='white', family='Arial', linespacing=0.95)

    # Scale bar stays at bottom-right, separated from the legend text.
    sb_x2 = 0.930
    sb_x1 = sb_x2 - 0.075
    sb_y = 0.086
    fig.lines.append(Line2D([sb_x1, sb_x2], [sb_y, sb_y], transform=fig.transFigure,
                            color='white', linewidth=2.8, solid_capstyle='butt', zorder=60))
    fig.text(sb_x2, sb_y + 0.012, f'{movie.scale_bar_um:.0f} μm', ha='right', va='bottom',
             fontsize=legend_fs, color='white', family='Arial')

    if save_path is not None:
        # Do not use bbox_inches="tight" here: Movie S12 has a fixed 2:1 canvas.
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        ext = save_path.suffix.lower().lstrip('.')
        kwargs = {"facecolor": fig.get_facecolor()}
        if ext in {"png", "jpg", "jpeg", "tif", "tiff"}:
            kwargs["dpi"] = fig.dpi
        fig.savefig(str(save_path), **kwargs)
    return fig

def _figure_to_rgb_array(fig):
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8).reshape(h, w, 4)
    return buf[..., :3].copy()


def export_movie_s12_mp4(movie: DemoMovie, save_path: Path):
    """Export exact 3000×1500 MP4 using an off-screen Agg canvas.

    A single high-resolution Figure is reused for every frame. No Tk pixmaps are
    created during export, which greatly reduces graphics-memory pressure.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    try:
        import imageio.v2 as imageio
    except Exception as exc:
        raise RuntimeError('imageio is required for MP4 export. Please install imageio and imageio-ffmpeg.') from exc

    save_path = Path(save_path)
    save_path.parent.mkdir(parents=True, exist_ok=True)

    dpi = 300
    fig = Figure(figsize=(movie.width_px / dpi, movie.height_px / dpi), dpi=dpi)
    canvas = FigureCanvasAgg(fig)

    # macro_block_size=1 prevents ImageIO from silently resizing 3000×1500
    # to 3008×1504. Both dimensions are even, so yuv420p/H.264 is valid.
    writer = imageio.get_writer(
        str(save_path),
        fps=movie.fps,
        codec='libx264',
        quality=8,
        format='FFMPEG',
        pixelformat='yuv420p',
        macro_block_size=1,
    )
    try:
        for i in range(movie.n_frames):
            draw_movie_s12_frame(movie, i, fig=fig)
            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            frame = np.ascontiguousarray(rgba[..., :3])
            writer.append_data(frame)
    finally:
        writer.close()
        fig.clear()


# -----------------------------
# GUI
# -----------------------------

class MTSimGUI:
    def __init__(self):
        ensure_matplotlib_backend(gui=True)
        import tkinter as tk
        from tkinter import ttk
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        self.tk = tk
        self.ttk = ttk
        self.root = tk.Tk()
        self.root.title("MT Directed Capture Simulation v10 - Python (enhanced GUI)")
        self.root.geometry("1580x900")
        self.root.minsize(1420, 820)

        self.params = SimParams()
        self.batch: Optional[BatchResult] = None
        self.batch_signature: Optional[str] = None
        self.batch_is_stale: bool = False
        self.last_trial: Optional[TrialResult] = None
        self.last_trial_with_gradient: bool = True
        self.movie_demo: Optional[DemoMovie] = None
        self.view_mode: str = "scene"
        self._movie_preview_fig = None
        self._movie_preview_canvas = None
        self.output_dir = tk.StringVar(value=str(Path.home() / "Desktop"))
        self.status = tk.StringVar(value="Ready. Legacy theta endpoint logic restored: first capture stops emission; already-emitted MTs continue for theta statistics.")
        self.plot_mode = tk.StringVar(value="Capture-time")
        self.frame_var = tk.IntVar(value=0)
        self.play_fps = tk.DoubleVar(value=8.0)
        self.frame_info = tk.StringVar(value="Frame: 0 / 0 | Time: 0.00 s")
        self.vars: Dict[str, object] = {}
        self.playing = False
        self.play_job = None
        self.batch_running = False
        self.export_running = False
        self.export_queue = queue.Queue()
        self.export_long_trace = tk.BooleanVar(value=False)
        self._slider_internal_update = False

        self._build_gui(Figure, FigureCanvasTkAgg, NavigationToolbar2Tk)
        self._draw_static_scene()

    def _build_gui(self, Figure, FigureCanvasTkAgg, NavigationToolbar2Tk):
        ttk = self.ttk
        tk = self.tk

        main = ttk.Frame(self.root, padding=8)
        main.pack(fill="both", expand=True)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        self.left = ttk.Frame(main)
        self.left.grid(row=0, column=0, sticky="nsw", padx=(0, 8))

        right = ttk.Frame(main)
        right.grid(row=0, column=1, sticky="nsew")
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        self._build_controls()

        topbar = ttk.Frame(right)
        topbar.grid(row=0, column=0, sticky="ew", pady=(0, 6))
        topbar.columnconfigure(0, weight=1)
        ttk.Label(topbar, text="Animation / plot view", font=("Arial", 12, "bold")).grid(row=0, column=0, sticky="w")
        ttk.Label(topbar, textvariable=self.frame_info, font=("Arial", 10)).grid(row=0, column=1, sticky="e")

        fig = Figure(figsize=(10.0, 5.0), dpi=100)
        fig.patch.set_facecolor(SCENE_BG)
        self.ax = fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(fig, master=right)
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.grid(row=1, column=0, sticky="nsew")

        toolbar_frame = ttk.Frame(right)
        toolbar_frame.grid(row=2, column=0, sticky="ew")
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side="left", fill="x")

        status_bar = ttk.Label(self.root, textvariable=self.status, relief="sunken", anchor="w")
        status_bar.pack(side="bottom", fill="x")

    def _add_slider(self, parent, row, key, label, minv, maxv, resolution, fmt):
        """Add one adjustable parameter with BOTH a slider and a numeric input box.

        The slider is kept for quick exploration.
        The numeric entry box is the precise control: type a value and press Enter,
        or leave the box / run a simulation to apply the typed value.
        """
        tk = self.tk
        var = tk.DoubleVar(value=getattr(self.params, key))
        entry_var = tk.StringVar(value=fmt(var.get()).replace("%", ""))
        self.vars[key] = var
        self.vars[f"{key}_entry_var"] = entry_var
        self.vars[f"{key}_limits"] = (float(minv), float(maxv), float(resolution))

        tk.Label(parent, text=label, anchor="w").grid(row=row, column=0, sticky="ew", pady=2)

        scale = tk.Scale(
            parent,
            from_=minv,
            to=maxv,
            resolution=resolution,
            orient="horizontal",
            variable=var,
            length=165,
            showvalue=False,
            command=lambda _=None, k=key: self._on_slider_changed(k),
        )
        scale.grid(row=row, column=1, sticky="ew", pady=2)

        entry = tk.Entry(parent, textvariable=entry_var, width=9, justify="right")
        entry.grid(row=row, column=2, sticky="e", padx=(4, 0), pady=2)
        entry.bind("<Return>", lambda _event, k=key: self._apply_numeric_entry(k))
        entry.bind("<FocusOut>", lambda _event, k=key: self._apply_numeric_entry(k))

        val = tk.Label(parent, text=fmt(var.get()), width=8, anchor="e")
        val.grid(row=row, column=3, sticky="e", padx=(4, 0), pady=2)
        self.vars[f"{key}_label"] = (val, fmt)

    def _build_controls(self):
        ttk = self.ttk
        tk = self.tk

        param_box = ttk.LabelFrame(self.left, text="Parameters", padding=8)
        param_box.grid(row=0, column=0, sticky="new")
        r = 0
        self._add_slider(param_box, r, "mt_rate", "MTs emitted per second", 0.01, 10.0, 0.01, lambda x: f"{x:.2f}"); r += 1
        self._add_slider(param_box, r, "angle_spread_pi", "Angle spread (×π)", 0.01, 0.50, 0.01, lambda x: f"{x:.2f}"); r += 1
        self._add_slider(param_box, r, "kt_diam", "KT diameter (μm)", 0.10, 3.00, 0.01, lambda x: f"{x:.2f}"); r += 1
        self._add_slider(param_box, r, "grad_r", "Gradient radius (μm)", 0.10, 10.00, 0.10, lambda x: f"{x:.1f}"); r += 1
        self._add_slider(param_box, r, "kt_dist", "KT-centrosome distance (μm)", 0.50, 20.00, 0.10, lambda x: f"{x:.1f}"); r += 1
        self._add_slider(param_box, r, "cycle_ms", "DI cycle interval (ms)", 10, 10000, 10, lambda x: f"{x:.0f}"); r += 1
        self._add_slider(param_box, r, "grow_dist", "Grow per cycle (μm)", 0.001, 0.500, 0.001, lambda x: f"{x:.3f}"); r += 1
        self._add_slider(param_box, r, "shrink_dist", "Shrink per cycle (μm)", 0.000, 0.500, 0.001, lambda x: f"{x:.3f}"); r += 1
        self._add_slider(param_box, r, "base_cat_percent", "Base catastrophe (%)", 0.0, 100.0, 0.1, lambda x: f"{x:.1f}%"); r += 1
        self._add_slider(param_box, r, "cont_bias", "Continuous bias (0-1)", 0.0, 1.0, 0.01, lambda x: f"{x:.2f}"); r += 1
        self._add_slider(param_box, r, "max_len", "Max MT length (μm)", 1.0, 30.0, 0.1, lambda x: f"{x:.1f}"); r += 1
        self._add_slider(param_box, r, "max_time", "Max time (s)", 20, 10000, 10, lambda x: f"{x:.0f}"); r += 1

        self.vars["mech_bias"] = tk.BooleanVar(value=self.params.mech_bias)
        self.vars["mech_inst"] = tk.BooleanVar(value=self.params.mech_inst)
        ttk.Checkbutton(param_box, text="Directional bias", variable=self.vars["mech_bias"], command=lambda: self._mark_batch_dirty("mech_bias")).grid(row=r, column=0, columnspan=4, sticky="w"); r += 1
        ttk.Checkbutton(param_box, text="Dynamic instability", variable=self.vars["mech_inst"], command=lambda: self._mark_batch_dirty("mech_inst")).grid(row=r, column=0, columnspan=4, sticky="w"); r += 1
        ttk.Button(param_box, text="Load confirmed defaults", command=self._set_confirmed_default_preset).grid(row=r, column=0, columnspan=4, sticky="ew", pady=(6, 2)); r += 1
        ttk.Label(param_box, text="Trials / condition").grid(row=r, column=0, sticky="w", pady=(8, 2))
        self.vars["n_trials"] = tk.IntVar(value=self.params.n_trials)
        ttk.Entry(param_box, textvariable=self.vars["n_trials"], width=10).grid(row=r, column=1, sticky="w", pady=(8, 2)); r += 1
        ttk.Label(param_box, text="Seed (blank=random)").grid(row=r, column=0, sticky="w")
        self.vars["seed_text"] = tk.StringVar(value="")
        ttk.Entry(param_box, textvariable=self.vars["seed_text"], width=12).grid(row=r, column=1, sticky="w")
        ttk.Button(param_box, text="Randomize", command=self.randomize_seed_gui).grid(row=r, column=2, columnspan=2, sticky="ew", padx=(4, 0)); r += 1

        run_box = ttk.LabelFrame(self.left, text="Run / plots", padding=8)
        run_box.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        ttk.Button(run_box, text="Single: Directed", command=lambda: self.run_single(True)).grid(row=0, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="Single: Random", command=lambda: self.run_single(False)).grid(row=0, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="Run batch", command=self.run_batch_gui).grid(row=1, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="Export batch", command=self.export_gui).grid(row=1, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="Show histogram", command=lambda: self.show_plot("Capture-time")).grid(row=2, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="Show θ plot", command=lambda: self.show_plot("Theta")).grid(row=2, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="Reset static scene", command=self._draw_static_scene).grid(row=3, column=0, columnspan=2, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="Prepare Movie S12 demo", command=self.prepare_movie_s12_demo).grid(row=4, column=0, columnspan=2, sticky="ew", padx=2, pady=2)
        ttk.Button(run_box, text="New random Movie S12 example", command=self.randomize_movie_s12_demo).grid(row=5, column=0, columnspan=2, sticky="ew", padx=2, pady=2)

        anim_box = ttk.LabelFrame(self.left, text="Animation controls", padding=8)
        anim_box.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        ttk.Button(anim_box, text="Play", command=self.play_animation).grid(row=0, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(anim_box, text="Pause", command=self.pause_animation).grid(row=0, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(anim_box, text="◀ Prev", command=self.prev_frame).grid(row=1, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(anim_box, text="Next ▶", command=self.next_frame).grid(row=1, column=1, sticky="ew", padx=2, pady=2)
        ttk.Label(anim_box, text="Animation FPS").grid(row=2, column=0, sticky="w", pady=(6, 0))
        ttk.Entry(anim_box, textvariable=self.play_fps, width=8).grid(row=2, column=1, sticky="w", pady=(6, 0))
        ttk.Label(anim_box, text="Frame").grid(row=3, column=0, sticky="w", pady=(6, 0))
        self.frame_scale = tk.Scale(anim_box, from_=0, to=0, orient="horizontal", variable=self.frame_var, resolution=1, length=240, command=self.on_frame_slider, showvalue=False)
        self.frame_scale.grid(row=4, column=0, columnspan=2, sticky="ew")

        export_box = ttk.LabelFrame(self.left, text="Current paused frame export", padding=8)
        export_box.grid(row=3, column=0, sticky="ew", pady=(8, 0))
        ttk.Button(export_box, text="Export paused SVG", command=lambda: self.export_current_frame("svg")).grid(row=0, column=0, sticky="ew", padx=2, pady=2)
        ttk.Button(export_box, text="Export paused PDF", command=lambda: self.export_current_frame("pdf")).grid(row=0, column=1, sticky="ew", padx=2, pady=2)
        ttk.Button(export_box, text="Export paused PNG", command=lambda: self.export_current_frame("png")).grid(row=1, column=0, columnspan=2, sticky="ew", padx=2, pady=2)
        ttk.Button(export_box, text="Export Movie S12 MP4", command=self.export_movie_s12_mp4_gui).grid(row=2, column=0, columnspan=2, sticky="ew", padx=2, pady=2)

    def _current_params_signature(self) -> str:
        return params_signature(self.get_params_from_gui())

    def _mark_batch_dirty(self, key: str = ""):
        """Mark the displayed batch plot as stale after any parameter change."""
        if getattr(self, "batch", None) is not None and not getattr(self, "batch_running", False):
            self.batch_is_stale = True
            suffix = f" ({key})" if key else ""
            self.status.set(
                "Parameters changed" + suffix + ". Existing batch plot is now OLD. Click Run batch again before interpreting plots."
            )

    def _ensure_batch_matches_current_params(self, action: str = "show this plot") -> bool:
        """Return True only if last batch was generated from current GUI params."""
        from tkinter import messagebox
        if not self.batch:
            self.status.set("Run batch first, then " + action + ".")
            return False
        current_sig = self._current_params_signature()
        if self.batch_signature is None:
            self.batch_signature = params_signature(self.batch.params)
        if current_sig != self.batch_signature or self.batch_is_stale:
            messagebox.showwarning(
                "Batch is outdated",
                "Current parameters differ from the last Run batch result.\n\n"
                "The plot on the right may be an OLD result.\n\n"
                "Please click Run batch again before showing/exporting plots."
            )
            self.status.set("Batch result is outdated. Click Run batch again.")
            return False
        return True

    def _set_confirmed_default_preset(self):
        """Load the user-confirmed default parameter set."""
        preset = {
            "mt_rate": 0.02,
            "angle_spread_pi": 0.30,
            "kt_diam": 0.80,
            "grad_r": 2.5,
            "kt_dist": 6.7,
            "cycle_ms": 500.0,
            "grow_dist": 0.110,
            "shrink_dist": 0.100,
            "base_cat_percent": 20.0,
            "cont_bias": 0.20,
            "max_len": 12.0,
            "max_time": 3000.0,
        }
        for key, value in preset.items():
            if key in self.vars:
                self.vars[key].set(value)
                entry_var = self.vars.get(f"{key}_entry_var")
                label_tuple = self.vars.get(f"{key}_label")
                if entry_var is not None and label_tuple is not None:
                    _label, fmt = label_tuple
                    entry_var.set(fmt(value).replace("%", ""))
        if "mech_bias" in self.vars:
            self.vars["mech_bias"].set(True)
        if "mech_inst" in self.vars:
            self.vars["mech_inst"].set(True)
        if "seed_text" in self.vars:
            self.vars["seed_text"].set("")
        self._update_value_labels()
        self._mark_batch_dirty("confirmed default preset")
        self.status.set("Confirmed default preset loaded. Seed left blank/random unless you type one. Click Run batch to generate a fresh θ plot.")

    def _on_slider_changed(self, key: str):
        """Synchronize numeric entry when the slider is moved."""
        label, fmt = self.vars[f"{key}_label"]
        var = self.vars[key]
        entry_var = self.vars.get(f"{key}_entry_var")
        value = float(var.get())
        label.configure(text=fmt(value))
        if entry_var is not None:
            entry_var.set(fmt(value).replace("%", ""))
        self._mark_batch_dirty(key)

    def _apply_numeric_entry(self, key: str):
        """Apply a typed value from the numeric box to the slider variable.

        Values outside the slider range are clamped to the allowed range.
        A trailing percent sign is allowed for Base catastrophe.
        """
        label, fmt = self.vars[f"{key}_label"]
        var = self.vars[key]
        entry_var = self.vars.get(f"{key}_entry_var")
        minv, maxv, _resolution = self.vars.get(f"{key}_limits", (-math.inf, math.inf, 0.0))

        if entry_var is None:
            return

        raw = str(entry_var.get()).strip().replace("%", "")
        try:
            value = float(raw)
        except Exception:
            value = float(var.get())

        if value < minv:
            value = minv
        if value > maxv:
            value = maxv

        old_value = float(var.get())
        var.set(value)
        label.configure(text=fmt(value))
        entry_var.set(fmt(value).replace("%", ""))
        if abs(old_value - value) > 1e-12:
            self._mark_batch_dirty(key)

    def _update_value_labels(self):
        for key, var in list(self.vars.items()):
            label_key = f"{key}_label"
            if label_key in self.vars:
                label, fmt = self.vars[label_key]
                try:
                    value = float(var.get())
                    label.configure(text=fmt(value))
                    entry_var = self.vars.get(f"{key}_entry_var")
                    if entry_var is not None:
                        entry_var.set(fmt(value).replace("%", ""))
                except Exception:
                    pass

    def get_params_from_gui(self) -> SimParams:
        # Apply all numeric text boxes before reading parameters.
        # This means the user can type a value and directly click Run batch without pressing Enter first.
        for _key in [
            "mt_rate", "angle_spread_pi", "kt_diam", "grad_r", "kt_dist", "cycle_ms",
            "grow_dist", "shrink_dist", "base_cat_percent", "cont_bias", "max_len", "max_time"
        ]:
            if _key in self.vars:
                self._apply_numeric_entry(_key)

        seed_text = self.vars["seed_text"].get().strip()
        seed = int(seed_text) if seed_text else None
        return SimParams(
            mt_rate=float(self.vars["mt_rate"].get()),
            angle_spread_pi=float(self.vars["angle_spread_pi"].get()),
            kt_diam=float(self.vars["kt_diam"].get()),
            grad_r=float(self.vars["grad_r"].get()),
            kt_dist=float(self.vars["kt_dist"].get()),
            cycle_ms=float(self.vars["cycle_ms"].get()),
            grow_dist=float(self.vars["grow_dist"].get()),
            shrink_dist=float(self.vars["shrink_dist"].get()),
            base_cat_percent=float(self.vars["base_cat_percent"].get()),
            cont_bias=float(self.vars["cont_bias"].get()),
            max_len=float(self.vars["max_len"].get()),
            max_time=float(self.vars["max_time"].get()),
            mech_bias=bool(self.vars["mech_bias"].get()),
            mech_inst=bool(self.vars["mech_inst"].get()),
            n_trials=int(self.vars["n_trials"].get()),
            seed=seed,
        )

    def _draw_static_scene(self):
        self.pause_animation()
        self.movie_demo = None
        self.view_mode = "scene"
        p = self.get_params_from_gui()
        geom = build_geometry(p)
        self.ax.clear()
        self.ax.set_facecolor(SCENE_BG)
        self.ax.figure.patch.set_facecolor(SCENE_BG)
        _draw_soft_gradient(self.ax, geom, p)
        self.ax.plot([geom.centro[0], geom.kt[0]], [geom.centro[1], geom.kt[1]], linestyle="--", color="#C3C3C3", linewidth=1.0, zorder=1)
        from matplotlib.patches import Circle
        self.ax.add_patch(Circle(geom.centro, 0.20, facecolor=CENTROSOME_GREEN, edgecolor="none", zorder=7))
        _draw_kinetochore_symbol(self.ax, geom.kt, p.kt_r)
        self.ax.text((geom.centro[0] + geom.kt[0]) / 2, 0.22, f"{p.kt_dist:.1f} μm", ha="center", va="bottom", fontsize=10, color="#666666", zorder=9)
        self.ax.text(0.02, 0.98, "Static geometry preview", transform=self.ax.transAxes, ha="left", va="top", fontsize=11, color="#202020", bbox=dict(boxstyle="round,pad=0.35", facecolor=(1,1,1,0.60), edgecolor="none"), zorder=10)
        self.ax.set_aspect("equal", adjustable="box")
        self.ax.set_xlim(-8, 8)
        self.ax.set_ylim(-4.5, 4.5)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        for spine in self.ax.spines.values():
            spine.set_visible(False)
        self.canvas.draw_idle()
        self.frame_info.set("Frame: 0 / 0 | Time: 0.00 s")
        self.status.set("Static preview reset. Updated color style and high-resolution canvas are active.")

    def _current_frame_count(self) -> int:
        if self.movie_demo is not None:
            return self.movie_demo.n_frames
        if self.last_trial and self.last_trial.history:
            return len(self.last_trial.history)
        return 0

    def _set_animation_range(self):
        n = max(0, self._current_frame_count() - 1)
        self.frame_scale.configure(to=n)

    def _set_frame_value_safely(self, value: int):
        """Set the frame slider without triggering on_frame_slider -> pause."""
        self._slider_internal_update = True
        try:
            self.frame_var.set(value)
        finally:
            self._slider_internal_update = False

    def _render_movie_preview_rgb(self, frame_index: int):
        """Render the same 2:1 Movie S12 composition off-screen for GUI preview."""
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        if self.movie_demo is None:
            return None
        if self._movie_preview_fig is None or self._movie_preview_canvas is None:
            # 1200×600 is light enough for interactive preview but keeps exact 2:1 layout.
            self._movie_preview_fig = Figure(figsize=(12.0, 6.0), dpi=100)
            self._movie_preview_canvas = FigureCanvasAgg(self._movie_preview_fig)

        draw_movie_s12_frame(self.movie_demo, frame_index, fig=self._movie_preview_fig)
        self._movie_preview_canvas.draw()
        rgba = np.asarray(self._movie_preview_canvas.buffer_rgba())
        return np.ascontiguousarray(rgba[..., :3])

    def draw_current_frame(self):
        if self.movie_demo is not None:
            idx = max(0, min(int(self.frame_var.get()), self.movie_demo.n_frames - 1))
            preview_rgb = self._render_movie_preview_rgb(idx)
            self.ax.clear()
            self.ax.set_facecolor('#000000')
            self.ax.figure.patch.set_facecolor('#000000')
            if preview_rgb is not None:
                self.ax.imshow(preview_rgb, interpolation='bilinear')
            self.ax.set_aspect('equal', adjustable='box')
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            for spine in self.ax.spines.values():
                spine.set_visible(False)
            self.canvas.draw_idle()
            time_s = 0.0
            right_hist = self.movie_demo.right_trial.history or []
            left_hist = self.movie_demo.left_trial.history or []
            if right_hist:
                time_s = float(right_hist[min(idx, len(right_hist) - 1)].get('time_s', 0.0))
            elif left_hist:
                time_s = float(left_hist[min(idx, len(left_hist) - 1)].get('time_s', 0.0))
            self.frame_info.set(f"Frame: {idx} / {self.movie_demo.n_frames - 1} | Time: {time_s:.2f} s")
            return

        if not self.last_trial:
            self._draw_static_scene()
            return
        idx = max(0, min(int(self.frame_var.get()), len(self.last_trial.history or [0]) - 1 if self.last_trial.history else 0))
        draw_scene_frame(self.last_trial, self.get_params_from_gui(), self.last_trial_with_gradient, frame_index=idx, ax=self.ax)
        self.canvas.draw_idle()
        if self.last_trial.history:
            snap = self.last_trial.history[idx]
            self.frame_info.set(f"Frame: {idx} / {len(self.last_trial.history)-1} | Time: {snap['time_s']:.2f} s")

    def _movie_fps_from_gui(self) -> float:
        """Read the Animation FPS field for both preview playback and MP4 export."""
        try:
            fps = float(self.play_fps.get())
        except Exception:
            fps = 10.0
        if not math.isfinite(fps) or fps <= 0:
            fps = 10.0
            self.play_fps.set(fps)
        return fps

    def randomize_seed_gui(self):
        """Generate and display a fresh reproducible seed for any subsequent run."""
        seed = random.SystemRandom().randint(1, 2147483647)
        self.vars["seed_text"].set(str(seed))
        self._mark_batch_dirty("seed")
        self.status.set(f"Random seed generated: {seed}. Run/prepare again to use this seed.")
        return seed

    def randomize_movie_s12_demo(self):
        """Prepare a genuinely new Movie S12 example while keeping left/right matched."""
        seed = self.randomize_seed_gui()
        self.prepare_movie_s12_demo()
        self.status.set(
            f"New random Movie S12 example prepared with seed {seed}. "
            "Left and right use the same seed/initial MT direction for a fair comparison."
        )

    def prepare_movie_s12_demo(self):
        self.pause_animation()
        self.last_trial = None
        # Blank now truly means random. Generate it visibly so the exact example
        # can be reproduced later by keeping/copying the displayed seed.
        if not self.vars["seed_text"].get().strip():
            self.vars["seed_text"].set(str(random.SystemRandom().randint(1, 2147483647)))
        base_p = self.get_params_from_gui()
        movie_fps = self._movie_fps_from_gui()
        self.movie_demo = build_movie_s12_demo(base_p, fps=movie_fps)
        self._movie_preview_fig = None
        self._movie_preview_canvas = None
        self.view_mode = "movie_s12"
        self._set_animation_range()
        self._set_frame_value_safely(0)
        self.draw_current_frame()
        self.status.set(
            f"Movie S12 demo prepared with seed {self.movie_demo.demo_seed} at {self.movie_demo.fps:g} FPS: "
            "one MT is emitted immediately at Time = 0 s; "
            "left = Without gradient (bias OFF, DI OFF); right = With gradient (bias ON, DI ON). "
            "Left/right share the same initial random seed."
        )

    def run_single(self, with_gradient: bool):
        self.pause_animation()
        self.movie_demo = None
        self.view_mode = "scene"
        p = self.get_params_from_gui()
        rng = random.Random(p.seed) if p.seed is not None else random.Random()
        trial = run_trial(with_gradient, p, trial_index=1, rng=rng, record_history=True)
        self.last_trial = trial
        self.last_trial_with_gradient = with_gradient
        self._set_animation_range()
        self.frame_var.set(0)
        self.draw_current_frame()
        status = "captured" if trial.captured else "not captured"
        self.status.set(f"Single {'Directed' if with_gradient else 'Random'} run prepared: {status}, first contact={trial.capture_time_s:.2f}s if captured, total MTs emitted={trial.total_mts_emitted}, theta stop={trial.stop_time_s:.2f}s. Use Play/Pause or frame slider.")

    def play_animation(self):
        if self.movie_demo is None and (not self.last_trial or not self.last_trial.history):
            self.prepare_movie_s12_demo()

        total_n = self._current_frame_count()
        if total_n <= 0:
            self.status.set("No animation frames available.")
            return

        if int(self.frame_var.get()) >= total_n - 1:
            self._set_frame_value_safely(0)

        self.playing = True
        self.status.set("Animation playing.")
        self._animate_step()

    def _animate_step(self):
        total_n = self._current_frame_count()
        if not self.playing or total_n <= 0:
            return
        idx = int(self.frame_var.get())
        self.draw_current_frame()
        if idx >= total_n - 1:
            self.playing = False
            self.play_job = None
            self.status.set("Animation finished. You can pause, scrub, or export the current frame.")
            return
        self._set_frame_value_safely(idx + 1)
        delay = max(20, int(round(1000.0 / max(0.1, float(self.play_fps.get())))))
        self.play_job = self.root.after(delay, self._animate_step)

    def pause_animation(self):
        self.playing = False
        if self.play_job is not None:
            try:
                self.root.after_cancel(self.play_job)
            except Exception:
                pass
        self.play_job = None

    def prev_frame(self):
        self.pause_animation()
        self.frame_var.set(max(0, int(self.frame_var.get()) - 1))
        self.draw_current_frame()

    def next_frame(self):
        self.pause_animation()
        mx = max(0, self._current_frame_count() - 1)
        self.frame_var.set(min(mx, int(self.frame_var.get()) + 1))
        self.draw_current_frame()

    def on_frame_slider(self, _value=None):
        if getattr(self, "_slider_internal_update", False):
            return
        self.pause_animation()
        if self.movie_demo is not None or (self.last_trial and self.last_trial.history):
            self.draw_current_frame()

    def run_batch_gui(self):
        """Run batch in a background thread so the GUI does not appear frozen."""
        self.pause_animation()
        if getattr(self, "batch_running", False):
            self.status.set("Batch is already running...")
            return

        p = self.get_params_from_gui()
        self.batch_running = True
        self.status.set("Batch running in background... You can wait without clicking again.")

        def worker():
            try:
                batch = run_batch(p, progress_cb=None)

                def finish_ok():
                    self.batch = batch
                    self.batch_signature = params_signature(batch.params)
                    self.batch_is_stale = False
                    self.batch_running = False
                    self.status.set("Batch finished with current parameters. Showing capture-time histogram. " + self._summary_text())
                    plot_capture_histogram(self.batch, ax=self.ax)
                    self.canvas.draw_idle()

                self.root.after(0, finish_ok)
            except Exception as exc:
                def finish_err():
                    self.batch_running = False
                    self.status.set(f"Batch failed: {exc}")

                self.root.after(0, finish_err)

        threading.Thread(target=worker, daemon=True).start()

    def _summary_text(self) -> str:
        if not self.batch:
            return ""
        rows = summary_rows(self.batch)
        parts = []
        for row in rows:
            ct = row["capture_time_mean_s"]
            parts.append(f"{row['condition']}: {row['captured_trials']}/{row['n_trials']} captured, mean time={ct:.1f}s" if ct is not None else f"{row['condition']}: no capture")
        return " | ".join(parts)

    def show_plot(self, mode: str):
        self.pause_animation()
        if not self._ensure_batch_matches_current_params("show plot"):
            return
        if mode == "Theta":
            plot_theta_vs_frame(self.batch, ax=self.ax)
        else:
            plot_capture_histogram(self.batch, ax=self.ax)
        self.canvas.draw_idle()
        self.status.set(f"Showing {mode} plot. Export will include the exact table used for this plot.")

    def choose_output_dir(self):
        from tkinter import filedialog
        d = filedialog.askdirectory(initialdir=self.output_dir.get())
        if d:
            self.output_dir.set(d)

    def export_current_frame(self, ext: str):
        """Export the current paused animation frame as SVG/PDF/PNG.

        Fixes in this version:
        - robust try/except with a visible error popup instead of silent failure;
        - forced correct suffix (.svg/.pdf/.png), even if the user selects "All files";
        - safe parent-folder creation;
        - vector font settings from ensure_matplotlib_backend are respected;
        - success popup confirms the exact saved file path.
        """
        self.pause_animation()

        from tkinter import filedialog, messagebox

        ext = ext.lower().lstrip(".")
        if ext not in {"svg", "pdf", "png"}:
            messagebox.showerror("Unsupported export format", f"Unsupported format: {ext}")
            return

        try:
            # Apply any typed numeric values before drawing/exporting.
            p = self.get_params_from_gui()

            initial_dir = Path(self.output_dir.get()).expanduser()
            if not initial_dir.exists():
                initial_dir = Path.home()

            mode = "directed" if self.last_trial_with_gradient else "random"
            frame_suffix = "static"
            if self.last_trial and self.last_trial.history:
                frame_suffix = f"frame_{int(self.frame_var.get()):04d}"
            default_name = f"current_paused_{mode}_{frame_suffix}.{ext}"

            path = filedialog.asksaveasfilename(
                initialdir=str(initial_dir),
                initialfile=default_name,
                defaultextension=f".{ext}",
                filetypes=[(ext.upper(), f"*.{ext}"), ("All files", "*.*")],
            )
            if not path:
                self.status.set("Current frame export canceled.")
                return

            save_path = normalize_export_save_path(path, ext)
            save_path.parent.mkdir(parents=True, exist_ok=True)

            if self.movie_demo is not None:
                frame_index = int(self.frame_var.get())
                draw_movie_s12_frame(self.movie_demo, frame_index, save_path=save_path)
            elif self.last_trial:
                frame_index = int(self.frame_var.get()) if self.last_trial.history else None
                draw_scene_frame(
                    self.last_trial,
                    p,
                    self.last_trial_with_gradient,
                    frame_index=frame_index,
                    save_path=save_path,
                )
            else:
                import matplotlib.pyplot as plt
                from matplotlib.patches import Circle

                fig, ax = plt.subplots(figsize=(8.6, 5.2), dpi=160)
                geom = build_geometry(p)
                ax.set_facecolor(SCENE_BG)
                fig.patch.set_facecolor(SCENE_BG)
                _draw_soft_gradient(ax, geom, p)
                ax.add_patch(Circle(geom.centro, 0.20, facecolor=CENTROSOME_GREEN, edgecolor="none", zorder=7))
                _draw_kinetochore_symbol(ax, geom.kt, p.kt_r)
                ax.set_aspect("equal", adjustable="box")
                ax.set_xlim(-8, 8)
                ax.set_ylim(-4.5, 4.5)
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                save_figure_robust(fig, save_path)
                plt.close(fig)

            if save_path.exists() and save_path.stat().st_size > 0:
                self.output_dir.set(str(save_path.parent))
                self.status.set(f"Exported current frame: {save_path}")
                messagebox.showinfo("Export finished", f"Current frame exported successfully.\n\n{save_path}")
            else:
                raise RuntimeError("Save command finished, but the output file was not created or is empty.")

        except Exception as exc:
            self.status.set(f"Current frame export failed: {exc}")
            messagebox.showerror(
                "Current frame export failed",
                "The current frame could not be exported.\n\n"
                f"Error:\n{exc}\n\n"
                "Try saving to a simple local folder such as Desktop, and close the target SVG/PDF if it is already open in Illustrator/Acrobat.",
            )

    def export_movie_s12_mp4_gui(self):
        self.pause_animation()
        from tkinter import filedialog, messagebox

        if self.movie_demo is None:
            self.prepare_movie_s12_demo()
        if self.movie_demo is None:
            messagebox.showwarning("No Movie S12 demo", "Movie S12 demo is not ready.")
            return

        initial_dir = Path(self.output_dir.get()).expanduser()
        if not initial_dir.exists():
            initial_dir = Path.home()

        path = filedialog.asksaveasfilename(
            initialdir=str(initial_dir),
            initialfile="Movie_S12_simulation_with_without_guidance_gradient.mp4",
            defaultextension=".mp4",
            filetypes=[("MP4", "*.mp4"), ("All files", "*.*")],
        )
        if not path:
            self.status.set("Movie S12 export canceled.")
            return

        save_path = normalize_export_save_path(path, "mp4")
        try:
            # Use the CURRENT Animation FPS field for MP4 export. This keeps
            # GUI preview playback and exported video speed synchronized.
            self.movie_demo.fps = self._movie_fps_from_gui()
            export_movie_s12_mp4(self.movie_demo, save_path)
            self.output_dir.set(str(save_path.parent))
            self.status.set(f"Movie S12 exported at {self.movie_demo.fps:g} FPS: {save_path}")
            messagebox.showinfo("Export finished", f"Movie S12 MP4 exported successfully at {self.movie_demo.fps:g} FPS.\n\n{save_path}")
        except Exception as exc:
            self.status.set(f"Movie S12 export failed: {exc}")
            messagebox.showerror(
                "Movie S12 export failed",
                "The Movie S12 MP4 could not be exported.\n\n"
                f"Error:\n{exc}\n\n"
                "If needed, install imageio and imageio-ffmpeg first."
            )

    def export_gui(self):
        """Export existing batch with visible feedback and folder selection.

        User-facing rule:
        - Click Export batch -> choose an output folder -> export starts.
        - No permanent path entry is needed in the GUI.
        """
        self.pause_animation()

        from tkinter import filedialog, messagebox

        if getattr(self, "batch_running", False):
            messagebox.showinfo("Batch still running", "Run batch is still running. Please wait until it finishes, then click Export batch.")
            self.status.set("Batch is still running. Wait for Batch finished before exporting.")
            return

        if getattr(self, "export_running", False):
            messagebox.showinfo("Export running", "Export is already running in the background.")
            self.status.set("Export is already running...")
            return

        if not self.batch:
            messagebox.showwarning("No batch data", "No batch data yet. Please click Run batch first, wait until it finishes, then click Export batch.")
            self.status.set("No batch data yet. Click Run batch first, then Export batch.")
            return

        if not self._ensure_batch_matches_current_params("export batch"):
            return

        initial_dir = self.output_dir.get() if self.output_dir.get() else str(Path.home() / "Desktop")
        chosen = filedialog.askdirectory(
            initialdir=initial_dir,
            title="Choose output folder for batch export",
        )
        if not chosen:
            self.status.set("Export canceled: no output folder selected.")
            return

        out = Path(chosen).expanduser()
        self.output_dir.set(str(out))

        try:
            out.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            messagebox.showerror("Output folder error", f"Cannot create output folder:\n{out}\n\n{exc}")
            self.status.set(f"Output folder error: {exc}")
            return

        export_long_trace_obj = getattr(self, "export_long_trace", None)
        include_long_trace = bool(export_long_trace_obj.get()) if export_long_trace_obj is not None else False

        batch_to_export = self.batch
        self.export_running = True
        self.status.set(
            f"Export started in background → {out} | "
            + "Article-ready workbook only: Capture_full/Four columns, Capture_sample500/Four columns, Theta_mean_curve, QC_summary, Parameters."
        )
        self.root.update_idletasks()

        # Clear stale queue items.
        try:
            while True:
                self.export_queue.get_nowait()
        except queue.Empty:
            pass

        def worker():
            try:
                exported = export_batch(batch_to_export, out, include_long_trace=include_long_trace)
                self.export_queue.put(("ok", exported))
            except Exception as exc:
                self.export_queue.put(("error", repr(exc)))

        threading.Thread(target=worker, daemon=True).start()
        messagebox.showinfo(
            "Export started",
            "Article-ready export has started in the background.\n\n"
            "Do not click Export repeatedly. Wait for the completion popup.\n\n"
            f"Output folder:\n{out}"
        )
        self.root.after(200, self._poll_export_queue)

    def _poll_export_queue(self):
        from tkinter import messagebox

        try:
            kind, payload = self.export_queue.get_nowait()
        except queue.Empty:
            if getattr(self, "export_running", False):
                self.root.after(200, self._poll_export_queue)
            return

        self.export_running = False

        if kind == "ok":
            exported = payload
            msg = "Export finished:\n" + "\n".join(f"{k}: {v}" for k, v in exported.items())
            self.status.set(msg)
            main_file = exported.get("article_ready_xlsx") or exported.get("article_ready_csv_dir") or next(iter(exported.values()), "")
            messagebox.showinfo("Export finished", f"Export finished successfully.\n\nMain file:\n{main_file}")
        else:
            self.status.set(f"Export failed: {payload}")
            messagebox.showerror("Export failed", f"Export failed:\n{payload}")

    def run(self):
        self.root.mainloop()


# -----------------------------
# Command line
# -----------------------------

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Python MT directed capture simulation v10")
    parser.add_argument("--nogui", action="store_true", help="Run without GUI and export batch results")
    parser.add_argument("--trials", type=int, default=3000, help="Trials per condition for --nogui")
    parser.add_argument("--out", type=str, default="./mt_sim_outputs", help="Output directory for --nogui")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed")
    parser.add_argument("--quick-test", action="store_true", help="Small no-GUI test: 3 trials per condition")
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    if args.nogui or args.quick_test:
        ensure_matplotlib_backend(gui=False)
        p = SimParams(n_trials=3 if args.quick_test else args.trials, seed=args.seed)
        print("Running no-GUI batch with parameters:")
        print(json.dumps(asdict(p), indent=2, ensure_ascii=False))
        batch = run_batch(p)
        out_dir = Path(args.out)
        exported = export_batch(batch, out_dir, include_long_trace=True)
        print("Summary:")
        for row in summary_rows(batch):
            print(row)
        print("Exported files:")
        for k, v in exported.items():
            print(f"  {k}: {v}")
        return 0

    app = MTSimGUI()
    app.run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
