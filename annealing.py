"""Polar simulated annealing for centrally symmetric orthogonal partitions.

The floating-point search follows the historical anneal.cpp scheme. It is
deliberately independent of the exact verifiers: zero energy is a candidate,
never a certificate. All compiled routines below are ordinary Python + NumPy
source accelerated with Numba; no C++ extension or compiler setup is needed.
"""
from dataclasses import dataclass
import json
import math
from pathlib import Path
import time

import numpy as np
from numba import njit


@njit(cache=True)
def sweep_energy(points, k):
    """Sweep [0, pi) using the central-symmetry reduction.

    For each representative, rank |p.n| and record sign((p.n)(p.Jn)).
    The first 2k ranks determine the overlap. Events either flip a sign or
    swap adjacent ranks. Numerically ambiguous sweeps get a large penalty.
    Like the historical implementation, the initial chamber is scored again
    after the final event; this harmless extra weight does not affect E=0.
    """
    m = len(points)
    angles = np.empty(m * (m + 1))
    left = np.empty(m * (m + 1), np.int64)
    right = np.empty(m * (m + 1), np.int64)
    signs = np.empty(m, np.int64)
    index = 0
    for i in range(m):
        x, y = points[i]
        signs[i] = int(x * y > 0)
        angles[index] = (math.atan2(y, x) + math.pi / 2) % math.pi
        left[index], right[index] = i, -1
        index += 1
        angles[index] = (math.atan2(x, -y) + math.pi / 2) % math.pi
        left[index], right[index] = i, -1
        index += 1
        for j in range(i):
            for s in (-1, 1):
                dx = x + s * points[j, 0]
                dy = y + s * points[j, 1]
                angles[index] = (math.atan2(dy, dx) + math.pi / 2) % math.pi
                left[index], right[index] = i, j
                index += 1

    order = np.argsort(-np.abs(points[:, 0]))
    rank = np.empty(m, np.int64)
    q = 0
    for a in range(m):
        rank[order[a]] = a
        if a < 2 * k:
            q += signs[order[a]]
    deficit = max(0, k + 1 - q)
    energy = deficit * deficit
    events = np.argsort(angles)
    previous = -1.0
    for a in events:
        if angles[a] == previous or angles[a] == 0.0:
            return 1e12
        previous = angles[a]
        i, j = left[a], right[a]
        if j < 0:
            if rank[i] < 2 * k:
                q += 1 - 2 * signs[i]
            signs[i] = 1 - signs[i]
        else:
            ri, rj = rank[i], rank[j]
            if abs(ri - rj) != 1:
                return 1e12
            if ri < 2 * k <= rj:
                q += signs[j] - signs[i]
            if rj < 2 * k <= ri:
                q += signs[i] - signs[j]
            rank[i], rank[j] = rj, ri
        deficit = max(0, k + 1 - q)
        energy += deficit * deficit
    return float(energy)


@njit(cache=True)
def initialize_ring(m, seed):
    """A fresh irregular ring, with no input from any saved counterexample."""
    np.random.seed(seed)
    points = np.empty((m, 2))
    for i in range(m):
        angle = (i + 0.3 * np.random.normal()) * math.pi / m
        radius = 0.7 + 0.6 * np.random.random()
        points[i, 0] = radius * math.cos(angle)
        points[i, 1] = radius * math.sin(angle)
    return points


@njit(cache=True)
def anneal_chunk(points, energy, best_points, best_energy, k, start, count, seed,
                 cycle_steps=20000):
    """Run a bounded chunk so Python can report progress and enforce a budget."""
    np.random.seed(seed)
    accepted = 0
    for offset in range(count):
        step = start + offset
        phase = (step % cycle_steps) / cycle_steps
        temperature = 12.0 * 0.015 ** phase
        scale = 0.2 * 0.01 ** phase
        if np.random.random() < 0.05:
            scale = 0.3
        i = np.random.randint(len(points))
        x, y = points[i]
        angle = math.atan2(y, x) + scale * np.random.normal()
        radius = math.hypot(x, y) * math.exp(scale * np.random.normal())
        radius = min(3.0, max(0.25, radius))
        points[i, 0] = radius * math.cos(angle)
        points[i, 1] = radius * math.sin(angle)
        proposal_energy = sweep_energy(points, k)
        if (proposal_energy <= energy or
                np.random.random() < math.exp((energy - proposal_energy) / temperature)):
            energy = proposal_energy
            accepted += 1
            if energy < best_energy:
                best_energy = energy
                best_points[:] = points
        else:
            points[i, 0], points[i, 1] = x, y
        if energy == 0:
            return energy, best_energy, offset + 1, accepted
    return energy, best_energy, count, accepted


@dataclass
class SearchResult:
    representatives: np.ndarray
    energy: float
    seed: int
    restart: int
    steps: int
    elapsed_seconds: float
    history: list
    parameters: dict

    def integer_candidate(self, k=8, scale=10_000_000):
        # Round representatives once, then negate: preserve central symmetry.
        reps = np.rint(self.representatives * scale).astype(np.int64).tolist()
        return {
            "status": "UNCERTIFIED", "k": k, "grid_energy": self.energy,
            "points": reps + [[-x, -y] for x, y in reps],
            "search": {"seed": self.seed, "restart": self.restart,
                       "proposals": self.steps, "elapsed_seconds": self.elapsed_seconds,
                       "integer_scale": scale, "initialization": "fresh irregular ring",
                       "parameters": self.parameters},
        }


def warm_up():
    """Compile kernels separately from the measured search budget."""
    p = initialize_ring(48, 0)
    e = sweep_energy(p, 8)
    anneal_chunk(p, e, p.copy(), e, 8, 0, 1, 0)


def search(seed=304, m=48, k=8, steps_per_restart=20000, restarts=6,
           max_seconds=240, chunk_size=1000, verbose=True, callback=None):
    """Search from fresh rings; return best candidate even on budget exhaustion.

    Chunk seeds are generated deterministically by NumPy's Generator. Keeping
    seed, parameters, chunk size, and dependency versions fixed fixes the
    proposal sequence. A time limit can stop different machines at different
    steps. No stored candidate is read by this function.
    """
    if not (m >= 2 and 0 < k and 2 * k < m):
        raise ValueError("Require 0 < 2*k < m for this one-sided search")
    if min(steps_per_restart, restarts, chunk_size, max_seconds) <= 0:
        raise ValueError("Search budgets must be positive")
    rng = np.random.default_rng(seed)
    started = time.monotonic()
    best_energy = float("inf")
    best_points = np.empty((m, 2))
    best_restart = 0
    history = []
    total = 0
    parameters = {"m": m, "k": k, "steps_per_restart": steps_per_restart,
                  "restarts": restarts, "chunk_size": chunk_size,
                  "max_seconds": max_seconds, "cycle_steps": 20000,
                  "temperature_start": 12.0, "temperature_ratio": 0.015,
                  "scale_start": 0.2, "scale_ratio": 0.01,
                  "large_move_probability": 0.05, "large_move_scale": 0.3,
                  "radius_bounds": [0.25, 3.0]}
    for restart in range(restarts):
        p = initialize_ring(m, int(rng.integers(0, 2**31 - 1)))
        energy = sweep_energy(p, k)
        if energy < best_energy:
            best_energy, best_points = energy, p.copy()
            best_restart = restart
        for start in range(0, steps_per_restart, chunk_size):
            old_best = best_energy
            energy, best_energy, used, accepted = anneal_chunk(
                p, energy, best_points, best_energy, k, start,
                min(chunk_size, steps_per_restart - start),
                int(rng.integers(0, 2**31 - 1)))
            if best_energy < old_best:
                best_restart = restart
            total += used
            elapsed = time.monotonic() - started
            row = {"restart": restart, "step": start + used, "proposals": total,
                   "energy": float(energy), "best_energy": float(best_energy),
                   "accepted": int(accepted), "elapsed_seconds": elapsed}
            history.append(row)
            if callback:
                callback(row)
            if verbose and (start == 0 or (start + used) % 5000 == 0 or best_energy == 0):
                print(f"restart {restart + 1}/{restarts}, step {start + used:5d}, "
                      f"best E={best_energy:g}, elapsed {elapsed:.1f}s", flush=True)
            if best_energy == 0 or elapsed >= max_seconds:
                return SearchResult(best_points.copy(), float(best_energy), seed,
                                    best_restart, total, elapsed, history, parameters)
    return SearchResult(best_points.copy(), float(best_energy), seed, best_restart,
                        total, time.monotonic() - started, history, parameters)


def save_candidate(result, path, k=8):
    data = result.integer_candidate(k=k)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")
    return data
