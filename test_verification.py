"""Small independent checks of the geometric reduction and exact verifiers.

Run: python -m unittest -v test_verification
These checks are quick; full certificates are produced by the notebook.
"""
import json
from pathlib import Path
import unittest
import numpy as np

from annealing import initialize_ring, sweep_energy
import exact_verify
import slope_verify


def direct_chamber_minimum(reps, k):
    """Brute force full projection orders at all full-set event midpoints.

    Does not use the central rank/sign updates of the annealer.
    """
    points = np.vstack((reps, -reps))
    angles = []
    for i in range(len(points)):
        for j in range(i):
            d = points[i] - points[j]
            t = np.arctan2(d[1], d[0])
            angles.extend([t % np.pi, (t + np.pi / 2) % np.pi])
    events = np.unique(angles)
    # Floating duplicates of the same central event can differ by tiny roundoff.
    events = events[np.r_[True, np.diff(events) > 1e-11]]
    following = np.r_[events[1:], events[0] + np.pi]
    theta = (events + following) / 2
    normals = np.array([np.cos(theta), np.sin(theta)])
    perpendicular = np.array([-np.sin(theta), np.cos(theta)])
    cap_order = np.argsort(-(points @ normals), axis=0)
    median_order = np.argsort(-(points @ perpendicular), axis=0)
    return min(len(set(cap_order[:2*k, i]) & set(median_order[:len(reps), i]))
               for i in range(len(theta)))


class GeometryTests(unittest.TestCase):
    def test_sweep_against_full_projection_orders(self):
        for seed in [7,19,113]:
            reps = initialize_ring(8, seed)
            minimum = direct_chamber_minimum(reps, 2)
            self.assertEqual(sweep_energy(reps, 2) == 0, minimum >= 3)

    def test_known_counterexample_and_changed_representatives(self):
        data = json.loads((Path(__file__).parent/'data/paper_example.json').read_text())
        reps = np.array(data['points'][:48], dtype=float) / 1e7
        self.assertEqual(direct_chamber_minimum(reps, 8), 9)
        self.assertEqual(sweep_energy(reps, 8), 0)
        reps[::3] *= -1
        self.assertEqual(sweep_energy(reps[::-1].copy(), 8), 0)

    def test_both_verifiers_find_a_partition(self):
        square = [[1,1],[-1,1],[-1,-1],[1,-1]]
        a = exact_verify.verify(square, 1, full=True)
        b = slope_verify.verify(square, 1)
        self.assertEqual(a['status'], 'partition_found')
        self.assertEqual(b['status'], 'partition_found')

    def test_critical_assignments_against_exhaustive_subsets(self):
        # Four tied points at the median and two at the cap cutoff.
        import itertools
        points = [[-2,0],[-1,0],[1,0],[2,0],[0,2],[0,-2]]
        n = (1,0)
        af, at, ac = exact_verify.top_options(points, (0,1), 3)
        bf, bt, bc = exact_verify.top_options(points, n, 2)
        overlaps = {len((af | set(a)) & (bf | set(b)))
                    for a in itertools.combinations(at, ac)
                    for b in itertools.combinations(bt, bc)}
        self.assertEqual(exact_verify.event_overlap_range(points, n, 1),
                         (min(overlaps), max(overlaps)))
        self.assertEqual(slope_verify.event_possible(points, n, 1), 1 in overlaps)


if __name__ == '__main__':
    unittest.main()
