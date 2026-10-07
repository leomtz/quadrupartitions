"""Readable wrapper around two independent, standard-library exact verifiers."""
import hashlib
import itertools
import json
import time

import exact_verify
import slope_verify


def geometry_checks(points):
    """Use Python integers (unbounded arithmetic), including all triples."""
    if any(len(p) != 2 or any(type(v) is not int for v in p) for p in points):
        raise ValueError("Coordinates must be Python integers")
    tuples = [tuple(p) for p in points]
    unique = set(tuples)
    if len(unique) != len(tuples):
        raise ValueError("Repeated points")
    central = all((-x, -y) in unique for x, y in tuples)
    collinear = 0
    for a, b, c in itertools.combinations(tuples, 3):
        collinear += ((b[0] - a[0]) * (c[1] - a[1]) ==
                      (b[1] - a[1]) * (c[0] - a[0]))
    return {"distinct_points": len(unique), "centrally_symmetric": central,
            "triples_checked": len(tuples) * (len(tuples)-1) * (len(tuples)-2) // 6,
            "collinear_triples": collinear}


def certify(points, k, verbose=True):
    """Certify the rounded coordinates; do not use the floating search score."""
    started = time.monotonic()
    geometry = geometry_checks(points)
    if len(points) % 2 or not 0 < k <= len(points) // 4:
        raise ValueError("Require even N and 0 < k <= N/4")
    if verbose:
        print("Exact ray/chamber verifier...", flush=True)
    rays = exact_verify.verify(points, k, full=True)
    if verbose:
        print("Independent rational-slope verifier...", flush=True)
    slopes = slope_verify.verify(points, k)
    agreement = all(rays[key] == slopes[key] for key in
                    ("status", "N", "k", "points_sha256", "chamber_count",
                     "defect_histogram", "weak_solution_event_count"))
    certified = (agreement and rays["status"] == "counterexample_certified" and
                 geometry["centrally_symmetric"] and geometry["collinear_triples"] == 0)
    report = {"certified": certified, "verifiers_agree": agreement,
              "geometry": geometry, "rays": rays, "slopes": slopes,
              "elapsed_seconds": time.monotonic() - started}
    if verbose:
        print(f"Certified: {certified}; {rays['chamber_count']:,} chambers; "
              f"q-k histogram: {rays['defect_histogram']}; "
              f"weak event solutions: {rays['weak_solution_event_count']}; "
              f"elapsed {report['elapsed_seconds']:.1f}s", flush=True)
    return report


def save_report(report, path):
    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file")
    parser.add_argument("--output")
    args = parser.parse_args()
    data = json.load(open(args.file))
    result = certify(data["points"], data["k"])
    if args.output:
        save_report(result, args.output)
    raise SystemExit(0 if result["certified"] else 1)
