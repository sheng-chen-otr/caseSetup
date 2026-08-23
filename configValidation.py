"""Config validation helpers (fail-fast checks with actionable messages).

Currently focused on suspension hardpoint sanity checks. These functions are
pure (they operate on already-parsed dictionaries) so they are easy to unit
test without a filesystem or OpenFOAM.

The key check is per-corner unit consistency: a recurring authoring mistake is
writing one corner (e.g. FL) in millimetres and the others in metres while a
single ``HARDPOINT_SCALE`` is applied uniformly. The metric corners then
collapse to a few millimetres and the kinematic solver diverges. Detecting the
~1000x spread between corners up front turns a confusing solver divergence into
a clear message.
"""
import math


# Points measured relative to the wheel centre, used to estimate a corner's
# physical size (a suspension corner spans on the order of 0.1-0.5 m / 100-500 mm).
_SPAN_KEYS = (
    "uca_f_inner",
    "uca_r_inner",
    "lca_f_inner",
    "lca_r_inner",
    "tie_inner",
    "uca_outer_static",
    "lca_outer_static",
    "tie_outer_static",
    "pushrod_outer_static",
)


def _distance(a, b):
    return math.sqrt(sum((float(a[i]) - float(b[i])) ** 2 for i in range(3)))


def cornerCharacteristicSpan(corner):
    """Return a characteristic physical span for one parsed corner.

    Defined as the mean distance from the wheel centre to the available inner /
    outer hardpoints. Returns ``None`` if there is not enough data to measure a
    span (so callers can skip rather than divide by zero).
    """
    center = corner.get("wheel_center_static")
    if center is None:
        return None
    dists = [
        _distance(center, corner[key])
        for key in _SPAN_KEYS
        if key in corner
    ]
    dists = [d for d in dists if d > 0.0]
    if not dists:
        return None
    return sum(dists) / len(dists)


def checkHardpointUnitConsistency(corners, ratioThreshold=50.0):
    """Check that all corners share a consistent coordinate unit.

    :param corners: mapping of corner name -> parsed hardpoint dict (each value
        a mapping of point name -> ``[x, y, z]``), as produced by the CFG loader
        *before* ``HARDPOINT_SCALE`` is applied.
    :param ratioThreshold: max allowed ratio between the largest and smallest
        per-corner characteristic span before the corners are considered to be
        in mismatched units.
    :returns: ``(ok, message)`` where ``ok`` is ``True`` when the spans are
        consistent (or there is insufficient data to judge), and ``message`` is
        an actionable diagnostic string when ``ok`` is ``False`` (else ``""``).
    """
    spans = {}
    for name, corner in corners.items():
        span = cornerCharacteristicSpan(corner)
        if span is not None:
            spans[name] = span

    if len(spans) < 2:
        return True, ""

    smallestCorner = min(spans, key=spans.get)
    largestCorner = max(spans, key=spans.get)
    smallest = spans[smallestCorner]
    largest = spans[largestCorner]

    if smallest <= 0.0:
        return True, ""

    ratio = largest / smallest
    if ratio < ratioThreshold:
        return True, ""

    spanReport = ", ".join(
        "%s=%.4g" % (name, spans[name]) for name in sorted(spans)
    )
    message = (
        "Suspension hardpoint corners appear to be in mismatched units "
        "(largest/smallest corner span ratio %.0fx, threshold %.0fx). "
        "Corner spans: %s. This usually means one corner is in millimetres "
        "and another in metres while a single HARDPOINT_SCALE is applied to "
        "all of them. Put every corner in the SAME unit (all mm or all m) and "
        "set HARDPOINT_SCALE to match (0.001 for mm, 1.0 for m)."
        % (ratio, ratioThreshold, spanReport)
    )
    return False, message
