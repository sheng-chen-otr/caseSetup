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

#mesh motion for internalSetup: crank-angle motion tables -> displacementInterpolation
#zones are flat faceZones normal to z (anchors fixed, body zones move with their body)
#a point's z displacement only depends on its z, so the zones must stay in the same z order

import os
import csv
import sys
import math

import numpy as np


#----- motion file

def readMotionFile(path):
    #csv with a theta column (deg) and one position column (m) per body, one cycle
    if not os.path.isfile(path):
        sys.exit('ERROR! Motion file %s not found!' % (path))
    with open(path, newline='') as f:
        rows = [r for r in csv.reader(f) if r and not r[0].strip().startswith('#')]
    header = [h.strip() for h in rows[0]]
    if 'theta' not in header:
        sys.exit('ERROR! Motion file %s needs a "theta" column (crank angle in deg)!' % (path))
    data = np.array([[float(v) for v in r] for r in rows[1:]], dtype=float)
    iTheta = header.index('theta')
    theta = data[:, iTheta]
    if np.any(np.diff(theta) <= 0):
        sys.exit('ERROR! Motion file %s: theta must increase strictly!' % (path))
    if theta[0] != 0.0 or theta[-1] > 360.0:
        sys.exit('ERROR! Motion file %s: theta must start at 0 and end at or below 360 deg!' % (path))
    columns = {h: data[:, i] for i, h in enumerate(header) if i != iTheta}
    return theta, columns


def periodicProfile(theta, pos, name, tol=1e-6):
    #closes the cycle at 360 deg; a given 360 value must match the 0 value
    stroke = pos.max() - pos.min()
    if theta[-1] == 360.0:
        if abs(pos[-1] - pos[0]) > tol + 1e-4*stroke:
            sys.exit('ERROR! Motion column %s is not periodic: %.6g at 0 deg, %.6g at 360 deg!'
                     % (name, pos[0], pos[-1]))
        return theta, pos
    return np.append(theta, 360.0), np.append(pos, pos[0])


def positionAt(angle, theta, pos):
    #periodic linear interpolation, angle in deg (any range)
    return np.interp(np.mod(angle, 360.0), theta, pos)


#----- zones

class Zone:
    #flat faceZone normal to z; body None = anchor (fixed)
    def __init__(self, name, patch, zRef, area, body=None):
        self.name = name
        self.patch = patch
        self.zRef = zRef
        self.area = area
        self.body = body


class Body:
    #rigid body moving along z with a crank-angle position profile
    #mode: how TDC_CLEARANCE is kept (soft / scale: the profile ends are adjusted, shift: the whole path is offset)
    def __init__(self, name, theta, pos, refAngle, clearance, mode='shift'):
        self.name = name
        self.theta = theta
        self.pos = np.array(pos, dtype=float)
        self.refAngle = refAngle
        self.clearance = clearance
        self.mode = mode
        self.offset = 0.0
        #the STL is drawn at the original profile value at REF_ANGLE, also after the profile is adjusted
        self.refPos = float(positionAt(refAngle, theta, self.pos))

    def shift(self, angle):
        #z shift of the body relative to its STL position (incl. clearance offset)
        return positionAt(angle, self.theta, self.pos) - self.refPos + self.offset


def zoneZ(zone, angle):
    if zone.body is None:
        return np.full_like(np.asarray(angle, dtype=float), zone.zRef)
    return zone.zRef + zone.body.shift(angle)


def flatZoneZ(name, tris, tol):
    #z of a flat zone from its triangles, error if it is not flat
    z = tris[:, :, 2]
    if z.max() - z.min() > tol:
        sys.exit('ERROR! Zone %s is not flat in z (%.6g to %.6g m): a motion zone must be a plane normal to z'
                 % (name, z.min(), z.max()))
    return float(z.mean())


def projectedArea(tris):
    #area projected on the xy plane
    a = tris[:, 1] - tris[:, 0]
    b = tris[:, 2] - tris[:, 0]
    return float(0.5*np.abs(a[:, 0]*b[:, 1] - a[:, 1]*b[:, 0]).sum())


def angleGrid(step):
    return np.arange(0.0, 360.0 + 0.5*step, step)


def sortedZones(zones, step):
    #z order of the zones: by their mean z over the cycle (a zone that briefly crosses
    #another is then still ordered sensibly and caught by checkZoneOrder)
    grid = angleGrid(step)
    return sorted(zones, key=lambda zn: float(np.mean(zoneZ(zn, grid))))


def zoneGaps(zones, step):
    #gaps between neighbouring zones over the cycle
    #returns list of (lower, upper, gap array over the grid)
    grid = angleGrid(step)
    order = sortedZones(zones, step)
    gaps = []
    for lower, upper in zip(order[:-1], order[1:]):
        gap = zoneZ(upper, grid) - zoneZ(lower, grid)
        gaps.append((lower, upper, gap))
    return gaps


def checkZoneOrder(zones, step):
    for lower, upper, gap in zoneGaps(zones, step):
        if gap.min() <= 0.0:
            iBad = int(np.argmin(gap))
            sys.exit('ERROR! Zones %s and %s cross at %.1f deg (gap %.6g m): the zone z order must not change '
                     'over the cycle' % (lower.name, upper.name, angleGrid(step)[iBad], gap.min()))


def applyClearance(zones, bodies, step, maxIter=10):
    #moves bodies away from fixed zones (or the other way round) until every gap that a body
    #closes stays >= its TDC_CLEARANCE over the cycle. returns {body: (offset, added volume)}
    report = {}
    for it in range(maxIter):
        changed = False
        for lower, upper, gap in zoneGaps(zones, step):
            if lower.body is upper.body:
                continue
            need = max(b.clearance for b in (lower.body, upper.body) if b is not None)
            short = need - gap.min()
            if short <= 1e-12:
                continue
            if lower.body is not None and upper.body is not None:
                sys.exit('ERROR! Bodies %s and %s come within %.6g m of each other (clearance %.6g m): '
                         'shift one of them in its motion column' % (lower.body.name, upper.body.name,
                                                                     gap.min(), need))
            body, zone = (upper.body, upper) if upper.body is not None else (lower.body, lower)
            sign = 1.0 if body is upper.body else -1.0
            body.offset += sign*short
            prev = report.get(body.name, (0.0, 0.0))
            report[body.name] = (prev[0] + sign*short, prev[1] + short*zone.area)
            changed = True
            break
        if not changed:
            return report
    sys.exit('ERROR! Clearance offsets did not settle (bodies pushed both ways?)')


#----- clearance by profile (CLEARANCE_MODE soft / scale)

def clearanceBounds(zones, body, step):
    #lowest and highest allowed shift of a body so every gap to a neighbouring anchor stays >= its TDC_CLEARANCE
    #returns (lower bound, upper bound), -inf / +inf where no anchor limits it
    lower, upper = -np.inf, np.inf
    for lo, up, gap in zoneGaps(zones, step):
        if up.body is body and lo.body is None:
            lower = max(lower, lo.zRef + body.clearance - up.zRef - body.offset)
        if lo.body is body and up.body is None:
            upper = min(upper, up.zRef - body.clearance - lo.zRef - body.offset)
    return lower, upper


def softFloor(p, low, raise_, high):
    #raises the lowest value low by raise_ with a smooth floor (x^4 + d^4)^(1/4), high stays where it is
    #the change dies out fast away from the end (d^4/(4 x^3)), the quartic correction keeps the far end exact and
    #the profile monotonic while the stroke is longer than raise_
    span = high - low
    corr = (span**4 + raise_**4)**0.25 - span
    return low + ((p - low)**4 + raise_**4)**0.25 - corr*((p - low)/span)**4


def softCeiling(p, high, lower_, low):
    #lowers the highest value high by lower_ with a smooth ceiling, low stays where it is
    span = high - low
    corr = (span**4 + lower_**4)**0.25 - span
    return high - ((high - p)**4 + lower_**4)**0.25 + corr*((high - p)/span)**4


def adjustProfiles(zones, bodies, step, margin=1e-9):
    #soft / scale bodies: the profile ends that come closer than TDC_CLEARANCE to an anchor are moved back, the
    #rest of the stroke is kept (soft: smooth floor/ceiling, only near that end; scale: the whole stroke scaled)
    #returns {body: report dict} for the bodies that were changed
    report = {}
    for body in bodies:
        if body.mode not in ('soft', 'scale'):
            continue
        lower, upper = clearanceBounds(zones, body, step)
        p = body.pos.copy()
        low, high = p.min(), p.max()
        lowTarget = lower + body.refPos + margin
        highTarget = upper + body.refPos - margin
        raise_ = lowTarget - low if lowTarget > low else 0.0
        lower_ = high - highTarget if highTarget < high else 0.0
        if raise_ == 0.0 and lower_ == 0.0:
            continue
        if high - lower_ <= low + raise_:
            sys.exit('ERROR! %s: its stroke does not fit between the anchors with TDC_CLEARANCE %g m on both sides'
                     % (body.name, body.clearance))
        if body.mode == 'scale':
            newLow, newHigh = low + raise_, high - lower_
            q = newLow + (p - low)*(newHigh - newLow)/(high - low)
        else:
            q = softFloor(p, low, raise_, high) if raise_ > 0 else p
            if lower_ > 0:
                q = softCeiling(q, high, lower_, q.min())
        dev = q - p
        k = int(np.argmax(np.abs(dev)))
        widths = np.gradient(body.theta)
        report[body.name] = {'mode': body.mode, 'raise': raise_, 'lower': lower_, 'strokeBefore': high - low,
                             'strokeAfter': q.max() - q.min(), 'maxDev': float(np.abs(dev).max()),
                             'devAngle': float(body.theta[k]),
                             'changedDeg': float(widths[np.abs(dev) > 2e-5].sum())}
        body.pos = q
    return report


#----- sealing of moving faces in a wall file (piston crown)

def sideWallRange(faceTris, wallTris, tol):
    #z range of the static walls around a moving face (the liner around a piston crown): walls that are not flat in
    #z and whose xy box covers most of the face's xy box. returns (zmin, zmax, names, wall points) or None
    fx, fy = faceTris[:, :, 0], faceTris[:, :, 1]
    faceBox = (fx.min(), fx.max(), fy.min(), fy.max())
    faceArea = max((faceBox[1] - faceBox[0])*(faceBox[3] - faceBox[2]), 1e-30)
    zmin, zmax, names, points = None, None, [], []
    for name, t in wallTris.items():
        z = t[:, :, 2]
        if z.max() - z.min() <= tol:
            continue
        box = (t[:, :, 0].min(), t[:, :, 0].max(), t[:, :, 1].min(), t[:, :, 1].max())
        dx = min(faceBox[1], box[1]) - max(faceBox[0], box[0])
        dy = min(faceBox[3], box[3]) - max(faceBox[2], box[2])
        #a liner has about the face's footprint (not a chamber wall that only encloses it)
        grow = 0.1*max(faceBox[1] - faceBox[0], faceBox[3] - faceBox[2])
        inside = (box[0] >= faceBox[0] - grow and box[1] <= faceBox[1] + grow and
                  box[2] >= faceBox[2] - grow and box[3] <= faceBox[3] + grow)
        if not inside or dx <= 0 or dy <= 0 or dx*dy < 0.5*faceArea:
            continue
        zmin = z.min() if zmin is None else min(zmin, z.min())
        zmax = z.max() if zmax is None else max(zmax, z.max())
        names.append(name)
        points.append(t.reshape(-1, 3))
    if not names:
        return None
    return zmin, zmax, names, np.vstack(points)


def pointsOutside2D(points, tris, tol):
    #xy points that lie outside every triangle (projected on xy) by more than tol
    a, b, c = tris[:, 0, :2], tris[:, 1, :2], tris[:, 2, :2]
    outside = []
    for point in points[:, :2]:
        #barycentric inside test for all triangles at once
        v0, v1, v2 = c - a, b - a, point - a
        d00, d01, d11 = (v0*v0).sum(1), (v0*v1).sum(1), (v1*v1).sum(1)
        d20, d21 = (v2*v0).sum(1), (v2*v1).sum(1)
        den = d00*d11 - d01*d01
        good = np.abs(den) > 1e-30
        u = np.where(good, (d11*d20 - d01*d21)/np.where(good, den, 1.0), -1.0)
        v = np.where(good, (d00*d21 - d01*d20)/np.where(good, den, 1.0), -1.0)
        if np.any((u >= -1e-12) & (v >= -1e-12) & (u + v <= 1.0 + 1e-12)):
            continue
        #distance to the triangle edges
        best = np.inf
        for e0, e1 in ((a, b), (b, c), (c, a)):
            ed = e1 - e0
            s = np.clip(((point - e0)*ed).sum(1)/np.maximum((ed*ed).sum(1), 1e-30), 0.0, 1.0)
            best = min(best, float(np.sqrt((((e0 + s[:, None]*ed) - point)**2).sum(1)).min()))
        if best > tol:
            outside.append(point)
    return outside


def sealProblems(faceTris, zFace, wallTris, tol):
    #a moving face in a wall file (piston crown) at height zFace must lie within its liner along z and cover the
    #bore across it, else snappy sees an open gap. returns a list of problem texts (empty = sealed, or no liner)
    liner = sideWallRange(faceTris, wallTris, tol)
    if liner is None:
        return []
    zmin, zmax, names, points = liner
    problems = []
    if zFace > zmax + tol or zFace < zmin - tol:
        problems.append('it sits at z %.6g m, outside its liner %s (z %.6g to %.6g m): extend the liner by at least '
                        '%.4g m, or mesh at an angle away from that end of the stroke'
                        % (zFace, ' '.join(names), zmin, zmax, max(zFace - zmax, zmin - zFace)))
    outside = pointsOutside2D(points, faceTris, tol)
    if outside:
        problems.append('it does not cover the bore of %s (%d liner points outside its outline): draw it slightly '
                        'larger than the bore' % (' '.join(names), len(outside)))
    return problems


def deformationRatio(zones, startAngle, step):
    #worst squash/stretch over the cycle of the gaps that change: max(g/g0, g0/g)
    worst, worstPair = 1.0, None
    grid = angleGrid(step)
    for lower, upper, gap in zoneGaps(zones, step):
        if lower.body is upper.body:
            continue
        g0 = float(zoneZ(upper, startAngle) - zoneZ(lower, startAngle))
        ratio = np.maximum(gap/g0, g0/gap).max()
        if ratio > worst:
            worst, worstPair = float(ratio), (lower.name, upper.name)
    return worst, worstPair


def segmentRatio(zones, angle0, angle1, step):
    #worst squash/stretch of the gaps that change between angle0 (meshing angle) and angle1: max(g/g0, g0/g)
    angles = np.append(np.arange(angle0, angle1, step), angle1)
    worst, worstPair = 1.0, None
    for lower, upper, gap in zoneGaps(zones, step):
        if lower.body is upper.body:
            continue
        g = zoneZ(upper, angles) - zoneZ(lower, angles)
        ratio = float(np.maximum(g/g[0], g[0]/g).max())
        if ratio > worst:
            worst, worstPair = ratio, (lower.name, upper.name)
    return worst, worstPair


def remeshAngles(zones, startAngle, step, limit):
    #remeshing: segments of the cycle, each meshed at its start angle and moved until any gap has been
    #squashed or stretched limit times. returns [(start, end)] in deg from startAngle to startAngle + 360
    segments = []
    angle0, cycleEnd = startAngle, startAngle + 360.0
    while angle0 < cycleEnd - 1e-9:
        angle1 = angle0
        while angle1 < cycleEnd - 1e-9:
            nxt = min(angle1 + step, cycleEnd)
            if segmentRatio(zones, angle0, nxt, step)[0] > limit:
                break
            angle1 = nxt
        if angle1 == angle0:
            sys.exit('ERROR! Remeshing: a gap is squashed or stretched more than REMESH_DEFORMATION %g within one '
                     'STEP_DEG step after %.1f deg: raise the limit or lower STEP_DEG' % (limit, angle0))
        segments.append((float(angle0), float(angle1)))
        angle0 = angle1
    return segments


def meshDisplacementZ(z, zones, angle):
    #z displacement of points from the geometry (STL) position to a crank angle as the motion solver moves the mesh:
    #linear in z between neighbouring zones (sorted at the STL position), constant beyond the outermost ones
    order = sorted(zones, key=lambda zn: zn.zRef)
    zRef = np.array([zn.zRef for zn in order])
    dz = np.array([float(zoneZ(zn, angle)) - zn.zRef for zn in order])
    return np.interp(np.asarray(z, dtype=float), zRef, dz)


def gapBands(zones, angle, step, cellOf, nCells, maxSplits=4):
    #z bands between neighbouring zones with fewer than nCells cells across at a meshing angle (end spaces at the
    #stroke ends), refined in z only. cellOf(zone) = cell size at the zone's surface level
    #returns [(lower, upper, zLo, zHi, splits, cells across before, after)]
    bands = []
    if nCells <= 0:
        return bands
    order = sortedZones(zones, step)
    for lower, upper in zip(order[:-1], order[1:]):
        if lower.body is upper.body and lower.body is not None:
            continue
        zLo, zHi = float(zoneZ(lower, angle)), float(zoneZ(upper, angle))
        cell = max(cellOf(lower), cellOf(upper))
        across = (zHi - zLo)/cell
        if across >= nCells:
            continue
        splits = min(int(np.ceil(np.log2(nCells/across))), maxSplits)
        bands.append((lower, upper, zLo, zHi, splits, across, across*2**splits))
    return bands


def locationAt(loc, zones, angle0, angle1, step, keep=(0.5, 0.5)):
    #LOC_IN_MESH at a meshing angle: in the same gap between zones as at angle0 (a point in an end space stays in it),
    #at the keep fraction of that gap (default the middle: snappy picks the region on the castellated mesh, before the
    #end-gap z refinement, so a point 1 mm from a face in a 3 mm gap of 2 mm cells landed in a one-cell pocket)
    order = sortedZones(zones, step)
    z = float(loc[2])
    for lower, upper in zip(order[:-1], order[1:]):
        lo0, hi0 = float(zoneZ(lower, angle0)), float(zoneZ(upper, angle0))
        if lo0 <= z <= hi0:
            frac = min(max((z - lo0)/(hi0 - lo0), keep[0]), keep[1])
            lo1, hi1 = float(zoneZ(lower, angle1)), float(zoneZ(upper, angle1))
            return [float(loc[0]), float(loc[1]), lo1 + frac*(hi1 - lo1)]
    return [float(v) for v in loc]


def autoStartAngle(zones, step):
    #start (= meshing) angle that minimises the worst squash/stretch ratio
    best = (None, None)
    for angle in angleGrid(step)[:-1]:
        ratio, pair = deformationRatio(zones, angle, step)
        if best[0] is None or ratio < best[0] - 1e-12:
            best = (ratio, angle)
    return float(best[1])


def bodyKinematics(body, freq, step):
    #peak velocity (m/s), peak acceleration (m/s2), largest velocity jump between samples
    grid = angleGrid(step)
    z = positionAt(grid, body.theta, body.pos)
    omega = 360.0*freq          #deg/s
    v = np.gradient(z, grid)*omega
    a = np.gradient(v, grid)*omega
    jump = np.abs(np.diff(v)).max()
    return {'stroke': float(z.max() - z.min()), 'vMax': float(np.abs(v).max()),
            'aMax': float(np.abs(a).max()), 'vJump': float(jump)}


#----- checks on the walls

def unanchoredHorizontalFaces(pidTris, zoneZs, zones, startAngle, step, tol):
    #horizontal wall faces strictly inside a stretching band move with the mesh: list them
    bands = [(float(zoneZ(lo, startAngle)), float(zoneZ(up, startAngle)))
             for lo, up, gap in zoneGaps(zones, step) if lo.body is not up.body]
    found = []
    for label, tris in pidTris.items():
        if len(tris) == 0:
            continue
        n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
        norm = np.linalg.norm(n, axis=1)
        horizontal = np.abs(n[:, 2]) > 0.99*np.maximum(norm, 1e-300)
        zc = tris[:, :, 2].mean(axis=1)
        for lo, up in bands:
            inside = horizontal & (zc > lo + tol) & (zc < up - tol)
            if inside.any() and not any(abs(zc[inside].mean() - z) < tol for z in zoneZs):
                found.append((label, float(zc[inside].min()), float(zc[inside].max())))
    return found


def pointTriangleDistance(points, tris, chunk=200):
    #minimum distance from each point to a set of triangles (exact, vectorised in chunks)
    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    ab, ac = b - a, c - a
    n = np.cross(ab, ac)
    nn = np.einsum('ij,ij->i', n, n)
    out = np.empty(len(points))
    for s in range(0, len(points), chunk):
        p = points[s:s+chunk, None, :]
        #distance to the plane where the projection falls inside the triangle
        ap = p - a
        t = np.einsum('pij,ij->pi', ap, n)/np.maximum(nn, 1e-300)
        q = p - t[..., None]*n
        v0, v1, v2 = ab, ac, q - a
        d00 = np.einsum('ij,ij->i', v0, v0)
        d01 = np.einsum('ij,ij->i', v0, v1)
        d11 = np.einsum('ij,ij->i', v1, v1)
        d20 = np.einsum('pij,ij->pi', v2, v0)
        d21 = np.einsum('pij,ij->pi', v2, v1)
        den = np.maximum(d00*d11 - d01*d01, 1e-300)
        bv = (d11*d20 - d01*d21)/den
        bw = (d00*d21 - d01*d20)/den
        inside = (bv >= 0) & (bw >= 0) & (bv + bw <= 1)
        dPlane = np.where(inside, np.abs(t)*np.sqrt(nn), np.inf)
        #otherwise the nearest edge
        dEdge = np.full(dPlane.shape, np.inf)
        for e0, e1 in ((a, b), (b, c), (c, a)):
            ed = e1 - e0
            lam = np.clip(np.einsum('pij,ij->pi', p - e0, ed)/np.maximum(np.einsum('ij,ij->i', ed, ed), 1e-300), 0, 1)
            proj = e0 + lam[..., None]*ed
            dEdge = np.minimum(dEdge, np.linalg.norm(p - proj, axis=2))
        out[s:s+chunk] = np.minimum(dPlane, dEdge).min(axis=1)
    return out


def gapWidth(bodyTris, wallTris, maxPoints=4000):
    #narrowest distance from the body gap faces (vertices + centroids) to the wall gap faces
    pts = np.vstack([bodyTris.reshape(-1, 3), bodyTris.mean(axis=1)])
    pts = np.unique(np.round(pts, 12), axis=0)
    if len(pts) > maxPoints:
        pts = pts[np.linspace(0, len(pts) - 1, maxPoints).astype(int)]
    return float(pointTriangleDistance(pts, wallTris).min())


def gapRefinement(gap, baseCell, gapCells):
    #snappy level so the gap has at least gapCells cells, and the distance of the refinement band
    level = max(0, int(math.ceil(math.log2(baseCell*gapCells/gap))))
    return level, 1.25*gap, baseCell/2**level


#----- tables

def writeTables(bodies, startAngle, freq, nCycles, step, tableDir):
    #constant/tables: zeroTable (anchors) + <body>Table, identical time points
    os.makedirs(tableDir, exist_ok=True)
    nSteps = int(round(nCycles*360.0/step))
    angles = startAngle + step*np.arange(nSteps + 1)
    times = (angles - startAngle)/(360.0*freq)
    tables = {'zeroTable': np.zeros_like(times)}
    for body in bodies:
        tables['%sTable' % (body.name)] = (positionAt(angles, body.theta, body.pos)
                                           - positionAt(startAngle, body.theta, body.pos))
    for name, dz in tables.items():
        with open(os.path.join(tableDir, name), 'w') as f:
            f.write('FoamFile { version 2.0; format ascii; class scalarVectorTable; object %s; }\n' % (name))
            f.write('%d\n(\n' % (len(times)))
            for t, d in zip(times, dz):
                f.write('    (%.12g (0 0 %.12g))\n' % (t, d))
            f.write(')\n')
    return list(tables.keys())


def writeCrankTables(bodies, startAngle, step, tableDir):
    #writing constant/tables/<body>Crank for DRIVE crank, displacement vs crank angle over one cycle
    os.makedirs(tableDir, exist_ok=True)
    #the table has to end at exactly 360 (outOfBounds repeat), the start angle is a grid point so it's zero there
    angles = np.linspace(0.0, 360.0, int(math.ceil(360.0/step - 1e-9)) + 1)
    angles = np.unique(np.append(angles, np.mod(startAngle, 360.0)))
    names = []
    for body in bodies:
        dz = positionAt(angles, body.theta, body.pos) - positionAt(startAngle, body.theta, body.pos)
        name = '%sCrank' % (body.name)
        with open(os.path.join(tableDir, name), 'w') as f:
            f.write('%d\n(\n' % (len(angles)))
            for a, d in zip(angles, dz):
                f.write('    (%.12g %.12g)\n' % (a, d))
            f.write(')\n')
        names.append(name)
    return names


def writeStlTriangles(path, pidTris):
    #ascii stl, one solid per PID
    with open(path, 'w') as f:
        for pid, tris in pidTris.items():
            f.write('solid %s\n' % (pid))
            for t in tris:
                n = np.cross(t[1] - t[0], t[2] - t[0])
                nl = np.linalg.norm(n)
                n = n/nl if nl > 0 else n
                f.write(' facet normal %g %g %g\n  outer loop\n' % tuple(n))
                for p in t:
                    f.write('   vertex %.9g %.9g %.9g\n' % tuple(p))
                f.write('  endloop\n endfacet\n')
            f.write('endsolid %s\n' % (pid))
