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
    def __init__(self, name, theta, pos, refAngle, clearance):
        self.name = name
        self.theta = theta
        self.pos = pos
        self.refAngle = refAngle
        self.clearance = clearance
        self.offset = 0.0

    def shift(self, angle):
        #z shift of the body relative to its STL position (incl. clearance offset)
        return positionAt(angle, self.theta, self.pos) - positionAt(self.refAngle, self.theta, self.pos) + self.offset


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
    #constant/tables/<body>Crank for DRIVE crank: (crank angle deg, displacement m) over one cycle,
    #0 to 360 deg (read with outOfBounds repeat), zero at the start angle
    os.makedirs(tableDir, exist_ok=True)
    #the period is the table range, so the grid must end at exactly 360; the start angle is a grid
    #point so the interpolated displacement is exactly zero there
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
