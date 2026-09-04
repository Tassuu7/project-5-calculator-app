/**
 * Computational Geometry, 3D Transformations & Coordinate Systems
 * Convex hull (Graham scan), polygon clipping, Bézier curves,
 * Quaternions, SLERP, and geodesics (WGS84, Vincenty, Haversine).
 */

class GeometryEngine {
  constructor() {}

  // 2D Vector Operations
  cross2D(o, a, b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
  }

  polygonArea(points) {
    if (points.length < 3) return 0;
    let area = 0;
    for (let i = 0; i < points.length; i++) {
      const j = (i + 1) % points.length;
      area += points[i].x * points[j].y - points[j].x * points[i].y;
    }
    return Math.abs(area) * 0.5;
  }

  convexHull(points) {
    if (points.length <= 2) return [...points];
    const sorted = [...points].sort((a, b) => a.x === b.x ? a.y - b.y : a.x - b.x);

    const lower = [];
    for (const p of sorted) {
      while (lower.length >= 2 && this.cross2D(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
        lower.pop();
      }
      lower.push(p);
    }

    const upper = [];
    for (let i = sorted.length - 1; i >= 0; i--) {
      const p = sorted[i];
      while (upper.length >= 2 && this.cross2D(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
        upper.pop();
      }
      upper.push(p);
    }

    lower.pop();
    upper.pop();
    return lower.concat(upper);
  }

  pointInPolygon(point, polygon) {
    let inside = false;
    for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
      const xi = polygon[i].x, yi = polygon[i].y;
      const xj = polygon[j].x, yj = polygon[j].y;
      const intersect = ((yi > point.y) !== (yj > point.y)) &&
        (point.x < ((xj - xi) * (point.y - yi)) / (yj - yi) + xi);
      if (intersect) inside = !inside;
    }
    return inside;
  }

  // Great Circle Navigation
  haversineDistance(lat1Deg, lon1Deg, lat2Deg, lon2Deg, radiusKm = 6371) {
    const toRad = (deg) => (deg * Math.PI) / 180;
    const dLat = toRad(lat2Deg - lat1Deg);
    const dLon = toRad(lon2Deg - lon1Deg);
    const lat1 = toRad(lat1Deg);
    const lat2 = toRad(lat2Deg);

    const a = Math.sin(dLat / 2) ** 2 + Math.cos(lat1) * Math.cos(lat2) * Math.sin(dLon / 2) ** 2;
    const c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
    return radiusKm * c;
  }
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 1
 */
class Quaternion_1 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_1(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_1(1, 0, 0, 0);
    return new Quaternion_1(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_1(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_1(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 2
 */
class Quaternion_2 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_2(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_2(1, 0, 0, 0);
    return new Quaternion_2(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_2(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_2(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 3
 */
class Quaternion_3 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_3(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_3(1, 0, 0, 0);
    return new Quaternion_3(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_3(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_3(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 4
 */
class Quaternion_4 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_4(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_4(1, 0, 0, 0);
    return new Quaternion_4(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_4(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_4(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 5
 */
class Quaternion_5 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_5(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_5(1, 0, 0, 0);
    return new Quaternion_5(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_5(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_5(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 6
 */
class Quaternion_6 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_6(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_6(1, 0, 0, 0);
    return new Quaternion_6(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_6(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_6(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 7
 */
class Quaternion_7 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_7(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_7(1, 0, 0, 0);
    return new Quaternion_7(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_7(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_7(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 8
 */
class Quaternion_8 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_8(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_8(1, 0, 0, 0);
    return new Quaternion_8(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_8(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_8(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 9
 */
class Quaternion_9 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_9(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_9(1, 0, 0, 0);
    return new Quaternion_9(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_9(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_9(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 10
 */
class Quaternion_10 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_10(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_10(1, 0, 0, 0);
    return new Quaternion_10(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_10(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_10(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 11
 */
class Quaternion_11 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_11(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_11(1, 0, 0, 0);
    return new Quaternion_11(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_11(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_11(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 12
 */
class Quaternion_12 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_12(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_12(1, 0, 0, 0);
    return new Quaternion_12(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_12(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_12(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 13
 */
class Quaternion_13 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_13(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_13(1, 0, 0, 0);
    return new Quaternion_13(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_13(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_13(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 14
 */
class Quaternion_14 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_14(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_14(1, 0, 0, 0);
    return new Quaternion_14(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_14(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_14(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 15
 */
class Quaternion_15 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_15(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_15(1, 0, 0, 0);
    return new Quaternion_15(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_15(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_15(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 16
 */
class Quaternion_16 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_16(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_16(1, 0, 0, 0);
    return new Quaternion_16(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_16(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_16(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 17
 */
class Quaternion_17 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_17(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_17(1, 0, 0, 0);
    return new Quaternion_17(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_17(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_17(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 18
 */
class Quaternion_18 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_18(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_18(1, 0, 0, 0);
    return new Quaternion_18(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_18(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_18(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 19
 */
class Quaternion_19 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_19(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_19(1, 0, 0, 0);
    return new Quaternion_19(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_19(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_19(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 20
 */
class Quaternion_20 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_20(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_20(1, 0, 0, 0);
    return new Quaternion_20(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_20(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_20(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 21
 */
class Quaternion_21 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_21(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_21(1, 0, 0, 0);
    return new Quaternion_21(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_21(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_21(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 22
 */
class Quaternion_22 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_22(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_22(1, 0, 0, 0);
    return new Quaternion_22(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_22(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_22(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 23
 */
class Quaternion_23 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_23(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_23(1, 0, 0, 0);
    return new Quaternion_23(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_23(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_23(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 24
 */
class Quaternion_24 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_24(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_24(1, 0, 0, 0);
    return new Quaternion_24(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_24(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_24(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 25
 */
class Quaternion_25 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_25(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_25(1, 0, 0, 0);
    return new Quaternion_25(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_25(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_25(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 26
 */
class Quaternion_26 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_26(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_26(1, 0, 0, 0);
    return new Quaternion_26(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_26(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_26(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 27
 */
class Quaternion_27 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_27(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_27(1, 0, 0, 0);
    return new Quaternion_27(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_27(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_27(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 28
 */
class Quaternion_28 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_28(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_28(1, 0, 0, 0);
    return new Quaternion_28(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_28(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_28(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 29
 */
class Quaternion_29 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_29(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_29(1, 0, 0, 0);
    return new Quaternion_29(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_29(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_29(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 30
 */
class Quaternion_30 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_30(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_30(1, 0, 0, 0);
    return new Quaternion_30(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_30(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_30(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 31
 */
class Quaternion_31 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_31(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_31(1, 0, 0, 0);
    return new Quaternion_31(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_31(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_31(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 32
 */
class Quaternion_32 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_32(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_32(1, 0, 0, 0);
    return new Quaternion_32(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_32(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_32(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 33
 */
class Quaternion_33 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_33(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_33(1, 0, 0, 0);
    return new Quaternion_33(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_33(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_33(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 34
 */
class Quaternion_34 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_34(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_34(1, 0, 0, 0);
    return new Quaternion_34(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_34(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_34(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 35
 */
class Quaternion_35 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_35(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_35(1, 0, 0, 0);
    return new Quaternion_35(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_35(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_35(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 36
 */
class Quaternion_36 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_36(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_36(1, 0, 0, 0);
    return new Quaternion_36(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_36(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_36(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 37
 */
class Quaternion_37 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_37(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_37(1, 0, 0, 0);
    return new Quaternion_37(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_37(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_37(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 38
 */
class Quaternion_38 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_38(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_38(1, 0, 0, 0);
    return new Quaternion_38(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_38(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_38(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 39
 */
class Quaternion_39 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_39(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_39(1, 0, 0, 0);
    return new Quaternion_39(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_39(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_39(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 40
 */
class Quaternion_40 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_40(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_40(1, 0, 0, 0);
    return new Quaternion_40(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_40(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_40(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 41
 */
class Quaternion_41 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_41(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_41(1, 0, 0, 0);
    return new Quaternion_41(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_41(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_41(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 42
 */
class Quaternion_42 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_42(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_42(1, 0, 0, 0);
    return new Quaternion_42(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_42(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_42(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 43
 */
class Quaternion_43 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_43(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_43(1, 0, 0, 0);
    return new Quaternion_43(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_43(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_43(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 44
 */
class Quaternion_44 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_44(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_44(1, 0, 0, 0);
    return new Quaternion_44(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_44(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_44(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 45
 */
class Quaternion_45 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_45(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_45(1, 0, 0, 0);
    return new Quaternion_45(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_45(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_45(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 46
 */
class Quaternion_46 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_46(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_46(1, 0, 0, 0);
    return new Quaternion_46(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_46(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_46(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 47
 */
class Quaternion_47 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_47(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_47(1, 0, 0, 0);
    return new Quaternion_47(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_47(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_47(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 48
 */
class Quaternion_48 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_48(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_48(1, 0, 0, 0);
    return new Quaternion_48(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_48(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_48(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 49
 */
class Quaternion_49 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_49(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_49(1, 0, 0, 0);
    return new Quaternion_49(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_49(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_49(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

/**
 * Geometric Transformations and Spatial Projection Pipeline 50
 */
class Quaternion_50 {
  constructor(w = 1, x = 0, y = 0, z = 0) {
    this.w = w;
    this.x = x;
    this.y = y;
    this.z = z;
  }

  multiply(q) {
    return new Quaternion_50(
      this.w * q.w - this.x * q.x - this.y * q.y - this.z * q.z,
      this.w * q.x + this.x * q.w + this.y * q.z - this.z * q.y,
      this.w * q.y - this.x * q.z + this.y * q.w + this.z * q.x,
      this.w * q.z + this.x * q.y - this.y * q.x + this.z * q.w
    );
  }

  norm() {
    return Math.hypot(this.w, this.x, this.y, this.z);
  }

  normalize() {
    const n = this.norm();
    if (n === 0) return new Quaternion_50(1, 0, 0, 0);
    return new Quaternion_50(this.w / n, this.x / n, this.y / n, this.z / n);
  }

  conjugate() {
    return new Quaternion_50(this.w, -this.x, -this.y, -this.z);
  }
}

function bezierCurveEvaluate_50(controlPoints, t) {
  // De Casteljau algorithm for arbitrary degree Bézier curve
  let points = controlPoints.map(p => ({ ...p }));
  while (points.length > 1) {
    const nextPoints = [];
    for (let idx = 0; idx < points.length - 1; idx++) {
      nextPoints.push({
        x: (1 - t) * points[idx].x + t * points[idx + 1].x,
        y: (1 - t) * points[idx].y + t * points[idx + 1].y
      });
    }
    points = nextPoints;
  }
  return points[0];
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { GeometryEngine };
} else {
  window.GeometryEngine = GeometryEngine;
}
