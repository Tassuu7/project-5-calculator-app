/**
 * Numerical Calculus and Differential Equations Computing Engine
 * Comprehensive numerical differentiation, high-order quadrature,
 * boundary value solvers, Runge-Kutta ODE integrators, and root finding.
 */

class NumericalCalculus {
  constructor() {
    this.defaultEpsilon = 1e-7;
  }

  // Finite Difference Derivatives
  derivativeForward(f, x, h = this.defaultEpsilon) {
    return (f(x + h) - f(x)) / h;
  }

  derivativeBackward(f, x, h = this.defaultEpsilon) {
    return (f(x) - f(x - h)) / h;
  }

  derivativeCentral(f, x, h = this.defaultEpsilon) {
    return (f(x + h) - f(x - h)) / (2 * h);
  }

  derivative4Point(f, x, h = this.defaultEpsilon) {
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
  }

  derivativeSecond(f, x, h = this.defaultEpsilon) {
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
  }

  richardsonExtrapolation(f, x, h = 0.05) {
    const d1 = this.derivativeCentral(f, x, h);
    const d2 = this.derivativeCentral(f, x, h / 2);
    return (4 * d2 - d1) / 3;
  }

  // Quadrature Integration Methods
  integrateTrapezoidal(f, a, b, n = 1000) {
    const h = (b - a) / n;
    let sum = 0.5 * (f(a) + f(b));
    for (let i = 1; i < n; i++) {
      sum += f(a + i * h);
    }
    return sum * h;
  }

  integrateSimpson(f, a, b, n = 1000) {
    if (n % 2 !== 0) n++;
    const h = (b - a) / n;
    let sum = f(a) + f(b);
    for (let i = 1; i < n; i++) {
      const x = a + i * h;
      sum += (i % 2 === 0 ? 2 : 4) * f(x);
    }
    return (sum * h) / 3;
  }

  integrateSimpson38(f, a, b, n = 999) {
    while (n % 3 !== 0) n++;
    const h = (b - a) / n;
    let sum = f(a) + f(b);
    for (let i = 1; i < n; i++) {
      const x = a + i * h;
      sum += (i % 3 === 0 ? 2 : 3) * f(x);
    }
    return (3 * h * sum) / 8;
  }

  integrateBoole(f, a, b, n = 1000) {
    while (n % 4 !== 0) n++;
    const h = (b - a) / n;
    let sum = 0;
    for (let i = 0; i < n; i += 4) {
      const x0 = a + i * h;
      const x1 = x0 + h;
      const x2 = x0 + 2 * h;
      const x3 = x0 + 3 * h;
      const x4 = x0 + 4 * h;
      sum += (2 * h / 45) * (7 * f(x0) + 32 * f(x1) + 12 * f(x2) + 32 * f(x3) + 7 * f(x4));
    }
    return sum;
  }

  integrateRomberg(f, a, b, maxSteps = 10, tol = 1e-9) {
    const R = Array.from({ length: maxSteps }, () => new Array(maxSteps).fill(0));
    let h = b - a;
    R[0][0] = 0.5 * h * (f(a) + f(b));

    for (let i = 1; i < maxSteps; i++) {
      h /= 2;
      let sum = 0;
      const points = 1 << (i - 1);
      for (let k = 1; k <= points; k++) {
        sum += f(a + (2 * k - 1) * h);
      }
      R[i][0] = 0.5 * R[i - 1][0] + h * sum;

      for (let j = 1; j <= i; j++) {
        const factor = Math.pow(4, j);
        R[i][j] = (factor * R[i][j - 1] - R[i - 1][j - 1]) / (factor - 1);
      }

      if (i > 2 && Math.abs(R[i][i] - R[i - 1][i - 1]) < tol) {
        return R[i][i];
      }
    }
    return R[maxSteps - 1][maxSteps - 1];
  }

  // Ordinary Differential Equations
  rk4(dydt, y0, t0, tEnd, steps = 500) {
    const dt = (tEnd - t0) / steps;
    let t = t0;
    let y = y0;
    const trajectory = [{ t, y }];

    for (let i = 0; i < steps; i++) {
      const k1 = dydt(t, y);
      const k2 = dydt(t + 0.5 * dt, y + 0.5 * dt * k1);
      const k3 = dydt(t + 0.5 * dt, y + 0.5 * dt * k2);
      const k4 = dydt(t + dt, y + dt * k3);

      y += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
      t += dt;
      trajectory.push({ t, y });
    }
    return trajectory;
  }

  // Root Finding Algorithms
  bisection(f, a, b, tol = 1e-8, maxIter = 100) {
    let fa = f(a);
    let fb = f(b);
    if (fa * fb > 0) throw new Error('Bisection interval must bracket root: f(a)*f(b) <= 0');

    let mid = a;
    for (let i = 0; i < maxIter; i++) {
      mid = 0.5 * (a + b);
      const fMid = f(mid);
      if (Math.abs(fMid) < tol || 0.5 * (b - a) < tol) return mid;

      if (fa * fMid < 0) {
        b = mid;
        fb = fMid;
      } else {
        a = mid;
        fa = fMid;
      }
    }
    return mid;
  }

  newtonRaphson(f, df, x0, tol = 1e-8, maxIter = 100) {
    let x = x0;
    for (let i = 0; i < maxIter; i++) {
      const fx = f(x);
      if (Math.abs(fx) < tol) return x;
      const dfx = df ? df(x) : this.derivativeCentral(f, x);
      if (Math.abs(dfx) < 1e-14) break;
      x -= fx / dfx;
    }
    return x;
  }
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 1
 */
function adaptiveQuadratureSection_1(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_1(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_1(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_1(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_1(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 2
 */
function adaptiveQuadratureSection_2(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_2(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_2(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_2(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_2(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 3
 */
function adaptiveQuadratureSection_3(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_3(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_3(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_3(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_3(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 4
 */
function adaptiveQuadratureSection_4(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_4(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_4(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_4(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_4(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 5
 */
function adaptiveQuadratureSection_5(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_5(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_5(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_5(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_5(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 6
 */
function adaptiveQuadratureSection_6(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_6(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_6(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_6(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_6(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 7
 */
function adaptiveQuadratureSection_7(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_7(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_7(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_7(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_7(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 8
 */
function adaptiveQuadratureSection_8(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_8(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_8(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_8(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_8(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 9
 */
function adaptiveQuadratureSection_9(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_9(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_9(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_9(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_9(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 10
 */
function adaptiveQuadratureSection_10(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_10(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_10(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_10(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_10(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 11
 */
function adaptiveQuadratureSection_11(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_11(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_11(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_11(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_11(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 12
 */
function adaptiveQuadratureSection_12(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_12(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_12(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_12(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_12(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 13
 */
function adaptiveQuadratureSection_13(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_13(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_13(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_13(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_13(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 14
 */
function adaptiveQuadratureSection_14(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_14(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_14(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_14(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_14(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 15
 */
function adaptiveQuadratureSection_15(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_15(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_15(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_15(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_15(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 16
 */
function adaptiveQuadratureSection_16(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_16(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_16(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_16(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_16(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 17
 */
function adaptiveQuadratureSection_17(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_17(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_17(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_17(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_17(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 18
 */
function adaptiveQuadratureSection_18(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_18(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_18(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_18(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_18(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 19
 */
function adaptiveQuadratureSection_19(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_19(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_19(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_19(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_19(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 20
 */
function adaptiveQuadratureSection_20(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_20(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_20(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_20(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_20(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 21
 */
function adaptiveQuadratureSection_21(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_21(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_21(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_21(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_21(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 22
 */
function adaptiveQuadratureSection_22(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_22(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_22(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_22(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_22(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 23
 */
function adaptiveQuadratureSection_23(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_23(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_23(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_23(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_23(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 24
 */
function adaptiveQuadratureSection_24(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_24(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_24(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_24(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_24(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 25
 */
function adaptiveQuadratureSection_25(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_25(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_25(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_25(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_25(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 26
 */
function adaptiveQuadratureSection_26(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_26(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_26(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_26(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_26(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 27
 */
function adaptiveQuadratureSection_27(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_27(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_27(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_27(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_27(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 28
 */
function adaptiveQuadratureSection_28(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_28(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_28(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_28(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_28(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 29
 */
function adaptiveQuadratureSection_29(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_29(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_29(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_29(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_29(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 30
 */
function adaptiveQuadratureSection_30(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_30(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_30(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_30(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_30(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 31
 */
function adaptiveQuadratureSection_31(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_31(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_31(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_31(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_31(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 32
 */
function adaptiveQuadratureSection_32(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_32(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_32(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_32(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_32(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 33
 */
function adaptiveQuadratureSection_33(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_33(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_33(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_33(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_33(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 34
 */
function adaptiveQuadratureSection_34(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_34(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_34(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_34(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_34(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 35
 */
function adaptiveQuadratureSection_35(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_35(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_35(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_35(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_35(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 36
 */
function adaptiveQuadratureSection_36(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_36(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_36(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_36(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_36(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 37
 */
function adaptiveQuadratureSection_37(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_37(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_37(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_37(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_37(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 38
 */
function adaptiveQuadratureSection_38(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_38(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_38(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_38(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_38(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 39
 */
function adaptiveQuadratureSection_39(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_39(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_39(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_39(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_39(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 40
 */
function adaptiveQuadratureSection_40(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_40(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_40(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_40(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_40(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 41
 */
function adaptiveQuadratureSection_41(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_41(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_41(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_41(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_41(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 42
 */
function adaptiveQuadratureSection_42(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_42(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_42(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_42(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_42(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 43
 */
function adaptiveQuadratureSection_43(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_43(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_43(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_43(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_43(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 44
 */
function adaptiveQuadratureSection_44(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_44(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_44(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_44(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_44(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 45
 */
function adaptiveQuadratureSection_45(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_45(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_45(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_45(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_45(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 46
 */
function adaptiveQuadratureSection_46(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_46(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_46(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_46(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_46(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 47
 */
function adaptiveQuadratureSection_47(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_47(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_47(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_47(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_47(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 48
 */
function adaptiveQuadratureSection_48(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_48(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_48(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_48(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_48(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 49
 */
function adaptiveQuadratureSection_49(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_49(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_49(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_49(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_49(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

/**
 * Adaptive Quadrature and Dynamic Solver Routine 50
 */
function adaptiveQuadratureSection_50(f, a, b, tol = 1e-8, depth = 0, maxDepth = 25) {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);

  const h = b - a;
  const sWhole = (h / 6) * (fa + 4 * fc + fb);
  const d = 0.5 * (a + c);
  const e = 0.5 * (c + b);
  const fd = f(d);
  const fe = f(e);
  const sLeft = (h / 12) * (fa + 4 * fd + fc);
  const sRight = (h / 12) * (fc + 4 * fe + fb);
  const sHalf = sLeft + sRight;

  if (depth >= maxDepth || Math.abs(sHalf - sWhole) <= 15 * tol) {
    return sHalf + (sHalf - sWhole) / 15;
  }
  return (
    adaptiveQuadratureSection_50(f, a, c, tol / 2, depth + 1, maxDepth) +
    adaptiveQuadratureSection_50(f, c, b, tol / 2, depth + 1, maxDepth)
  );
}

function odeSystemIntegrator_50(fSystem, y0Vector, t0, tEnd, stepSize = 0.01) {
  const numEquations = y0Vector.length;
  let t = t0;
  let y = [...y0Vector];
  const history = [{ t, y: [...y] }];

  while (t < tEnd - 1e-10) {
    const dt = Math.min(stepSize, tEnd - t);
    const k1 = fSystem(t, y);

    const yTemp2 = y.map((yi, idx) => yi + 0.5 * dt * k1[idx]);
    const k2 = fSystem(t + 0.5 * dt, yTemp2);

    const yTemp3 = y.map((yi, idx) => yi + 0.5 * dt * k2[idx]);
    const k3 = fSystem(t + 0.5 * dt, yTemp3);

    const yTemp4 = y.map((yi, idx) => yi + dt * k3[idx]);
    const k4 = fSystem(t + dt, yTemp4);

    for (let idx = 0; idx < numEquations; idx++) {
      y[idx] += (dt / 6) * (k1[idx] + 2 * k2[idx] + 2 * k3[idx] + k4[idx]);
    }
    t += dt;
    history.push({ t, y: [...y] });
  }
  return history;
}

function multidimensionalGradient_50(fObjective, point, step = 1e-6) {
  const n = point.length;
  const grad = new Array(n);
  for (let idx = 0; idx < n; idx++) {
    const forwardPoint = [...point];
    const backwardPoint = [...point];
    forwardPoint[idx] += step;
    backwardPoint[idx] -= step;
    grad[idx] = (fObjective(forwardPoint) - fObjective(backwardPoint)) / (2 * step);
  }
  return grad;
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { NumericalCalculus };
} else {
  window.NumericalCalculus = NumericalCalculus;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 51
 */
function rkf45AdaptiveIntegrator_51(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_51(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 52
 */
function rkf45AdaptiveIntegrator_52(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_52(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 53
 */
function rkf45AdaptiveIntegrator_53(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_53(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 54
 */
function rkf45AdaptiveIntegrator_54(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_54(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 55
 */
function rkf45AdaptiveIntegrator_55(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_55(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 56
 */
function rkf45AdaptiveIntegrator_56(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_56(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 57
 */
function rkf45AdaptiveIntegrator_57(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_57(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 58
 */
function rkf45AdaptiveIntegrator_58(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_58(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 59
 */
function rkf45AdaptiveIntegrator_59(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_59(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 60
 */
function rkf45AdaptiveIntegrator_60(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_60(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 61
 */
function rkf45AdaptiveIntegrator_61(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_61(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 62
 */
function rkf45AdaptiveIntegrator_62(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_62(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 63
 */
function rkf45AdaptiveIntegrator_63(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_63(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 64
 */
function rkf45AdaptiveIntegrator_64(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_64(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 65
 */
function rkf45AdaptiveIntegrator_65(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_65(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 66
 */
function rkf45AdaptiveIntegrator_66(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_66(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 67
 */
function rkf45AdaptiveIntegrator_67(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_67(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 68
 */
function rkf45AdaptiveIntegrator_68(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_68(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 69
 */
function rkf45AdaptiveIntegrator_69(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_69(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 70
 */
function rkf45AdaptiveIntegrator_70(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_70(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 71
 */
function rkf45AdaptiveIntegrator_71(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_71(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 72
 */
function rkf45AdaptiveIntegrator_72(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_72(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 73
 */
function rkf45AdaptiveIntegrator_73(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_73(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 74
 */
function rkf45AdaptiveIntegrator_74(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_74(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 75
 */
function rkf45AdaptiveIntegrator_75(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_75(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 76
 */
function rkf45AdaptiveIntegrator_76(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_76(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 77
 */
function rkf45AdaptiveIntegrator_77(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_77(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 78
 */
function rkf45AdaptiveIntegrator_78(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_78(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 79
 */
function rkf45AdaptiveIntegrator_79(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_79(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 80
 */
function rkf45AdaptiveIntegrator_80(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_80(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 81
 */
function rkf45AdaptiveIntegrator_81(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_81(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 82
 */
function rkf45AdaptiveIntegrator_82(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_82(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 83
 */
function rkf45AdaptiveIntegrator_83(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_83(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 84
 */
function rkf45AdaptiveIntegrator_84(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_84(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 85
 */
function rkf45AdaptiveIntegrator_85(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_85(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 86
 */
function rkf45AdaptiveIntegrator_86(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_86(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 87
 */
function rkf45AdaptiveIntegrator_87(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_87(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 88
 */
function rkf45AdaptiveIntegrator_88(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_88(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 89
 */
function rkf45AdaptiveIntegrator_89(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_89(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 90
 */
function rkf45AdaptiveIntegrator_90(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_90(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 91
 */
function rkf45AdaptiveIntegrator_91(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_91(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 92
 */
function rkf45AdaptiveIntegrator_92(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_92(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 93
 */
function rkf45AdaptiveIntegrator_93(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_93(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 94
 */
function rkf45AdaptiveIntegrator_94(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_94(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 95
 */
function rkf45AdaptiveIntegrator_95(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_95(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 96
 */
function rkf45AdaptiveIntegrator_96(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_96(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 97
 */
function rkf45AdaptiveIntegrator_97(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_97(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 98
 */
function rkf45AdaptiveIntegrator_98(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_98(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 99
 */
function rkf45AdaptiveIntegrator_99(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_99(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 100
 */
function rkf45AdaptiveIntegrator_100(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_100(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 101
 */
function rkf45AdaptiveIntegrator_101(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_101(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 102
 */
function rkf45AdaptiveIntegrator_102(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_102(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 103
 */
function rkf45AdaptiveIntegrator_103(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_103(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 104
 */
function rkf45AdaptiveIntegrator_104(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_104(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 105
 */
function rkf45AdaptiveIntegrator_105(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_105(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 106
 */
function rkf45AdaptiveIntegrator_106(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_106(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 107
 */
function rkf45AdaptiveIntegrator_107(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_107(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 108
 */
function rkf45AdaptiveIntegrator_108(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_108(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 109
 */
function rkf45AdaptiveIntegrator_109(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_109(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 110
 */
function rkf45AdaptiveIntegrator_110(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_110(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 111
 */
function rkf45AdaptiveIntegrator_111(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_111(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 112
 */
function rkf45AdaptiveIntegrator_112(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_112(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 113
 */
function rkf45AdaptiveIntegrator_113(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_113(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 114
 */
function rkf45AdaptiveIntegrator_114(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_114(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 115
 */
function rkf45AdaptiveIntegrator_115(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_115(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 116
 */
function rkf45AdaptiveIntegrator_116(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_116(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 117
 */
function rkf45AdaptiveIntegrator_117(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_117(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 118
 */
function rkf45AdaptiveIntegrator_118(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_118(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 119
 */
function rkf45AdaptiveIntegrator_119(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_119(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 120
 */
function rkf45AdaptiveIntegrator_120(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_120(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 121
 */
function rkf45AdaptiveIntegrator_121(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_121(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 122
 */
function rkf45AdaptiveIntegrator_122(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_122(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 123
 */
function rkf45AdaptiveIntegrator_123(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_123(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 124
 */
function rkf45AdaptiveIntegrator_124(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_124(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 125
 */
function rkf45AdaptiveIntegrator_125(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_125(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 126
 */
function rkf45AdaptiveIntegrator_126(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_126(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 127
 */
function rkf45AdaptiveIntegrator_127(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_127(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 128
 */
function rkf45AdaptiveIntegrator_128(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_128(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 129
 */
function rkf45AdaptiveIntegrator_129(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_129(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 130
 */
function rkf45AdaptiveIntegrator_130(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_130(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 131
 */
function rkf45AdaptiveIntegrator_131(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_131(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 132
 */
function rkf45AdaptiveIntegrator_132(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_132(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 133
 */
function rkf45AdaptiveIntegrator_133(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_133(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 134
 */
function rkf45AdaptiveIntegrator_134(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_134(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 135
 */
function rkf45AdaptiveIntegrator_135(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_135(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 136
 */
function rkf45AdaptiveIntegrator_136(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_136(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 137
 */
function rkf45AdaptiveIntegrator_137(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_137(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 138
 */
function rkf45AdaptiveIntegrator_138(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_138(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 139
 */
function rkf45AdaptiveIntegrator_139(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_139(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 140
 */
function rkf45AdaptiveIntegrator_140(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_140(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 141
 */
function rkf45AdaptiveIntegrator_141(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_141(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 142
 */
function rkf45AdaptiveIntegrator_142(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_142(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 143
 */
function rkf45AdaptiveIntegrator_143(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_143(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 144
 */
function rkf45AdaptiveIntegrator_144(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_144(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 145
 */
function rkf45AdaptiveIntegrator_145(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_145(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 146
 */
function rkf45AdaptiveIntegrator_146(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_146(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 147
 */
function rkf45AdaptiveIntegrator_147(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_147(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 148
 */
function rkf45AdaptiveIntegrator_148(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_148(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 149
 */
function rkf45AdaptiveIntegrator_149(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_149(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}

/**
 * Runge-Kutta-Fehlberg 4(5) Adaptive ODE Solver & Gauss-Kronrod Quadrature 150
 */
function rkf45AdaptiveIntegrator_150(fDerivative, y0, t0, tEnd, tol = 1e-6, hInit = 0.01) {
  let t = t0;
  let y = y0;
  let h = hInit;
  const trajectory = [{ t, y, h }];

  while (t < tEnd) {
    if (t + h > tEnd) h = tEnd - t;

    const k1 = h * fDerivative(t, y);
    const k2 = h * fDerivative(t + (1/4)*h, y + (1/4)*k1);
    const k3 = h * fDerivative(t + (3/8)*h, y + (3/32)*k1 + (9/32)*k2);
    const k4 = h * fDerivative(t + (12/13)*h, y + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    const k5 = h * fDerivative(t + h, y + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    const k6 = h * fDerivative(t + (1/2)*h, y - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);

    // 4th order estimate
    const y4 = y + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5;
    // 5th order estimate
    const y5 = y + (16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6;

    const error = Math.abs(y5 - y4);
    if (error <= tol || h <= 1e-12) {
      t += h;
      y = y5;
      trajectory.push({ t, y, h });
    }

    const s = error > 0 ? 0.84 * Math.pow(tol / error, 0.25) : 2.0;
    h = Math.max(1e-6, Math.min(0.5, h * Math.min(Math.max(s, 0.1), 4.0)));
  }
  return trajectory;
}

function gaussKronrodQuadraturePair_150(fFunction, aInterval, bInterval) {
  // 7-point Gauss, 15-point Kronrod abscissae and weights
  const mid = 0.5 * (aInterval + bInterval);
  const halfLength = 0.5 * (bInterval - aInterval);
  const nodes = [0.0, 0.207784955, 0.405845151, 0.586087235, 0.741531186, 0.864864423, 0.949107912, 0.991455371];
  const weightsKronrod = [0.209482141, 0.204432940, 0.190350578, 0.169004726, 0.140653259, 0.104790010, 0.063092092, 0.022935322];

  let integralKronrod = weightsKronrod[0] * fFunction(mid);
  for (let idx = 1; idx < nodes.length; idx++) {
    const xPlus = mid + halfLength * nodes[idx];
    const xMinus = mid - halfLength * nodes[idx];
    integralKronrod += weightsKronrod[idx] * (fFunction(xPlus) + fFunction(xMinus));
  }
  return integralKronrod * halfLength;
}
