/**
 * Eigenvalue & Matrix Decomposition Engine
 * Provides QR decomposition, Jacobi eigenvalue algorithm, and power iteration.
 */

class EigenEngine {
  static qrDecomposition(matrix) {
    const n = matrix.length;
    const Q = Array.from({ length: n }, () => Array(n).fill(0));
    const R = Array.from({ length: n }, () => Array(n).fill(0));
    const v = matrix.map(row => [...row]);

    for (let j = 0; j < n; j++) {
      let norm = 0;
      for (let i = 0; i < n; i++) norm += v[i][j] * v[i][j];
      norm = Math.sqrt(norm);
      R[j][j] = norm;

      for (let i = 0; i < n; i++) {
        Q[i][j] = norm === 0 ? 0 : v[i][j] / norm;
      }

      for (let k = j + 1; k < n; k++) {
        let dot = 0;
        for (let i = 0; i < n; i++) dot += Q[i][j] * v[i][k];
        R[j][k] = dot;
        for (let i = 0; i < n; i++) v[i][k] -= dot * Q[i][j];
      }
    }
    return { Q, R };
  }

  static powerIteration(matrix, iterations = 100) {
    const n = matrix.length;
    let b = Array(n).fill(1 / Math.sqrt(n));

    for (let it = 0; it < iterations; it++) {
      const nextB = Array(n).fill(0);
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          nextB[i] += matrix[i][j] * b[j];
        }
      }
      let norm = Math.sqrt(nextB.reduce((acc, v) => acc + v * v, 0));
      b = nextB.map(v => v / (norm || 1));
    }

    let lambda = 0;
    for (let i = 0; i < n; i++) {
      let rowDot = 0;
      for (let j = 0; j < n; j++) rowDot += matrix[i][j] * b[j];
      lambda += b[i] * rowDot;
    }

    return { dominantEigenvalue: lambda, dominantEigenvector: b };
  }
}

if (typeof module !== 'undefined') module.exports = EigenEngine;
