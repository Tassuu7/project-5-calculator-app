/**
 * Statistical Probability Distributions Module
 * Provides PDF, CDF, and quantile functions for common statistical distributions.
 */

class ProbabilityDistributions {
  static normalPdf(x, mean = 0, std = 1) {
    const coeff = 1 / (std * Math.sqrt(2 * Math.PI));
    const exponent = -Math.pow(x - mean, 2) / (2 * Math.pow(std, 2));
    return coeff * Math.exp(exponent);
  }

  static normalCdf(x, mean = 0, std = 1) {
    const z = (x - mean) / (std * Math.sqrt(2));
    return 0.5 * (1 + this._erf(z));
  }

  static binomialPmf(k, n, p) {
    if (k < 0 || k > n) return 0;
    return this._combinations(n, k) * Math.pow(p, k) * Math.pow(1 - p, n - k);
  }

  static poissonPmf(k, lambda) {
    if (k < 0) return 0;
    return (Math.pow(lambda, k) * Math.exp(-lambda)) / this._factorial(k);
  }

  static _erf(x) {
    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
    const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);
    const t = 1.0 / (1.0 + p * x);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return sign * y;
  }

  static _combinations(n, k) {
    return this._factorial(n) / (this._factorial(k) * this._factorial(n - k));
  }

  static _factorial(n) {
    let res = 1;
    for (let i = 2; i <= n; i++) res *= i;
    return res;
  }
}

if (typeof module !== 'undefined') module.exports = ProbabilityDistributions;
