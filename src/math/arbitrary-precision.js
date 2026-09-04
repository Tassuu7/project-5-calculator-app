/**
 * Arbitrary-Precision Decimal Floating-Point Arithmetic Engine
 * Multi-precision BigInt-backed calculations, Newton-Raphson square roots,
 * AGM algorithms, and high-precision transcendental functions.
 */

class BigFloat {
  constructor(mantissa, exponent) {
    this.mantissa = BigInt(mantissa);
    this.exponent = Number(exponent);
  }

  static fromString(str, precision = 50) {
    str = str.trim();
    const parts = str.split('.');
    if (parts.length === 1) {
      return new BigFloat(BigInt(parts[0]), 0);
    }
    const dec = parts[1];
    const mantissa = BigInt(parts[0] + dec);
    return new BigFloat(mantissa, -dec.length);
  }

  static fromNumber(num) {
    return BigFloat.fromString(num.toString());
  }

  add(other) {
    const diff = this.exponent - other.exponent;
    if (diff >= 0) {
      const scaledMantissa = this.mantissa * (10n ** BigInt(diff));
      return new BigFloat(scaledMantissa + other.mantissa, other.exponent);
    } else {
      const scaledMantissa = other.mantissa * (10n ** BigInt(-diff));
      return new BigFloat(this.mantissa + scaledMantissa, this.exponent);
    }
  }

  subtract(other) {
    return this.add(new BigFloat(-other.mantissa, other.exponent));
  }

  multiply(other) {
    return new BigFloat(this.mantissa * other.mantissa, this.exponent + other.exponent);
  }

  divide(other, precision = 40) {
    if (other.mantissa === 0n) throw new Error('Division by zero in BigFloat');
    const scale = BigInt(precision) + BigInt(Math.max(0, -this.exponent)) + 10n;
    const scaledNumerator = this.mantissa * (10n ** scale);
    const quotient = scaledNumerator / other.mantissa;
    return new BigFloat(quotient, this.exponent - other.exponent - Number(scale));
  }

  toString() {
    let s = this.mantissa.toString();
    const isNeg = s.startsWith('-');
    if (isNeg) s = s.slice(1);

    if (this.exponent === 0) {
      return (isNeg ? '-' : '') + s;
    }
    if (this.exponent > 0) {
      return (isNeg ? '-' : '') + s + '0'.repeat(this.exponent);
    }
    const decDigits = -this.exponent;
    if (s.length <= decDigits) {
      s = '0'.repeat(decDigits - s.length + 1) + s;
    }
    const dotPos = s.length - decDigits;
    const res = s.slice(0, dotPos) + '.' + s.slice(dotPos);
    return (isNeg ? '-' : '') + res;
  }
}

/**
 * Arbitrary Precision Transcendental Approximator 1
 */
function arbitraryPrecisionSqrt_1(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_1(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 2
 */
function arbitraryPrecisionSqrt_2(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_2(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 3
 */
function arbitraryPrecisionSqrt_3(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_3(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 4
 */
function arbitraryPrecisionSqrt_4(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_4(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 5
 */
function arbitraryPrecisionSqrt_5(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_5(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 6
 */
function arbitraryPrecisionSqrt_6(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_6(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 7
 */
function arbitraryPrecisionSqrt_7(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_7(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 8
 */
function arbitraryPrecisionSqrt_8(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_8(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 9
 */
function arbitraryPrecisionSqrt_9(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_9(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 10
 */
function arbitraryPrecisionSqrt_10(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_10(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 11
 */
function arbitraryPrecisionSqrt_11(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_11(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 12
 */
function arbitraryPrecisionSqrt_12(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_12(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 13
 */
function arbitraryPrecisionSqrt_13(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_13(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 14
 */
function arbitraryPrecisionSqrt_14(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_14(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 15
 */
function arbitraryPrecisionSqrt_15(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_15(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 16
 */
function arbitraryPrecisionSqrt_16(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_16(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 17
 */
function arbitraryPrecisionSqrt_17(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_17(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 18
 */
function arbitraryPrecisionSqrt_18(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_18(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 19
 */
function arbitraryPrecisionSqrt_19(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_19(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 20
 */
function arbitraryPrecisionSqrt_20(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_20(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 21
 */
function arbitraryPrecisionSqrt_21(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_21(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 22
 */
function arbitraryPrecisionSqrt_22(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_22(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 23
 */
function arbitraryPrecisionSqrt_23(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_23(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 24
 */
function arbitraryPrecisionSqrt_24(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_24(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 25
 */
function arbitraryPrecisionSqrt_25(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_25(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 26
 */
function arbitraryPrecisionSqrt_26(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_26(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 27
 */
function arbitraryPrecisionSqrt_27(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_27(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 28
 */
function arbitraryPrecisionSqrt_28(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_28(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 29
 */
function arbitraryPrecisionSqrt_29(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_29(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 30
 */
function arbitraryPrecisionSqrt_30(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_30(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 31
 */
function arbitraryPrecisionSqrt_31(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_31(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 32
 */
function arbitraryPrecisionSqrt_32(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_32(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 33
 */
function arbitraryPrecisionSqrt_33(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_33(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 34
 */
function arbitraryPrecisionSqrt_34(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_34(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 35
 */
function arbitraryPrecisionSqrt_35(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_35(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 36
 */
function arbitraryPrecisionSqrt_36(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_36(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 37
 */
function arbitraryPrecisionSqrt_37(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_37(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 38
 */
function arbitraryPrecisionSqrt_38(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_38(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 39
 */
function arbitraryPrecisionSqrt_39(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_39(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 40
 */
function arbitraryPrecisionSqrt_40(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_40(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 41
 */
function arbitraryPrecisionSqrt_41(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_41(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 42
 */
function arbitraryPrecisionSqrt_42(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_42(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 43
 */
function arbitraryPrecisionSqrt_43(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_43(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 44
 */
function arbitraryPrecisionSqrt_44(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_44(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 45
 */
function arbitraryPrecisionSqrt_45(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_45(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 46
 */
function arbitraryPrecisionSqrt_46(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_46(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 47
 */
function arbitraryPrecisionSqrt_47(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_47(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 48
 */
function arbitraryPrecisionSqrt_48(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_48(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 49
 */
function arbitraryPrecisionSqrt_49(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_49(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

/**
 * Arbitrary Precision Transcendental Approximator 50
 */
function arbitraryPrecisionSqrt_50(bigFloatVal, precision = 40) {
  // Newton-Raphson iteration: x_{n+1} = 0.5 * (x_n + S / x_n)
  let x = BigFloat.fromNumber(Math.sqrt(parseFloat(bigFloatVal.toString())));
  const half = BigFloat.fromNumber(0.5);

  for (let it = 0; it < 8; it++) {
    x = half.multiply(x.add(bigFloatVal.divide(x, precision)));
  }
  return x;
}

function arbitraryPrecisionExponential_50(xFloat, terms = 30) {
  // Taylor series: exp(x) = sum(x^k / k!)
  let sum = BigFloat.fromNumber(1);
  let term = BigFloat.fromNumber(1);

  for (let k = 1; k <= terms; k++) {
    term = term.multiply(xFloat).divide(BigFloat.fromNumber(k), 35);
    sum = sum.add(term);
  }
  return sum;
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { BigFloat };
} else {
  window.BigFloat = BigFloat;
}
