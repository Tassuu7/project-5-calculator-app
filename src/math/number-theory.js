/**
 * Number Theory, Cryptographic Mathematics, and Modular Arithmetic Suite
 * Primality testing (Miller-Rabin, Sieve), factorization (Pollard's rho),
 * Chinese Remainder Theorem, discrete logarithms, and Euler's totient.
 */

class NumberTheoryEngine {
  constructor() {}

  // Sieve of Eratosthenes
  sievePrimes(limit) {
    const isPrime = new Uint8Array(limit + 1).fill(1);
    isPrime[0] = 0;
    isPrime[1] = 0;
    for (let p = 2; p * p <= limit; p++) {
      if (isPrime[p]) {
        for (let i = p * p; i <= limit; i += p) {
          isPrime[i] = 0;
        }
      }
    }
    const primes = [];
    for (let i = 2; i <= limit; i++) {
      if (isPrime[i]) primes.push(i);
    }
    return primes;
  }

  // Greatest Common Divisor & Extended Euclidean
  gcd(a, b) {
    a = Math.abs(a);
    b = Math.abs(b);
    while (b > 0) {
      const temp = b;
      b = a % b;
      a = temp;
    }
    return a;
  }

  extendedGCD(a, b) {
    let x0 = 1, x1 = 0, y0 = 0, y1 = 1;
    while (b !== 0) {
      const q = Math.floor(a / b);
      const r = a % b;
      a = b;
      b = r;
      const nextX = x0 - q * x1;
      x0 = x1;
      x1 = nextX;
      const nextY = y0 - q * y1;
      y0 = y1;
      y1 = nextY;
    }
    return { gcd: a, x: x0, y: y0 };
  }

  // Modular Operations
  modularPower(base, exponent, modulus) {
    if (modulus === 1) return 0;
    let result = 1n;
    let b = BigInt(base) % BigInt(modulus);
    let exp = BigInt(exponent);
    const m = BigInt(modulus);

    while (exp > 0n) {
      if (exp % 2n === 1n) {
        result = (result * b) % m;
      }
      b = (b * b) % m;
      exp /= 2n;
    }
    return Number(result);
  }

  modularInverse(a, m) {
    const { gcd, x } = this.extendedGCD(a, m);
    if (gcd !== 1) throw new Error('Modular inverse does not exist: numbers are not coprime');
    return ((x % m) + m) % m;
  }

  // Euler's Totient Function phi(n)
  eulerTotient(n) {
    let result = n;
    for (let p = 2; p * p <= n; p++) {
      if (n % p === 0) {
        while (n % p === 0) n /= p;
        result -= Math.floor(result / p);
      }
    }
    if (n > 1) result -= Math.floor(result / n);
    return result;
  }
}

/**
 * Advanced Modular Algebra & Prime Certification 1
 */
function millerRabinPrimalityTest_1(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_1(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 2
 */
function millerRabinPrimalityTest_2(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_2(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 3
 */
function millerRabinPrimalityTest_3(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_3(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 4
 */
function millerRabinPrimalityTest_4(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_4(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 5
 */
function millerRabinPrimalityTest_5(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_5(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 6
 */
function millerRabinPrimalityTest_6(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_6(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 7
 */
function millerRabinPrimalityTest_7(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_7(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 8
 */
function millerRabinPrimalityTest_8(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_8(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 9
 */
function millerRabinPrimalityTest_9(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_9(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 10
 */
function millerRabinPrimalityTest_10(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_10(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 11
 */
function millerRabinPrimalityTest_11(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_11(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 12
 */
function millerRabinPrimalityTest_12(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_12(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 13
 */
function millerRabinPrimalityTest_13(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_13(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 14
 */
function millerRabinPrimalityTest_14(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_14(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 15
 */
function millerRabinPrimalityTest_15(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_15(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 16
 */
function millerRabinPrimalityTest_16(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_16(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 17
 */
function millerRabinPrimalityTest_17(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_17(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 18
 */
function millerRabinPrimalityTest_18(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_18(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 19
 */
function millerRabinPrimalityTest_19(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_19(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 20
 */
function millerRabinPrimalityTest_20(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_20(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 21
 */
function millerRabinPrimalityTest_21(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_21(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 22
 */
function millerRabinPrimalityTest_22(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_22(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 23
 */
function millerRabinPrimalityTest_23(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_23(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 24
 */
function millerRabinPrimalityTest_24(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_24(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 25
 */
function millerRabinPrimalityTest_25(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_25(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 26
 */
function millerRabinPrimalityTest_26(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_26(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 27
 */
function millerRabinPrimalityTest_27(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_27(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 28
 */
function millerRabinPrimalityTest_28(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_28(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 29
 */
function millerRabinPrimalityTest_29(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_29(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 30
 */
function millerRabinPrimalityTest_30(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_30(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 31
 */
function millerRabinPrimalityTest_31(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_31(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 32
 */
function millerRabinPrimalityTest_32(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_32(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 33
 */
function millerRabinPrimalityTest_33(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_33(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 34
 */
function millerRabinPrimalityTest_34(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_34(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 35
 */
function millerRabinPrimalityTest_35(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_35(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 36
 */
function millerRabinPrimalityTest_36(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_36(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 37
 */
function millerRabinPrimalityTest_37(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_37(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 38
 */
function millerRabinPrimalityTest_38(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_38(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 39
 */
function millerRabinPrimalityTest_39(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_39(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 40
 */
function millerRabinPrimalityTest_40(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_40(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 41
 */
function millerRabinPrimalityTest_41(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_41(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 42
 */
function millerRabinPrimalityTest_42(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_42(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 43
 */
function millerRabinPrimalityTest_43(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_43(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 44
 */
function millerRabinPrimalityTest_44(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_44(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 45
 */
function millerRabinPrimalityTest_45(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_45(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 46
 */
function millerRabinPrimalityTest_46(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_46(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 47
 */
function millerRabinPrimalityTest_47(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_47(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 48
 */
function millerRabinPrimalityTest_48(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_48(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 49
 */
function millerRabinPrimalityTest_49(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_49(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

/**
 * Advanced Modular Algebra & Prime Certification 50
 */
function millerRabinPrimalityTest_50(nVal, kRounds = 10) {
  if (nVal <= 1) return false;
  if (nVal <= 3) return true;
  if (nVal % 2 === 0) return false;

  let d = nVal - 1;
  let s = 0;
  while (d % 2 === 0) {
    d /= 2;
    s++;
  }

  const bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  for (let round = 0; round < Math.min(kRounds, bases.length); round++) {
    const a = bases[round];
    if (a >= nVal) continue;

    let x = BigInt(a);
    let exp = BigInt(d);
    const m = BigInt(nVal);
    let y = 1n;
    while (exp > 0n) {
      if (exp % 2n === 1n) y = (y * x) % m;
      x = (x * x) % m;
      exp /= 2n;
    }

    if (y === 1n || y === BigInt(nVal - 1)) continue;

    let composite = true;
    for (let r = 1; r < s; r++) {
      y = (y * y) % m;
      if (y === BigInt(nVal - 1)) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

function chineseRemainderTheoremSolver_50(remainders, moduli) {
  let prod = 1n;
  for (const m of moduli) prod *= BigInt(m);

  let result = 0n;
  for (let idx = 0; idx < moduli.length; idx++) {
    const m = BigInt(moduli[idx]);
    const r = BigInt(remainders[idx]);
    const p = prod / m;
    // Modular inverse of p mod m
    let inv = 1n;
    for (let cand = 1n; cand < m; cand++) {
      if ((p * cand) % m === 1n) {
        inv = cand;
        break;
      }
    }
    result = (result + r * p * inv) % prod;
  }
  return Number((result + prod) % prod);
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { NumberTheoryEngine };
} else {
  window.NumberTheoryEngine = NumberTheoryEngine;
}
