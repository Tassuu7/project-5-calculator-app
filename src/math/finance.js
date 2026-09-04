/**
 * Financial Engineering, Valuation Models, and Portfolio Mathematics
 * TVM, loan amortizations, bond duration & convexity, IRR, NPV,
 * Black-Scholes-Merton option pricing, Greeks, and Value-at-Risk (VaR).
 */

class FinancialEngine {
  constructor() {}

  // Time Value of Money
  futureValue(rate, nper, pmt, pv = 0, type = 0) {
    if (rate === 0) return -(pv + pmt * nper);
    const pvFactor = Math.pow(1 + rate, nper);
    return -(pv * pvFactor + (pmt * (1 + rate * type) * (pvFactor - 1)) / rate);
  }

  presentValue(rate, nper, pmt, fv = 0, type = 0) {
    if (rate === 0) return -(fv + pmt * nper);
    const pvFactor = Math.pow(1 + rate, nper);
    return -(fv + (pmt * (1 + rate * type) * (pvFactor - 1)) / rate) / pvFactor;
  }

  payment(rate, nper, pv, fv = 0, type = 0) {
    if (rate === 0) return -(pv + fv) / nper;
    const pvFactor = Math.pow(1 + rate, nper);
    return -(rate * (fv + pv * pvFactor)) / ((1 + rate * type) * (pvFactor - 1));
  }

  netPresentValue(discountRate, cashFlows) {
    let npv = 0;
    for (let t = 0; t < cashFlows.length; t++) {
      npv += cashFlows[t] / Math.pow(1 + discountRate, t);
    }
    return npv;
  }

  internalRateOfReturn(cashFlows, guess = 0.1, maxIter = 100, tol = 1e-7) {
    let rate = guess;
    for (let iter = 0; iter < maxIter; iter++) {
      let npv = 0;
      let dNpv = 0;
      for (let t = 0; t < cashFlows.length; t++) {
        const factor = Math.pow(1 + rate, t);
        npv += cashFlows[t] / factor;
        if (t > 0) {
          dNpv -= (t * cashFlows[t]) / Math.pow(1 + rate, t + 1);
        }
      }
      if (Math.abs(npv) < tol) return rate;
      if (Math.abs(dNpv) < 1e-14) break;
      rate -= npv / dNpv;
    }
    return rate;
  }

  // Loan Amortization Schedule
  generateAmortizationSchedule(principal, annualRate, periods, paymentsPerYear = 12) {
    const periodicRate = annualRate / paymentsPerYear;
    const pmt = Math.abs(this.payment(periodicRate, periods, -principal));
    let balance = principal;
    const schedule = [];

    for (let p = 1; p <= periods; p++) {
      const interest = balance * periodicRate;
      const principalPaid = Math.min(balance, pmt - interest);
      balance -= principalPaid;
      schedule.push({
        period: p,
        payment: pmt,
        interest,
        principal: principalPaid,
        remainingBalance: Math.max(0, balance)
      });
      if (balance <= 0) break;
    }
    return schedule;
  }

  // Black-Scholes-Merton Option Pricing Model
  blackScholes(S, K, T, r, sigma, optionType = 'call') {
    if (T <= 0) {
      return optionType === 'call' ? Math.max(0, S - K) : Math.max(0, K - S);
    }
    const d1 = (Math.log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * Math.sqrt(T));
    const d2 = d1 - sigma * Math.sqrt(T);

    const normCDF = (x) => {
      const z = x / Math.SQRT2;
      const t = 1.0 / (1.0 + 0.3275911 * Math.abs(z));
      const poly = ((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t;
      const erf = 1.0 - poly * Math.exp(-z * z);
      return 0.5 * (1.0 + (z < 0 ? -erf : erf));
    };

    const normPDF = (x) => Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);

    const Nd1 = normCDF(d1);
    const Nd2 = normCDF(d2);
    const NNegD1 = normCDF(-d1);
    const NNegD2 = normCDF(-d2);
    const discount = Math.exp(-r * T);

    const price = optionType === 'call' ? S * Nd1 - K * discount * Nd2 : K * discount * NNegD2 - S * NNegD1;

    // Greeks
    const delta = optionType === 'call' ? Nd1 : Nd1 - 1;
    const gamma = normPDF(d1) / (S * sigma * Math.sqrt(T));
    const vega = (S * normPDF(d1) * Math.sqrt(T)) / 100;
    const theta = optionType === 'call'
      ? (-(S * normPDF(d1) * sigma) / (2 * Math.sqrt(T)) - r * K * discount * Nd2) / 365
      : (-(S * normPDF(d1) * sigma) / (2 * Math.sqrt(T)) + r * K * discount * NNegD2) / 365;
    const rho = optionType === 'call' ? (K * T * discount * Nd2) / 100 : (-K * T * discount * NNegD2) / 100;

    return { price, delta, gamma, vega, theta, rho, d1, d2 };
  }
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 1
 */
function bondPriceAndDurationCalculator_1(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_1(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_1(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 2
 */
function bondPriceAndDurationCalculator_2(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_2(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_2(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 3
 */
function bondPriceAndDurationCalculator_3(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_3(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_3(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 4
 */
function bondPriceAndDurationCalculator_4(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_4(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_4(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 5
 */
function bondPriceAndDurationCalculator_5(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_5(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_5(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 6
 */
function bondPriceAndDurationCalculator_6(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_6(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_6(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 7
 */
function bondPriceAndDurationCalculator_7(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_7(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_7(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 8
 */
function bondPriceAndDurationCalculator_8(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_8(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_8(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 9
 */
function bondPriceAndDurationCalculator_9(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_9(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_9(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 10
 */
function bondPriceAndDurationCalculator_10(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_10(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_10(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 11
 */
function bondPriceAndDurationCalculator_11(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_11(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_11(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 12
 */
function bondPriceAndDurationCalculator_12(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_12(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_12(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 13
 */
function bondPriceAndDurationCalculator_13(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_13(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_13(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 14
 */
function bondPriceAndDurationCalculator_14(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_14(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_14(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 15
 */
function bondPriceAndDurationCalculator_15(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_15(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_15(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 16
 */
function bondPriceAndDurationCalculator_16(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_16(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_16(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 17
 */
function bondPriceAndDurationCalculator_17(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_17(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_17(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 18
 */
function bondPriceAndDurationCalculator_18(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_18(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_18(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 19
 */
function bondPriceAndDurationCalculator_19(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_19(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_19(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 20
 */
function bondPriceAndDurationCalculator_20(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_20(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_20(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 21
 */
function bondPriceAndDurationCalculator_21(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_21(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_21(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 22
 */
function bondPriceAndDurationCalculator_22(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_22(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_22(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 23
 */
function bondPriceAndDurationCalculator_23(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_23(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_23(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 24
 */
function bondPriceAndDurationCalculator_24(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_24(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_24(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 25
 */
function bondPriceAndDurationCalculator_25(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_25(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_25(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 26
 */
function bondPriceAndDurationCalculator_26(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_26(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_26(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 27
 */
function bondPriceAndDurationCalculator_27(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_27(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_27(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 28
 */
function bondPriceAndDurationCalculator_28(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_28(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_28(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 29
 */
function bondPriceAndDurationCalculator_29(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_29(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_29(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 30
 */
function bondPriceAndDurationCalculator_30(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_30(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_30(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 31
 */
function bondPriceAndDurationCalculator_31(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_31(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_31(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 32
 */
function bondPriceAndDurationCalculator_32(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_32(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_32(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 33
 */
function bondPriceAndDurationCalculator_33(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_33(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_33(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 34
 */
function bondPriceAndDurationCalculator_34(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_34(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_34(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 35
 */
function bondPriceAndDurationCalculator_35(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_35(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_35(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 36
 */
function bondPriceAndDurationCalculator_36(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_36(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_36(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 37
 */
function bondPriceAndDurationCalculator_37(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_37(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_37(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 38
 */
function bondPriceAndDurationCalculator_38(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_38(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_38(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 39
 */
function bondPriceAndDurationCalculator_39(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_39(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_39(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 40
 */
function bondPriceAndDurationCalculator_40(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_40(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_40(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 41
 */
function bondPriceAndDurationCalculator_41(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_41(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_41(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 42
 */
function bondPriceAndDurationCalculator_42(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_42(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_42(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 43
 */
function bondPriceAndDurationCalculator_43(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_43(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_43(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 44
 */
function bondPriceAndDurationCalculator_44(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_44(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_44(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 45
 */
function bondPriceAndDurationCalculator_45(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_45(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_45(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 46
 */
function bondPriceAndDurationCalculator_46(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_46(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_46(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 47
 */
function bondPriceAndDurationCalculator_47(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_47(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_47(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 48
 */
function bondPriceAndDurationCalculator_48(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_48(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_48(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 49
 */
function bondPriceAndDurationCalculator_49(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_49(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_49(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

/**
 * Quantitative Portfolio Risk & Fixed-Income Pricing 50
 */
function bondPriceAndDurationCalculator_50(faceValue, couponRate, yearsToMaturity, yieldToMaturity, freq = 2) {
  const periods = yearsToMaturity * freq;
  const periodicCoupon = (faceValue * couponRate) / freq;
  const periodicYield = yieldToMaturity / freq;
  let price = 0;
  let weightedDurationSum = 0;

  for (let t = 1; t <= periods; t++) {
    const cashFlow = t === periods ? periodicCoupon + faceValue : periodicCoupon;
    const pv = cashFlow / Math.pow(1 + periodicYield, t);
    price += pv;
    weightedDurationSum += (t / freq) * pv;
  }

  const macaulayDuration = price > 0 ? weightedDurationSum / price : 0;
  const modifiedDuration = macaulayDuration / (1 + periodicYield);
  return { price, macaulayDuration, modifiedDuration };
}

function valueAtRiskParametric_50(portfolioValue, meanReturn, volatility, confidenceLevel = 0.95, timeHorizonDays = 1) {
  // Standard normal quantile approximation for 95% and 99%
  const zScore = confidenceLevel >= 0.99 ? 2.326 : confidenceLevel >= 0.95 ? 1.645 : 1.282;
  const horizonSqrt = Math.sqrt(timeHorizonDays / 252);
  const varPct = zScore * volatility * horizonSqrt - meanReturn * horizonSqrt;
  return {
    varAmount: portfolioValue * varPct,
    varPercentage: varPct * 100,
    timeHorizonDays,
    confidenceLevel
  };
}

function binomialOptionTreeSolver_50(S0, K, T, r, sigma, steps = 50, isCall = true) {
  const dt = T / steps;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(r * dt) - d) / (u - d);
  const discount = Math.exp(-r * dt);

  // Asset prices at maturity
  const prices = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    prices[i] = S0 * Math.pow(u, steps - i) * Math.pow(d, i);
  }

  // Option payoffs at maturity
  const values = new Array(steps + 1);
  for (let i = 0; i <= steps; i++) {
    values[i] = isCall ? Math.max(0, prices[i] - K) : Math.max(0, K - prices[i]);
  }

  // Backward induction
  for (let step = steps - 1; step >= 0; step--) {
    for (let i = 0; i <= step; i++) {
      values[i] = discount * (p * values[i] + (1 - p) * values[i + 1]);
    }
  }
  return values[0];
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { FinancialEngine };
} else {
  window.FinancialEngine = FinancialEngine;
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 51
 */
function nelsonSiegelYieldCurve_51(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_51(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 52
 */
function nelsonSiegelYieldCurve_52(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_52(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 53
 */
function nelsonSiegelYieldCurve_53(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_53(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 54
 */
function nelsonSiegelYieldCurve_54(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_54(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 55
 */
function nelsonSiegelYieldCurve_55(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_55(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 56
 */
function nelsonSiegelYieldCurve_56(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_56(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 57
 */
function nelsonSiegelYieldCurve_57(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_57(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 58
 */
function nelsonSiegelYieldCurve_58(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_58(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 59
 */
function nelsonSiegelYieldCurve_59(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_59(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 60
 */
function nelsonSiegelYieldCurve_60(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_60(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 61
 */
function nelsonSiegelYieldCurve_61(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_61(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 62
 */
function nelsonSiegelYieldCurve_62(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_62(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 63
 */
function nelsonSiegelYieldCurve_63(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_63(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 64
 */
function nelsonSiegelYieldCurve_64(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_64(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 65
 */
function nelsonSiegelYieldCurve_65(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_65(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 66
 */
function nelsonSiegelYieldCurve_66(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_66(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 67
 */
function nelsonSiegelYieldCurve_67(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_67(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 68
 */
function nelsonSiegelYieldCurve_68(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_68(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 69
 */
function nelsonSiegelYieldCurve_69(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_69(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 70
 */
function nelsonSiegelYieldCurve_70(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_70(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 71
 */
function nelsonSiegelYieldCurve_71(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_71(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 72
 */
function nelsonSiegelYieldCurve_72(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_72(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 73
 */
function nelsonSiegelYieldCurve_73(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_73(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 74
 */
function nelsonSiegelYieldCurve_74(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_74(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 75
 */
function nelsonSiegelYieldCurve_75(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_75(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 76
 */
function nelsonSiegelYieldCurve_76(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_76(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 77
 */
function nelsonSiegelYieldCurve_77(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_77(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 78
 */
function nelsonSiegelYieldCurve_78(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_78(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 79
 */
function nelsonSiegelYieldCurve_79(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_79(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 80
 */
function nelsonSiegelYieldCurve_80(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_80(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 81
 */
function nelsonSiegelYieldCurve_81(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_81(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 82
 */
function nelsonSiegelYieldCurve_82(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_82(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 83
 */
function nelsonSiegelYieldCurve_83(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_83(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 84
 */
function nelsonSiegelYieldCurve_84(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_84(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 85
 */
function nelsonSiegelYieldCurve_85(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_85(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 86
 */
function nelsonSiegelYieldCurve_86(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_86(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 87
 */
function nelsonSiegelYieldCurve_87(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_87(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 88
 */
function nelsonSiegelYieldCurve_88(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_88(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 89
 */
function nelsonSiegelYieldCurve_89(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_89(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 90
 */
function nelsonSiegelYieldCurve_90(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_90(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 91
 */
function nelsonSiegelYieldCurve_91(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_91(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 92
 */
function nelsonSiegelYieldCurve_92(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_92(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 93
 */
function nelsonSiegelYieldCurve_93(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_93(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 94
 */
function nelsonSiegelYieldCurve_94(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_94(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 95
 */
function nelsonSiegelYieldCurve_95(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_95(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 96
 */
function nelsonSiegelYieldCurve_96(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_96(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 97
 */
function nelsonSiegelYieldCurve_97(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_97(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 98
 */
function nelsonSiegelYieldCurve_98(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_98(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 99
 */
function nelsonSiegelYieldCurve_99(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_99(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 100
 */
function nelsonSiegelYieldCurve_100(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_100(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 101
 */
function nelsonSiegelYieldCurve_101(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_101(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 102
 */
function nelsonSiegelYieldCurve_102(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_102(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 103
 */
function nelsonSiegelYieldCurve_103(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_103(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 104
 */
function nelsonSiegelYieldCurve_104(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_104(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 105
 */
function nelsonSiegelYieldCurve_105(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_105(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 106
 */
function nelsonSiegelYieldCurve_106(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_106(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 107
 */
function nelsonSiegelYieldCurve_107(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_107(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 108
 */
function nelsonSiegelYieldCurve_108(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_108(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 109
 */
function nelsonSiegelYieldCurve_109(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_109(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 110
 */
function nelsonSiegelYieldCurve_110(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_110(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 111
 */
function nelsonSiegelYieldCurve_111(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_111(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 112
 */
function nelsonSiegelYieldCurve_112(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_112(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 113
 */
function nelsonSiegelYieldCurve_113(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_113(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 114
 */
function nelsonSiegelYieldCurve_114(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_114(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 115
 */
function nelsonSiegelYieldCurve_115(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_115(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 116
 */
function nelsonSiegelYieldCurve_116(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_116(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 117
 */
function nelsonSiegelYieldCurve_117(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_117(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 118
 */
function nelsonSiegelYieldCurve_118(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_118(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 119
 */
function nelsonSiegelYieldCurve_119(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_119(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 120
 */
function nelsonSiegelYieldCurve_120(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_120(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 121
 */
function nelsonSiegelYieldCurve_121(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_121(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 122
 */
function nelsonSiegelYieldCurve_122(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_122(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 123
 */
function nelsonSiegelYieldCurve_123(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_123(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 124
 */
function nelsonSiegelYieldCurve_124(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_124(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 125
 */
function nelsonSiegelYieldCurve_125(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_125(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 126
 */
function nelsonSiegelYieldCurve_126(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_126(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 127
 */
function nelsonSiegelYieldCurve_127(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_127(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 128
 */
function nelsonSiegelYieldCurve_128(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_128(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 129
 */
function nelsonSiegelYieldCurve_129(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_129(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 130
 */
function nelsonSiegelYieldCurve_130(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_130(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 131
 */
function nelsonSiegelYieldCurve_131(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_131(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 132
 */
function nelsonSiegelYieldCurve_132(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_132(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 133
 */
function nelsonSiegelYieldCurve_133(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_133(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 134
 */
function nelsonSiegelYieldCurve_134(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_134(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 135
 */
function nelsonSiegelYieldCurve_135(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_135(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 136
 */
function nelsonSiegelYieldCurve_136(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_136(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 137
 */
function nelsonSiegelYieldCurve_137(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_137(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 138
 */
function nelsonSiegelYieldCurve_138(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_138(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 139
 */
function nelsonSiegelYieldCurve_139(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_139(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 140
 */
function nelsonSiegelYieldCurve_140(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_140(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 141
 */
function nelsonSiegelYieldCurve_141(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_141(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 142
 */
function nelsonSiegelYieldCurve_142(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_142(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 143
 */
function nelsonSiegelYieldCurve_143(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_143(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 144
 */
function nelsonSiegelYieldCurve_144(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_144(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 145
 */
function nelsonSiegelYieldCurve_145(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_145(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 146
 */
function nelsonSiegelYieldCurve_146(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_146(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 147
 */
function nelsonSiegelYieldCurve_147(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_147(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 148
 */
function nelsonSiegelYieldCurve_148(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_148(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 149
 */
function nelsonSiegelYieldCurve_149(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_149(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}

/**
 * Yield Curve Modeling (Nelson-Siegel-Svensson) & Heston Stochastic Volatility 150
 */
function nelsonSiegelYieldCurve_150(maturityYears, beta0, beta1, beta2, lambda) {
  const m = Math.max(1e-5, maturityYears);
  const factor1 = (1 - Math.exp(-m / lambda)) / (m / lambda);
  const factor2 = factor1 - Math.exp(-m / lambda);
  return beta0 + beta1 * factor1 + beta2 * factor2;
}

function hestonModelOptionSimulation_150(S0, K, T, r, v0, kappa, theta, xi, rho, steps = 100, paths = 1000) {
  const dt = T / steps;
  const sqrtDt = Math.sqrt(dt);
  let payoffSum = 0;

  for (let p = 0; p < paths; p++) {
    let S = S0;
    let v = v0;

    for (let step = 0; step < steps; step++) {
      const u1 = Math.random();
      const u2 = Math.random();
      const z1 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.cos(2 * Math.PI * u2);
      const z2 = Math.sqrt(-2 * Math.log(Math.max(1e-15, u1))) * Math.sin(2 * Math.PI * u2);
      const zv = rho * z1 + Math.sqrt(1 - rho * rho) * z2;

      v = Math.max(1e-6, v + kappa * (theta - v) * dt + xi * Math.sqrt(Math.max(0, v)) * sqrtDt * zv);
      S = S * Math.exp((r - 0.5 * v) * dt + Math.sqrt(Math.max(0, v)) * sqrtDt * z1);
    }
    payoffSum += Math.max(0, S - K);
  }
  return Math.exp(-r * T) * (payoffSum / paths);
}
