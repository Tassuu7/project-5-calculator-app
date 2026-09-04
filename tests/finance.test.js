const { FinancialEngine } = require('../src/math/finance.js');

describe('Financial Engineering Engine', () => {
  let fin;

  beforeEach(() => {
    fin = new FinancialEngine();
  });

  test('Time Value of Money (FV, PV, PMT)', () => {
    const fv = fin.futureValue(0.05, 10, 0, -1000);
    expect(fv).toBeCloseTo(1628.89, 1);

    const pv = fin.presentValue(0.05, 10, 0, fv);
    expect(pv).toBeCloseTo(-1000, 1);
  });

  test('Net Present Value (NPV) & IRR', () => {
    const cashFlows = [-1000, 300, 420, 680];
    const npv = fin.netPresentValue(0.10, cashFlows);
    expect(npv).toBeGreaterThan(0);

    const irr = fin.internalRateOfReturn(cashFlows);
    expect(irr).toBeGreaterThan(0.10);
  });

  test('Black-Scholes Option Pricing', () => {
    const call = fin.blackScholes(100, 100, 1, 0.05, 0.2, 'call');
    expect(call.price).toBeGreaterThan(10);
    expect(call.price).toBeLessThan(11);
    expect(call.delta).toBeGreaterThan(0.6);
  });
});
