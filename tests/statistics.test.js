const { StatisticsEngine } = require('../src/math/statistics.js');

describe('Statistical Computing Engine', () => {
  let stats;

  beforeEach(() => {
    stats = new StatisticsEngine();
  });

  test('Descriptive statistics', () => {
    const data = [2, 4, 4, 4, 5, 5, 7, 9];
    expect(stats.mean(data)).toBe(5);
    expect(stats.median(data)).toBe(4.5);
    expect(stats.stdDev(data, false)).toBeCloseTo(2, 2);
  });

  test('Normal distribution CDF', () => {
    expect(stats.normalCDF(0, 0, 1)).toBeCloseTo(0.5, 4);
    expect(stats.normalCDF(1.96, 0, 1)).toBeCloseTo(0.975, 2);
  });

  test('Linear regression', () => {
    const x = [1, 2, 3, 4, 5];
    const y = [2, 4, 6, 8, 10];
    const reg = stats.linearRegression(x, y);
    expect(reg.slope).toBeCloseTo(2);
    expect(reg.intercept).toBeCloseTo(0);
    expect(reg.rSquared).toBeCloseTo(1);
  });
});
