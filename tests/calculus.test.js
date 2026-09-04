const { NumericalCalculus } = require('../src/math/calculus.js');

describe('Numerical Calculus Engine', () => {
  let calc;

  beforeEach(() => {
    calc = new NumericalCalculus();
  });

  test('Numerical differentiation', () => {
    const f = (x) => x * x * x; // f'(x) = 3x^2 => f'(2) = 12
    const deriv = calc.derivativeCentral(f, 2);
    expect(deriv).toBeCloseTo(12, 4);
  });

  test('Numerical integration via Simpson rule', () => {
    const f = (x) => x * x; // integral_0^3 x^2 dx = 9
    const area = calc.integrateSimpson(f, 0, 3, 100);
    expect(area).toBeCloseTo(9, 4);
  });

  test('Root finding via Newton-Raphson', () => {
    const f = (x) => x * x - 25; // root = 5
    const root = calc.newtonRaphson(f, null, 10);
    expect(root).toBeCloseTo(5, 5);
  });
});
