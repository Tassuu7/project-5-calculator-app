const { Vector, Matrix } = require('../src/math/algebra.js');

describe('Linear Algebra & Matrix Engine', () => {
  test('Vector addition, scaling, and dot product', () => {
    const v1 = new Vector([1, 2, 3]);
    const v2 = new Vector([4, 5, 6]);

    const sum = v1.add(v2);
    expect(sum.toArray()).toEqual([5, 7, 9]);

    const dot = v1.dot(v2);
    expect(dot).toBe(32);

    expect(v1.norm()).toBeCloseTo(Math.sqrt(14));
  });

  test('Matrix multiplication and determinant', () => {
    const m1 = new Matrix(2, 2, [[1, 2], [3, 4]]);
    const m2 = new Matrix(2, 2, [[2, 0], [1, 2]]);

    const prod = m1.multiply(m2);
    expect(prod.toArray()).toEqual([[4, 4], [10, 8]]);

    const det = m1.determinant();
    expect(det).toBeCloseTo(-2);
  });

  test('Matrix inverse', () => {
    const m = new Matrix(2, 2, [[4, 7], [2, 6]]);
    const inv = m.inverse();
    const identity = m.multiply(inv);

    expect(identity.get(0, 0)).toBeCloseTo(1);
    expect(identity.get(0, 1)).toBeCloseTo(0);
    expect(identity.get(1, 0)).toBeCloseTo(0);
    expect(identity.get(1, 1)).toBeCloseTo(1);
  });
});
