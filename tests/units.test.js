const { UnitConverter } = require('../units.js');

describe('Unit Conversion Suite', () => {
  let converter;

  beforeEach(() => {
    converter = new UnitConverter();
  });

  test('Length conversions', () => {
    expect(converter.convert(1, 'length', 'km', 'm')).toBe(1000);
    expect(converter.convert(1, 'length', 'in', 'cm')).toBeCloseTo(2.54, 4);
  });

  test('Temperature conversions', () => {
    expect(converter.convert(0, 'temperature', 'c', 'f')).toBe(32);
    expect(converter.convert(100, 'temperature', 'c', 'f')).toBe(212);
    expect(converter.convert(300, 'temperature', 'k', 'c')).toBeCloseTo(26.85, 2);
  });

  test('Digital storage conversions', () => {
    expect(converter.convert(1, 'digital', 'GB', 'MB')).toBe(1000);
    expect(converter.convert(1, 'digital', 'GiB', 'MiB')).toBe(1024);
  });
});
