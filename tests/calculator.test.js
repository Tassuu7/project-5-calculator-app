const ScientificEngine = require('../scientific.js');
const MemoryManager = require('../memory.js');
const SCIENTIFIC_CONSTANTS = require('../constants.js');

describe('Standard & Scientific Calculator Engine', () => {
  let engine;

  beforeEach(() => {
    engine = new ScientificEngine();
  });

  test('Basic arithmetic precedence', () => {
    expect(engine.evaluate('2 + 3 * 4')).toBe(14);
    expect(engine.evaluate('(2 + 3) * 4')).toBe(20);
    expect(engine.evaluate('100 / 4 - 5')).toBe(20);
  });

  test('Scientific powers and roots', () => {
    expect(engine.evaluate('2 ^ 3')).toBe(8);
    expect(engine.evaluate('sqrt(81)')).toBe(9);
    expect(engine.evaluate('cbrt(64)')).toBe(4);
  });

  test('Trigonometry degrees and radians', () => {
    engine.setAngleMode('DEG');
    expect(Math.abs(engine.evaluate('sin(30)') - 0.5)).toBeLessThan(1e-6);
    expect(Math.abs(engine.evaluate('tan(45)') - 1.0)).toBeLessThan(1e-6);

    engine.setAngleMode('RAD');
    expect(Math.abs(engine.evaluate('sin(pi / 2)') - 1.0)).toBeLessThan(1e-6);
  });

  test('Factorials and logarithms', () => {
    expect(engine.evaluate('5!')).toBe(120);
    expect(engine.evaluate('0!')).toBe(1);
    expect(Math.abs(engine.evaluate('log(100)') - 2)).toBeLessThan(1e-6);
    expect(Math.abs(engine.evaluate('ln(e)') - 1)).toBeLessThan(1e-6);
  });
});
