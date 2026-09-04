/**
 * Comprehensive Verification Test Suite
 * Tests ScientificEngine, UnitConverter, MemoryManager, and Constants.
 */

const ScientificEngine = require('./scientific.js');
const { UnitConverter } = require('./units.js');
const MemoryManager = require('./memory.js');
const SCIENTIFIC_CONSTANTS = require('./constants.js');

let passedCount = 0;
let totalCount = 0;

function assert(condition, message) {
  totalCount++;
  if (!condition) {
    console.error(`❌ FAIL: ${message}`);
    process.exit(1);
  } else {
    passedCount++;
    console.log(`✅ PASS: ${message}`);
  }
}

function assertClose(actual, expected, tolerance = 1e-6, message = '') {
  totalCount++;
  if (Math.abs(actual - expected) > tolerance) {
    console.error(`❌ FAIL: ${message} (Expected ~${expected}, got ${actual})`);
    process.exit(1);
  } else {
    passedCount++;
    console.log(`✅ PASS: ${message}`);
  }
}

console.log('--- 1. Testing Scientific Engine ---');
const engine = new ScientificEngine();

// Basic expressions
assert(engine.evaluate('2 + 3 * 4') === 14, 'Precedence: 2 + 3 * 4 = 14');
assert(engine.evaluate('(2 + 3) * 4') === 20, 'Parentheses: (2 + 3) * 4 = 20');
assert(engine.evaluate('10 / 2 - 3') === 2, 'Division and subtraction: 10 / 2 - 3 = 2');
assert(engine.evaluate('2 ^ 3') === 8, 'Power: 2 ^ 3 = 8');
assert(engine.evaluate('2 ^ 3 ^ 2') === 512, 'Right associativity of power: 2 ^ (3 ^ 2) = 512');

// Trig in DEG
engine.setAngleMode('DEG');
assertClose(engine.evaluate('sin(30)'), 0.5, 1e-6, 'sin(30 deg) = 0.5');
assertClose(engine.evaluate('cos(60)'), 0.5, 1e-6, 'cos(60 deg) = 0.5');
assertClose(engine.evaluate('tan(45)'), 1.0, 1e-6, 'tan(45 deg) = 1.0');

// Trig in RAD
engine.setAngleMode('RAD');
assertClose(engine.evaluate('sin(pi / 2)'), 1.0, 1e-6, 'sin(pi / 2 rad) = 1.0');
assertClose(engine.evaluate('cos(pi)'), -1.0, 1e-6, 'cos(pi rad) = -1.0');

// Functions
assert(engine.evaluate('sqrt(144)') === 12, 'sqrt(144) = 12');
assert(engine.evaluate('cbrt(27)') === 3, 'cbrt(27) = 3');
assert(engine.evaluate('abs(-42)') === 42, 'abs(-42) = 42');
assert(engine.evaluate('5!') === 120, '5! = 120');
assert(engine.evaluate('0!') === 1, '0! = 1');
assertClose(engine.evaluate('log(1000)'), 3.0, 1e-6, 'log10(1000) = 3.0');
assertClose(engine.evaluate('ln(e)'), 1.0, 1e-6, 'ln(e) = 1.0');

console.log('\n--- 2. Testing Unit Converter ---');
const converter = new UnitConverter();

// Length
assert(converter.convert(1, 'length', 'km', 'm') === 1000, '1 km = 1000 m');
assert(converter.convert(1, 'length', 'mi', 'km') > 1.609, '1 mile > 1.609 km');
assertClose(converter.convert(1, 'length', 'in', 'cm'), 2.54, 1e-4, '1 inch = 2.54 cm');

// Mass
assert(converter.convert(1, 'mass', 'kg', 'g') === 1000, '1 kg = 1000 g');
assertClose(converter.convert(1, 'mass', 'lb', 'kg'), 0.453592, 1e-4, '1 lb ~ 0.4536 kg');

// Temperature
assertClose(converter.convert(0, 'temperature', 'c', 'f'), 32, 1e-4, '0 C = 32 F');
assertClose(converter.convert(100, 'temperature', 'c', 'f'), 212, 1e-4, '100 C = 212 F');
assertClose(converter.convert(300, 'temperature', 'k', 'c'), 26.85, 1e-4, '300 K = 26.85 C');

// Digital
assert(converter.convert(1, 'digital', 'B', 'b') === 8, '1 Byte = 8 bits');
assert(converter.convert(1, 'digital', 'GB', 'MB') === 1000, '1 GB = 1000 MB');
assert(converter.convert(1, 'digital', 'GiB', 'MiB') === 1024, '1 GiB = 1024 MiB');

// Speed
assertClose(converter.convert(100, 'speed', 'km_h', 'm_s'), 27.7778, 1e-3, '100 km/h ~ 27.78 m/s');

console.log('\n--- 3. Testing Memory Manager ---');
global.localStorage = {
  store: {},
  getItem(k) { return this.store[k] || null; },
  setItem(k, v) { this.store[k] = v; }
};
const mem = new MemoryManager();
mem.clear();
assert(!mem.hasMemory(), 'Memory should be clear');
mem.store(50);
assert(mem.recall() === 50, 'Stored 50, recalled 50');
mem.add(25);
assert(mem.recall() === 75, 'M+ 25 gives 75');
mem.subtract(15);
assert(mem.recall() === 60, 'M- 15 gives 60');
mem.clear();
assert(mem.recall() === 0, 'Cleared memory is 0');

console.log('\n--- 4. Testing Constants ---');
assert(SCIENTIFIC_CONSTANTS.length >= 20, 'Constants database has >= 20 items');
const cConst = SCIENTIFIC_CONSTANTS.find(c => c.symbol === 'c');
assert(cConst && cConst.value === 299792458, 'Speed of light constant verified');

console.log(`\n==============================`);
console.log(`All ${passedCount}/${totalCount} tests passed successfully!`);
console.log(`==============================`);
