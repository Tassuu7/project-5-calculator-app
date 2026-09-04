/**
 * Build and Verification Pipeline
 * Validates module integrity, checks syntax across all engine modules,
 * and confirms asset readiness for production deployment.
 */

const fs = require('fs');
const path = require('path');

console.log('=== Starting Dark Calculator Suite Build ===');
const startTime = Date.now();

const modulesToCheck = [
  'scientific.js',
  'units.js',
  'memory.js',
  'constants.js',
  'src/server.js',
  'src/math/algebra.js',
  'src/math/calculus.js',
  'src/math/statistics.js',
  'src/math/finance.js',
  'src/math/physics.js',
  'src/math/geometry.js',
  'src/math/number-theory.js',
  'src/math/signal.js',
  'src/math/symbolic.js',
  'src/math/arbitrary-precision.js'
];

let errors = 0;
for (const mod of modulesToCheck) {
  const modPath = path.resolve(__dirname, mod);
  if (!fs.existsSync(modPath)) {
    console.error(`❌ Missing module: ${mod}`);
    errors++;
    continue;
  }
  try {
    require(modPath);
    console.log(`  ✓ Verified: ${mod}`);
  } catch (err) {
    console.error(`❌ Error validating ${mod}: ${err.message}`);
    errors++;
  }
}

if (errors > 0) {
  console.error(`\nBuild failed with ${errors} error(s).`);
  process.exit(1);
}

const elapsed = Date.now() - startTime;
console.log(`\n=== Build Completed Successfully in ${elapsed}ms ===\n`);
