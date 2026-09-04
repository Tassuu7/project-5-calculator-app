/**
 * Dark Calculator & Scientific Computing Suite - Application Server
 * Serves static frontend assets and provides comprehensive mathematical REST API.
 */

const http = require('http');
const fs = require('fs');
const path = require('path');

const PORT = process.env.PORT || 3000;
const ROOT_DIR = path.resolve(__dirname, '..');

// Import core modules
let ScientificEngine, UnitConverter, SCIENTIFIC_CONSTANTS;
try {
  ScientificEngine = require('../scientific.js');
  const unitsModule = require('../units.js');
  UnitConverter = unitsModule.UnitConverter;
  SCIENTIFIC_CONSTANTS = require('../constants.js');
} catch (e) {
  console.warn('Warning: Some modules loaded with fallback:', e.message);
}

const MIME_TYPES = {
  '.html': 'text/html',
  '.css': 'text/css',
  '.js': 'application/javascript',
  '.json': 'application/json',
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
  '.txt': 'text/plain'
};

function parseJsonBody(req) {
  return new Promise((resolve, reject) => {
    let body = '';
    req.on('data', chunk => {
      body += chunk.toString();
      if (body.length > 1e6) {
        req.destroy();
        reject(new Error('Payload too large'));
      }
    });
    req.on('end', () => {
      if (!body) return resolve({});
      try {
        resolve(JSON.parse(body));
      } catch (err) {
        reject(new Error('Invalid JSON payload'));
      }
    });
  });
}

function sendJson(res, statusCode, data) {
  res.writeHead(statusCode, {
    'Content-Type': 'application/json',
    'Access-Control-Allow-Origin': '*',
    'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
    'Access-Control-Allow-Headers': 'Content-Type'
  });
  res.end(JSON.stringify(data));
}

const server = http.createServer(async (req, res) => {
  const parsedUrl = new URL(req.url, `http://localhost:${PORT}`);
  const pathname = parsedUrl.pathname;

  // Handle CORS preflight
  if (req.method === 'OPTIONS') {
    res.writeHead(204, {
      'Access-Control-Allow-Origin': '*',
      'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
      'Access-Control-Allow-Headers': 'Content-Type'
    });
    res.end();
    return;
  }

  // --- API ROUTING ---
  if (pathname === '/api/health') {
    return sendJson(res, 200, {
      status: 'healthy',
      application: 'Multi-Mode Dark Calculator & Scientific Engine',
      version: '2.0.0',
      uptimeSeconds: Math.floor(process.uptime()),
      timestamp: new Date().toISOString()
    });
  }

  if (pathname === '/api/constants' && req.method === 'GET') {
    return sendJson(res, 200, {
      status: 'success',
      count: (SCIENTIFIC_CONSTANTS || []).length,
      constants: SCIENTIFIC_CONSTANTS || []
    });
  }

  if (pathname === '/api/calculate' && req.method === 'POST') {
    try {
      const body = await parseJsonBody(req);
      const expr = body.expression || '';
      const mode = body.angleMode || 'DEG';

      if (!expr) {
        return sendJson(res, 400, { error: 'Missing expression parameter' });
      }

      const engine = new ScientificEngine();
      engine.setAngleMode(mode);
      const result = engine.evaluate(expr);

      return sendJson(res, 200, {
        status: 'success',
        expression: expr,
        angleMode: mode,
        result
      });
    } catch (err) {
      return sendJson(res, 400, { error: err.message });
    }
  }

  if (pathname === '/api/convert' && req.method === 'POST') {
    try {
      const body = await parseJsonBody(req);
      const { value, category, fromUnit, toUnit } = body;

      if (value === undefined || !category || !fromUnit || !toUnit) {
        return sendJson(res, 400, { error: 'Missing parameters: value, category, fromUnit, toUnit' });
      }

      const converter = new UnitConverter();
      const result = converter.convert(Number(value), category, fromUnit, toUnit);

      return sendJson(res, 200, {
        status: 'success',
        input: { value, category, fromUnit, toUnit },
        result
      });
    } catch (err) {
      return sendJson(res, 400, { error: err.message });
    }
  }

  // --- STATIC ASSET SERVING ---
  let safePath = pathname === '/' ? '/index.html' : pathname;
  safePath = path.normalize(safePath).replace(/^(\.\.[\/\\])+/, '');
  const filePath = path.join(ROOT_DIR, safePath);

  if (fs.existsSync(filePath) && fs.statSync(filePath).isFile()) {
    const ext = path.extname(filePath).toLowerCase();
    const contentType = MIME_TYPES[ext] || 'application/octet-stream';
    res.writeHead(200, {
      'Content-Type': contentType,
      'Cache-Control': 'no-cache, no-store, must-revalidate',
      'Pragma': 'no-cache',
      'Expires': '0'
    });
    fs.createReadStream(filePath).pipe(res);
    return;
  }

  // 404 Not Found
  sendJson(res, 404, { error: 'Not Found', path: pathname });
});

if (require.main === module) {
  server.listen(PORT, () => {
    console.log(`\n==================================================`);
    console.log(`  Dark Calculator & Scientific Engine Server v2.0.0`);
    console.log(`  Server listening on http://localhost:${PORT}`);
    console.log(`  API Health check: http://localhost:${PORT}/api/health`);
    console.log(`==================================================\n`);
  });
}

module.exports = server;
