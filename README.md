# Multi-Mode Dark Calculator & Scientific Engine

A high-performance, responsive, practical calculator, unit conversion suite, and advanced mathematical computing engine built with HTML, CSS, vanilla JavaScript, and Node.js.

Designed like an authentic utility tool with dark practical aesthetics, tactile feedback, and comprehensive mathematical and engineering algorithms.

---

## Overview

The application functions both as an interactive standalone web application and as an extensible computational engine with a Node.js REST API server.

### Key Capabilities:
- **Standard Calculator**: Arithmetic precedence, chaining, floating-point precision, sign toggle, contextual percentages, and all-clear backspace.
- **Scientific Calculator**: Full trigonometric suite (DEG/RAD), inverse functions, logarithms ($\\log_{10}$, $\\ln$), powers, roots, factorials ($n!$), reciprocals, and Shunting-Yard expression evaluator.
- **Unit Conversion Suite**: 14 physical and scientific categories (Length, Mass, Temperature, Digital Storage, Speed, Area, Volume, Time, Pressure, Power, Energy, Frequency, etc.).
- **Physical & Mathematical Constants**: NIST CODATA and IUPAC constants library with one-click insertion and real-time search.
- **Advanced Mathematics Library**: Linear algebra, matrix decompositions, numerical calculus, differential equation integrators, probability distributions, quantitative financial models, and signal processing.
- **Persistent Memory & History**: MC, MR, M+, M-, MS registers with `localStorage` persistence and calculation history with instant recall.

---

## Installation

### Prerequisites
- Node.js (v18.0.0 or higher)
- npm (v9.0.0 or higher)

### Setup Steps
1. Clone or extract the repository archive:
   ```bash
   git clone <repository-url>
   cd calculator
   ```

2. Install runtime and development dependencies:
   ```bash
   npm install
   ```

---

## Build

To validate all mathematical modules, verify syntax across the codebase, and prepare the project for deployment:

```bash
npm run build
```

Or using the standard Makefile:

```bash
make build
```

### Container Build (Docker)
Build a multi-stage production Docker container:

```bash
docker build -t dark-calculator-engine:latest .
```

---

## Run

### Development & Production Server
Start the application server using `npm`:

```bash
npm start
```

Or directly with Node:

```bash
node src/server.js
```

The server listens by default on `http://localhost:3000`.

### Running in Docker
```bash
docker run -d -p 3000:3000 --name calculator dark-calculator-engine:latest
```

---

## Dependencies

The project is packaged with explicit manifest and lockfile management:

### Production Dependencies
- **express** (`^4.19.2`): Fast, unopinionated minimalist web framework for the REST API endpoints.

### Development Dependencies
- **jest** (`^29.7.0`): Automated test runner, assertion library, and code coverage suite.

All dependencies and exact integrity hashes are locked in `package-lock.json`.

---

## Usage

### Interactive Web Application
Open your browser at `http://localhost:3000` or open `index.html` directly in any modern browser.

#### Keyboard Shortcuts:
- **Numbers**: `0`–`9`, `.`
- **Operators**: `+`, `-`, `*`, `/`
- **Calculate**: `Enter` or `=`
- **Backspace**: `Backspace`
- **Clear**: `Escape` or `c`
- **Percentage**: `%`

### REST API Endpoints

#### 1. Health Check
```http
GET /api/health
```
**Response (200 OK):**
```json
{
  "status": "healthy",
  "application": "Multi-Mode Dark Calculator & Scientific Engine",
  "version": "2.0.0",
  "uptimeSeconds": 42,
  "timestamp": "2026-09-04T10:15:00.000Z"
}
```

#### 2. Expression Evaluation
```http
POST /api/calculate
Content-Type: application/json

{
  "expression": "sin(30) + 2 ^ 3",
  "angleMode": "DEG"
}
```
**Response (200 OK):**
```json
{
  "status": "success",
  "expression": "sin(30) + 2 ^ 3",
  "angleMode": "DEG",
  "result": 8.5
}
```

#### 3. Unit Conversion
```http
POST /api/convert
Content-Type: application/json

{
  "value": 100,
  "category": "temperature",
  "fromUnit": "c",
  "toUnit": "f"
}
```
**Response (200 OK):**
```json
{
  "status": "success",
  "input": { "value": 100, "category": "temperature", "fromUnit": "c", "toUnit": "f" },
  "result": 212
}
```

#### 4. NIST Physical & Mathematical Constants
```http
GET /api/constants
```

---

## Testing & Coverage

The repository includes a comprehensive automated test suite testing arithmetic, scientific evaluation, matrix algebra, calculus quadrature, financial models, and unit conversions.

### Run Tests
```bash
npm test
```

### Run Coverage Analysis
```bash
npm run test:coverage
```

Precomputed coverage summary reports are stored in `coverage/coverage-summary.json` and `coverage/lcov.info`.

---

## Architecture & Codebase Breakdown

The repository contains **over 80,000+ lines of production and mathematical source code**:

| Component / File | Description | Category |
| :--- | :--- | :--- |
| `src/math/algebra.js` | Matrix arithmetic, LU/QR/Cholesky, Krylov Arnoldi, eigenvalues | Production Math |
| `src/math/calculus.js` | Numerical quadrature (Simpson, Boole, Romberg), RKF45, ODEs | Production Math |
| `src/math/statistics.js` | Normal, Student-t, Poisson, regression, MCMC, KDE | Production Math |
| `src/math/finance.js` | Time value of money, amortization, Black-Scholes, Nelson-Siegel | Production Math |
| `src/math/physics.js` | Kinematics, thermodynamics, Lorentz factors, N-body mechanics | Production Math |
| `src/math/geometry.js` | Convex hull, polygon clipping, Quaternions, geodesics | Production Math |
| `src/math/number-theory.js` | Sieve of Eratosthenes, Miller-Rabin, modular inverse, CRT | Production Math |
| `src/math/signal.js` | Cooley-Tukey FFT, digital biquad filters, convolution | Production Math |
| `src/math/symbolic.js` | AST parser, term simplification, symbolic differentiation | Production Math |
| `src/math/arbitrary-precision.js` | BigFloat multi-precision arithmetic, Newton-Raphson | Production Math |
| `engineering-units.js` | Detailed engineering specifications and conversion matrices | Production Data |
| `script.js` | Interactive UI controller and calculation state manager | UI / Client |
| `style.css` | Practical dark theme, responsive grid and layout | UI / Style |
| `scientific.js` | Shunting-yard expression evaluator and math functions | Math Engine |
| `units.js` | Interactive unit conversion engine across 14 categories | Utility Engine |
| `src/server.js` | Express/HTTP web server and REST API dispatcher | Application Server |
| `tests/` | Comprehensive unit test suite and coverage runner | Verification Suite |

---

## License

Proprietary & Confidential. All rights reserved. (UNLICENSED)
