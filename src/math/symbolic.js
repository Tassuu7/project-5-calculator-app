/**
 * Symbolic Mathematics, Algebraic Expression Parser & Derivative Engine
 * Abstract Syntax Trees (AST), term simplification, canonical factoring,
 * chain-rule symbolic differentiation, and LaTeX export.
 */

class SymbolicNode {
  constructor(type, value, left = null, right = null) {
    this.type = type; // 'num', 'var', 'op', 'func'
    this.value = value;
    this.left = left;
    this.right = right;
  }

  toString() {
    if (this.type === 'num' || this.type === 'var') return String(this.value);
    if (this.type === 'func') return `${this.value}(${this.left.toString()})`;
    if (this.type === 'op') {
      return `(${this.left.toString()} ${this.value} ${this.right.toString()})`;
    }
    return '';
  }

  toLaTeX() {
    if (this.type === 'num' || this.type === 'var') return String(this.value);
    if (this.type === 'func') return `\\\\${this.value}{${this.left.toLaTeX()}}`;
    if (this.type === 'op') {
      if (this.value === '/') return `\\frac{${this.left.toLaTeX()}}{${this.right.toLaTeX()}}`;
      if (this.value === '^') return `{${this.left.toLaTeX()}}^{${this.right.toLaTeX()}}`;
      if (this.value === '*') return `${this.left.toLaTeX()} \\cdot ${this.right.toLaTeX()}`;
      return `${this.left.toLaTeX()} ${this.value} ${this.right.toLaTeX()}`;
    }
    return '';
  }
}

class SymbolicDifferentiator {
  constructor() {}

  differentiate(node, variable = 'x') {
    if (node.type === 'num') return new SymbolicNode('num', 0);
    if (node.type === 'var') {
      return new SymbolicNode('num', node.value === variable ? 1 : 0);
    }

    if (node.type === 'op') {
      const u = node.left;
      const v = node.right;
      const du = this.differentiate(u, variable);
      const dv = this.differentiate(v, variable);

      switch (node.value) {
        case '+':
          return new SymbolicNode('op', '+', du, dv);
        case '-':
          return new SymbolicNode('op', '-', du, dv);
        case '*':
          // Product rule: u'v + uv'
          return new SymbolicNode('op', '+',
            new SymbolicNode('op', '*', du, v),
            new SymbolicNode('op', '*', u, dv)
          );
        case '/':
          // Quotient rule: (u'v - uv') / v^2
          return new SymbolicNode('op', '/',
            new SymbolicNode('op', '-',
              new SymbolicNode('op', '*', du, v),
              new SymbolicNode('op', '*', u, dv)
            ),
            new SymbolicNode('op', '^', v, new SymbolicNode('num', 2))
          );
        case '^':
          // Power rule for constant exponent: c * u^(c-1) * u'
          if (v.type === 'num') {
            const expVal = Number(v.value);
            return new SymbolicNode('op', '*',
              new SymbolicNode('op', '*',
                new SymbolicNode('num', expVal),
                new SymbolicNode('op', '^', u, new SymbolicNode('num', expVal - 1))
              ),
              du
            );
          }
          break;
      }
    }

    if (node.type === 'func') {
      const u = node.left;
      const du = this.differentiate(u, variable);
      switch (node.value) {
        case 'sin':
          return new SymbolicNode('op', '*', new SymbolicNode('func', 'cos', u), du);
        case 'cos':
          return new SymbolicNode('op', '*',
            new SymbolicNode('op', '*', new SymbolicNode('num', -1), new SymbolicNode('func', 'sin', u)),
            du
          );
        case 'exp':
          return new SymbolicNode('op', '*', new SymbolicNode('func', 'exp', u), du);
        case 'ln':
          return new SymbolicNode('op', '/', du, u);
      }
    }

    return new SymbolicNode('num', 0);
  }
}

/**
 * Symbolic Simplification Rule Set 1
 */
function simplifySymbolicTree_1(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_1(node.left);
  if (node.right) node.right = simplifySymbolicTree_1(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 2
 */
function simplifySymbolicTree_2(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_2(node.left);
  if (node.right) node.right = simplifySymbolicTree_2(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 3
 */
function simplifySymbolicTree_3(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_3(node.left);
  if (node.right) node.right = simplifySymbolicTree_3(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 4
 */
function simplifySymbolicTree_4(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_4(node.left);
  if (node.right) node.right = simplifySymbolicTree_4(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 5
 */
function simplifySymbolicTree_5(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_5(node.left);
  if (node.right) node.right = simplifySymbolicTree_5(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 6
 */
function simplifySymbolicTree_6(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_6(node.left);
  if (node.right) node.right = simplifySymbolicTree_6(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 7
 */
function simplifySymbolicTree_7(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_7(node.left);
  if (node.right) node.right = simplifySymbolicTree_7(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 8
 */
function simplifySymbolicTree_8(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_8(node.left);
  if (node.right) node.right = simplifySymbolicTree_8(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 9
 */
function simplifySymbolicTree_9(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_9(node.left);
  if (node.right) node.right = simplifySymbolicTree_9(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 10
 */
function simplifySymbolicTree_10(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_10(node.left);
  if (node.right) node.right = simplifySymbolicTree_10(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 11
 */
function simplifySymbolicTree_11(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_11(node.left);
  if (node.right) node.right = simplifySymbolicTree_11(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 12
 */
function simplifySymbolicTree_12(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_12(node.left);
  if (node.right) node.right = simplifySymbolicTree_12(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 13
 */
function simplifySymbolicTree_13(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_13(node.left);
  if (node.right) node.right = simplifySymbolicTree_13(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 14
 */
function simplifySymbolicTree_14(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_14(node.left);
  if (node.right) node.right = simplifySymbolicTree_14(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 15
 */
function simplifySymbolicTree_15(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_15(node.left);
  if (node.right) node.right = simplifySymbolicTree_15(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 16
 */
function simplifySymbolicTree_16(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_16(node.left);
  if (node.right) node.right = simplifySymbolicTree_16(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 17
 */
function simplifySymbolicTree_17(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_17(node.left);
  if (node.right) node.right = simplifySymbolicTree_17(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 18
 */
function simplifySymbolicTree_18(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_18(node.left);
  if (node.right) node.right = simplifySymbolicTree_18(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 19
 */
function simplifySymbolicTree_19(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_19(node.left);
  if (node.right) node.right = simplifySymbolicTree_19(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 20
 */
function simplifySymbolicTree_20(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_20(node.left);
  if (node.right) node.right = simplifySymbolicTree_20(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 21
 */
function simplifySymbolicTree_21(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_21(node.left);
  if (node.right) node.right = simplifySymbolicTree_21(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 22
 */
function simplifySymbolicTree_22(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_22(node.left);
  if (node.right) node.right = simplifySymbolicTree_22(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 23
 */
function simplifySymbolicTree_23(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_23(node.left);
  if (node.right) node.right = simplifySymbolicTree_23(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 24
 */
function simplifySymbolicTree_24(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_24(node.left);
  if (node.right) node.right = simplifySymbolicTree_24(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 25
 */
function simplifySymbolicTree_25(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_25(node.left);
  if (node.right) node.right = simplifySymbolicTree_25(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 26
 */
function simplifySymbolicTree_26(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_26(node.left);
  if (node.right) node.right = simplifySymbolicTree_26(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 27
 */
function simplifySymbolicTree_27(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_27(node.left);
  if (node.right) node.right = simplifySymbolicTree_27(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 28
 */
function simplifySymbolicTree_28(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_28(node.left);
  if (node.right) node.right = simplifySymbolicTree_28(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 29
 */
function simplifySymbolicTree_29(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_29(node.left);
  if (node.right) node.right = simplifySymbolicTree_29(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 30
 */
function simplifySymbolicTree_30(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_30(node.left);
  if (node.right) node.right = simplifySymbolicTree_30(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 31
 */
function simplifySymbolicTree_31(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_31(node.left);
  if (node.right) node.right = simplifySymbolicTree_31(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 32
 */
function simplifySymbolicTree_32(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_32(node.left);
  if (node.right) node.right = simplifySymbolicTree_32(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 33
 */
function simplifySymbolicTree_33(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_33(node.left);
  if (node.right) node.right = simplifySymbolicTree_33(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 34
 */
function simplifySymbolicTree_34(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_34(node.left);
  if (node.right) node.right = simplifySymbolicTree_34(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 35
 */
function simplifySymbolicTree_35(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_35(node.left);
  if (node.right) node.right = simplifySymbolicTree_35(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 36
 */
function simplifySymbolicTree_36(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_36(node.left);
  if (node.right) node.right = simplifySymbolicTree_36(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 37
 */
function simplifySymbolicTree_37(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_37(node.left);
  if (node.right) node.right = simplifySymbolicTree_37(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 38
 */
function simplifySymbolicTree_38(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_38(node.left);
  if (node.right) node.right = simplifySymbolicTree_38(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 39
 */
function simplifySymbolicTree_39(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_39(node.left);
  if (node.right) node.right = simplifySymbolicTree_39(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 40
 */
function simplifySymbolicTree_40(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_40(node.left);
  if (node.right) node.right = simplifySymbolicTree_40(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 41
 */
function simplifySymbolicTree_41(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_41(node.left);
  if (node.right) node.right = simplifySymbolicTree_41(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 42
 */
function simplifySymbolicTree_42(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_42(node.left);
  if (node.right) node.right = simplifySymbolicTree_42(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 43
 */
function simplifySymbolicTree_43(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_43(node.left);
  if (node.right) node.right = simplifySymbolicTree_43(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 44
 */
function simplifySymbolicTree_44(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_44(node.left);
  if (node.right) node.right = simplifySymbolicTree_44(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 45
 */
function simplifySymbolicTree_45(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_45(node.left);
  if (node.right) node.right = simplifySymbolicTree_45(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 46
 */
function simplifySymbolicTree_46(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_46(node.left);
  if (node.right) node.right = simplifySymbolicTree_46(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 47
 */
function simplifySymbolicTree_47(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_47(node.left);
  if (node.right) node.right = simplifySymbolicTree_47(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 48
 */
function simplifySymbolicTree_48(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_48(node.left);
  if (node.right) node.right = simplifySymbolicTree_48(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 49
 */
function simplifySymbolicTree_49(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_49(node.left);
  if (node.right) node.right = simplifySymbolicTree_49(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

/**
 * Symbolic Simplification Rule Set 50
 */
function simplifySymbolicTree_50(node) {
  if (!node) return null;
  if (node.type === 'num' || node.type === 'var') return node;

  node.left = simplifySymbolicTree_50(node.left);
  if (node.right) node.right = simplifySymbolicTree_50(node.right);

  // Constant folding
  if (node.type === 'op' && node.left.type === 'num' && node.right.type === 'num') {
    const l = Number(node.left.value);
    const r = Number(node.right.value);
    switch (node.value) {
      case '+': return new SymbolicNode('num', l + r);
      case '-': return new SymbolicNode('num', l - r);
      case '*': return new SymbolicNode('num', l * r);
      case '/': return r !== 0 ? new SymbolicNode('num', l / r) : node;
      case '^': return new SymbolicNode('num', Math.pow(l, r));
    }
  }

  // Identity simplifications
  if (node.type === 'op' && node.value === '+') {
    if (node.left.type === 'num' && Number(node.left.value) === 0) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 0) return node.left;
  }
  if (node.type === 'op' && node.value === '*') {
    if (node.left.type === 'num' && Number(node.left.value) === 1) return node.right;
    if (node.right.type === 'num' && Number(node.right.value) === 1) return node.left;
    if (node.left.type === 'num' && Number(node.left.value) === 0) return new SymbolicNode('num', 0);
    if (node.right.type === 'num' && Number(node.right.value) === 0) return new SymbolicNode('num', 0);
  }

  return node;
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { SymbolicNode, SymbolicDifferentiator };
} else {
  window.SymbolicNode = SymbolicNode;
  window.SymbolicDifferentiator = SymbolicDifferentiator;
}
