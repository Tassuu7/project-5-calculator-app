/**
 * Scientific Math Engine
 * Implements Shunting-Yard algorithm and Reverse Polish Notation (RPN) evaluator
 * Supports operator precedence, parentheses, unary operators, trig (DEG/RAD), and functions.
 */

class ScientificEngine {
  constructor() {
    this.angleMode = 'DEG'; // 'DEG' or 'RAD'
  }

  setAngleMode(mode) {
    if (mode === 'DEG' || mode === 'RAD') {
      this.angleMode = mode;
    }
  }

  toRadians(angle) {
    return this.angleMode === 'DEG' ? (angle * Math.PI) / 180 : angle;
  }

  fromRadians(rad) {
    return this.angleMode === 'DEG' ? (rad * 180) / Math.PI : rad;
  }

  factorial(n) {
    if (n < 0 || !Number.isInteger(n)) return NaN;
    if (n === 0 || n === 1) return 1;
    if (n > 170) return Infinity; // JS overflow
    let result = 1;
    for (let i = 2; i <= n; i++) {
      result *= i;
    }
    return result;
  }

  /**
   * Tokenize mathematical expression string
   */
  tokenize(expression) {
    const tokens = [];
    let i = 0;
    const str = expression.replace(/\s+/g, '');

    const isDigit = (c) => /[0-9]/.test(c);
    const isLetter = (c) => /[a-zA-Z]/.test(c);

    while (i < str.length) {
      const char = str[i];

      // Number (including decimal)
      if (isDigit(char) || (char === '.' && (i + 1 < str.length && isDigit(str[i + 1])))) {
        let numStr = '';
        while (i < str.length && (isDigit(str[i]) || str[i] === '.')) {
          numStr += str[i];
          i++;
        }
        tokens.push({ type: 'NUMBER', value: parseFloat(numStr) });
        continue;
      }

      // Identifiers: functions & constants
      if (isLetter(char) || char === 'π') {
        let name = '';
        while (i < str.length && (isLetter(str[i]) || str[i] === 'π' || str[i] === '²' || str[i] === '⁻' || str[i] === '¹')) {
          name += str[i];
          i++;
        }

        if (name === 'π' || name.toLowerCase() === 'pi') {
          tokens.push({ type: 'NUMBER', value: Math.PI });
        } else if (name === 'e') {
          tokens.push({ type: 'NUMBER', value: Math.E });
        } else {
          tokens.push({ type: 'FUNCTION', value: name.toLowerCase() });
        }
        continue;
      }

      // Operators and parentheses
      if ('+-*/×÷^%()!'.includes(char)) {
        // Standardize operators
        let op = char;
        if (op === '×') op = '*';
        if (op === '÷') op = '/';

        // Check for unary minus vs binary minus
        if (op === '-') {
          const prevToken = tokens[tokens.length - 1];
          if (!prevToken || prevToken.type === 'OPERATOR' || (prevToken.type === 'PAREN' && prevToken.value === '(')) {
            tokens.push({ type: 'UNARY_MINUS', value: 'u-' });
            i++;
            continue;
          }
        }

        if (op === '(' || op === ')') {
          tokens.push({ type: 'PAREN', value: op });
        } else if (op === '!') {
          tokens.push({ type: 'POSTFIX', value: '!' });
        } else {
          tokens.push({ type: 'OPERATOR', value: op });
        }
        i++;
        continue;
      }

      i++;
    }

    return tokens;
  }

  /**
   * Shunting-Yard Algorithm to convert Infix to RPN (Postfix)
   */
  toRPN(tokens) {
    const outputQueue = [];
    const operatorStack = [];

    const precedence = {
      '+': 1,
      '-': 1,
      '*': 2,
      '/': 2,
      '%': 2,
      '^': 3,
      'u-': 4,
      '!': 5
    };

    const isRightAssociative = (op) => op === '^' || op === 'u-';

    for (let i = 0; i < tokens.length; i++) {
      const token = tokens[i];

      if (token.type === 'NUMBER') {
        outputQueue.push(token);
      } else if (token.type === 'FUNCTION') {
        operatorStack.push(token);
      } else if (token.type === 'POSTFIX') {
        outputQueue.push(token);
      } else if (token.type === 'UNARY_MINUS' || token.type === 'OPERATOR') {
        const op1 = token.value;
        while (operatorStack.length > 0) {
          const top = operatorStack[operatorStack.length - 1];
          if (top.type === 'FUNCTION') {
            outputQueue.push(operatorStack.pop());
          } else if (top.type === 'OPERATOR' || top.type === 'UNARY_MINUS') {
            const op2 = top.value;
            const p1 = precedence[op1] || 0;
            const p2 = precedence[op2] || 0;

            if ((!isRightAssociative(op1) && p1 <= p2) || (isRightAssociative(op1) && p1 < p2)) {
              outputQueue.push(operatorStack.pop());
            } else {
              break;
            }
          } else {
            break;
          }
        }
        operatorStack.push(token);
      } else if (token.type === 'PAREN' && token.value === '(') {
        operatorStack.push(token);
      } else if (token.type === 'PAREN' && token.value === ')') {
        while (operatorStack.length > 0 && !(operatorStack[operatorStack.length - 1].type === 'PAREN' && operatorStack[operatorStack.length - 1].value === '(')) {
          outputQueue.push(operatorStack.pop());
        }

        if (operatorStack.length === 0) {
          throw new Error('Mismatched parentheses');
        }
        // Discard the left paren '('
        operatorStack.pop();

        // If the token at the top of the stack is a function, pop it onto the output queue
        if (operatorStack.length > 0 && operatorStack[operatorStack.length - 1].type === 'FUNCTION') {
          outputQueue.push(operatorStack.pop());
        }
      }
    }

    while (operatorStack.length > 0) {
      const top = operatorStack.pop();
      if (top.type === 'PAREN') {
        throw new Error('Mismatched parentheses');
      }
      outputQueue.push(top);
    }

    return outputQueue;
  }

  /**
   * Evaluate RPN queue
   */
  evaluateRPN(rpnQueue) {
    const stack = [];

    for (const token of rpnQueue) {
      if (token.type === 'NUMBER') {
        stack.push(token.value);
      } else if (token.type === 'UNARY_MINUS') {
        if (stack.length < 1) throw new Error('Invalid syntax');
        const val = stack.pop();
        stack.push(-val);
      } else if (token.type === 'POSTFIX' && token.value === '!') {
        if (stack.length < 1) throw new Error('Invalid syntax');
        const val = stack.pop();
        stack.push(this.factorial(val));
      } else if (token.type === 'OPERATOR') {
        if (stack.length < 2) throw new Error('Invalid syntax');
        const b = stack.pop();
        const a = stack.pop();

        switch (token.value) {
          case '+':
            stack.push(a + b);
            break;
          case '-':
            stack.push(a - b);
            break;
          case '*':
            stack.push(a * b);
            break;
          case '/':
            if (b === 0) throw new Error('Cannot divide by 0');
            stack.push(a / b);
            break;
          case '%':
            stack.push(a % b);
            break;
          case '^':
            stack.push(Math.pow(a, b));
            break;
          default:
            throw new Error(`Unknown operator ${token.value}`);
        }
      } else if (token.type === 'FUNCTION') {
        if (stack.length < 1) throw new Error('Invalid syntax');
        const arg = stack.pop();
        const fn = token.value;

        switch (fn) {
          case 'sin':
            stack.push(Math.sin(this.toRadians(arg)));
            break;
          case 'cos':
            stack.push(Math.cos(this.toRadians(arg)));
            break;
          case 'tan': {
            // Check for asymptote tan(90 deg)
            const rad = this.toRadians(arg);
            const cosVal = Math.cos(rad);
            if (Math.abs(cosVal) < 1e-15) throw new Error('Undefined (asymptote)');
            stack.push(Math.tan(rad));
            break;
          }
          case 'asin':
          case 'arcsin':
            if (arg < -1 || arg > 1) throw new Error('Domain error for asin');
            stack.push(this.fromRadians(Math.asin(arg)));
            break;
          case 'acos':
          case 'arccos':
            if (arg < -1 || arg > 1) throw new Error('Domain error for acos');
            stack.push(this.fromRadians(Math.acos(arg)));
            break;
          case 'atan':
          case 'arctan':
            stack.push(this.fromRadians(Math.atan(arg)));
            break;
          case 'log':
          case 'log10':
            if (arg <= 0) throw new Error('Domain error (log <= 0)');
            stack.push(Math.log10(arg));
            break;
          case 'ln':
            if (arg <= 0) throw new Error('Domain error (ln <= 0)');
            stack.push(Math.log(arg));
            break;
          case 'sqrt':
          case '√':
            if (arg < 0) throw new Error('Domain error (negative square root)');
            stack.push(Math.sqrt(arg));
            break;
          case 'cbrt':
            stack.push(Math.cbrt(arg));
            break;
          case 'abs':
            stack.push(Math.abs(arg));
            break;
          case 'sqr':
            stack.push(arg * arg);
            break;
          case 'cube':
            stack.push(arg * arg * arg);
            break;
          case 'inv':
          case 'recip':
            if (arg === 0) throw new Error('Cannot divide by 0');
            stack.push(1 / arg);
            break;
          case 'exp':
            stack.push(Math.exp(arg));
            break;
          default:
            throw new Error(`Unknown function: ${fn}`);
        }
      }
    }

    if (stack.length !== 1) {
      throw new Error('Invalid expression');
    }

    const result = stack[0];
    if (isNaN(result)) throw new Error('Result is NaN');
    if (!isFinite(result)) throw new Error('Infinity / Overflow');

    // Precision rounding
    return parseFloat(result.toPrecision(12));
  }

  evaluate(expression) {
    const tokens = this.tokenize(expression);
    if (tokens.length === 0) return 0;
    const rpn = this.toRPN(tokens);
    return this.evaluateRPN(rpn);
  }
}

// Export for node testing or window browser global
if (typeof module !== 'undefined' && module.exports) {
  module.exports = ScientificEngine;
} else {
  window.ScientificEngine = ScientificEngine;
}
