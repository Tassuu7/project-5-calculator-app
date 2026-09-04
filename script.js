/**
 * Main Calculator Controller
 * Manages Standard Mode, Scientific Mode, Unit Converter, Memory, History, and Constants.
 */

document.addEventListener('DOMContentLoaded', () => {
  // DOM Elements - Shell & Modes
  const appContainer = document.querySelector('.calculator-app');
  const tabButtons = document.querySelectorAll('.tab-btn');
  const calcView = document.getElementById('calc-view');
  const converterView = document.getElementById('converter-view');
  const scientificGrid = document.getElementById('scientific-grid');

  // DOM Elements - Display & Keypad
  const previousOperandElement = document.getElementById('previous-operand');
  const currentOperandElement = document.getElementById('current-operand');
  const memoryIndicator = document.getElementById('memory-indicator');
  const angleIndicator = document.getElementById('angle-indicator');
  const angleToggleBtn = document.getElementById('angle-toggle-btn');

  // DOM Elements - Modals
  const historyToggleBtn = document.getElementById('history-toggle-btn');
  const historyPanel = document.getElementById('history-panel');
  const closeHistoryBtn = document.getElementById('close-history-btn');
  const clearHistoryBtn = document.getElementById('clear-history-btn');
  const exportHistoryBtn = document.getElementById('export-history-btn');
  const historyList = document.getElementById('history-list');

  const constantsToggleBtn = document.getElementById('constants-toggle-btn');
  const constantsPanel = document.getElementById('constants-panel');
  const closeConstantsBtn = document.getElementById('close-constants-btn');
  const constantsList = document.getElementById('constants-list');
  const constantsSearch = document.getElementById('constants-search');

  // DOM Elements - Converter
  const converterCategory = document.getElementById('converter-category');
  const converterFromValue = document.getElementById('converter-from-value');
  const converterFromUnit = document.getElementById('converter-from-unit');
  const converterToValue = document.getElementById('converter-to-value');
  const converterToUnit = document.getElementById('converter-to-unit');
  const converterSwapBtn = document.getElementById('converter-swap-btn');
  const converterFormula = document.getElementById('converter-formula');

  // Engines
  const sciEngine = new ScientificEngine();
  const unitConverter = new UnitConverter();
  let currentMode = 'standard'; // 'standard', 'scientific', 'converter'

  // Calculator State
  let currentOperand = '0';
  let previousOperand = '';
  let operation = undefined;
  let shouldResetScreen = false;
  let hasError = false;
  let scientificExpression = ''; // used for continuous scientific expressions with parentheses

  // Memory Module
  const memory = new MemoryManager((hasMem, val) => {
    if (memoryIndicator) {
      memoryIndicator.hidden = !hasMem;
      memoryIndicator.title = `Memory: ${val}`;
    }
  });

  // History State
  let history = [];
  const maxHistoryItems = 30;

  function loadHistory() {
    try {
      const saved = localStorage.getItem('calc_full_history');
      if (saved) {
        history = JSON.parse(saved);
        renderHistory();
      }
    } catch (e) {
      history = [];
    }
  }

  function saveHistory() {
    try {
      localStorage.setItem('calc_full_history', JSON.stringify(history));
    } catch (e) {}
  }

  function addHistory(expression, result) {
    const timeStr = new Date().toLocaleTimeString([], { hour: '2-digit', minute: '2-digit', second: '2-digit' });
    history.unshift({ expression, result, time: timeStr });
    if (history.length > maxHistoryItems) history.pop();
    saveHistory();
    renderHistory();
  }

  function clearHistory() {
    history = [];
    saveHistory();
    renderHistory();
  }

  function renderHistory() {
    if (!historyList) return;
    if (history.length === 0) {
      historyList.innerHTML = '<li class="history-empty">No calculations yet</li>';
      return;
    }

    historyList.innerHTML = '';
    history.forEach(item => {
      const li = document.createElement('li');
      li.className = 'history-item';
      li.tabIndex = 0;
      li.setAttribute('role', 'button');
      li.setAttribute('aria-label', `Use result ${item.result}`);
      li.innerHTML = `
        <span class="history-time">${escapeHtml(item.time)}</span>
        <span class="history-expression">${escapeHtml(item.expression)} =</span>
        <span class="history-result">${escapeHtml(item.result)}</span>
      `;

      const selectItem = () => {
        currentOperand = item.result.toString();
        previousOperand = '';
        operation = undefined;
        shouldResetScreen = true;
        hasError = false;
        updateDisplay();
        historyPanel.hidden = true;
        historyToggleBtn.classList.remove('active');
      };

      li.addEventListener('click', selectItem);
      li.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') {
          e.preventDefault();
          selectItem();
        }
      });

      historyList.appendChild(li);
    });
  }

  function exportHistory() {
    if (history.length === 0) {
      alert('History is empty.');
      return;
    }
    const textLines = history.map(h => `[${h.time}] ${h.expression} = ${h.result}`);
    const blob = new Blob([textLines.join('\n')], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `calculator_history_${Date.now()}.txt`;
    a.click();
    URL.revokeObjectURL(url);
  }

  // Display Updates
  function formatDisplayNumber(numberStr) {
    if (numberStr === '' || numberStr === undefined || numberStr === null) return '';
    if (numberStr === 'Cannot divide by 0' || numberStr === 'Error' || String(numberStr).includes('Error')) {
      return numberStr;
    }
    if (String(numberStr).includes('e')) return numberStr;

    const parts = numberStr.toString().split('.');
    const integerPart = parseFloat(parts[0]);
    const decimalPart = parts[1];

    let integerDisplay;
    if (isNaN(integerPart)) {
      integerDisplay = '';
    } else {
      integerDisplay = integerPart.toLocaleString('en', { maximumFractionDigits: 0 });
      if (parts[0].startsWith('-') && !integerDisplay.startsWith('-')) {
        integerDisplay = '-' + integerDisplay;
      }
    }

    if (decimalPart != null) {
      return `${integerDisplay}.${decimalPart}`;
    }
    return integerDisplay;
  }

  function updateDisplay() {
    if (hasError) {
      currentOperandElement.textContent = currentOperand;
      previousOperandElement.textContent = '';
      return;
    }

    currentOperandElement.textContent = formatDisplayNumber(currentOperand) || '0';

    if (operation != null && previousOperand !== '') {
      previousOperandElement.textContent = `${formatDisplayNumber(previousOperand)} ${operation}`;
    } else {
      previousOperandElement.textContent = '';
    }

    // Highlight active operator
    document.querySelectorAll('.btn-operator').forEach(btn => {
      if (operation && btn.dataset.operator === operation && currentOperand === '') {
        btn.classList.add('active-operator');
      } else {
        btn.classList.remove('active-operator');
      }
    });
  }

  function clearAll() {
    currentOperand = '0';
    previousOperand = '';
    operation = undefined;
    shouldResetScreen = false;
    hasError = false;
    scientificExpression = '';
    updateDisplay();
  }

  function deleteLast() {
    if (hasError) {
      clearAll();
      return;
    }
    if (shouldResetScreen) {
      currentOperand = '0';
      shouldResetScreen = false;
      updateDisplay();
      return;
    }
    if (currentOperand === '0') return;

    if (currentOperand.length === 1 || (currentOperand.length === 2 && currentOperand.startsWith('-'))) {
      currentOperand = '0';
    } else {
      currentOperand = currentOperand.slice(0, -1);
    }
    updateDisplay();
  }

  function appendNumber(number) {
    if (hasError) {
      clearAll();
    }
    if (shouldResetScreen) {
      currentOperand = '';
      shouldResetScreen = false;
    }
    if (number === '.' && currentOperand.includes('.')) return;

    if (currentOperand === '0' && number !== '.') {
      currentOperand = number.toString();
    } else if (currentOperand === '' && number === '.') {
      currentOperand = '0.';
    } else {
      if (currentOperand.replace(/[^0-9]/g, '').length >= 15) return;
      currentOperand += number.toString();
    }
    updateDisplay();
  }

  function negate() {
    if (hasError || currentOperand === '0') return;
    if (currentOperand.startsWith('-')) {
      currentOperand = currentOperand.slice(1);
    } else {
      currentOperand = '-' + currentOperand;
    }
    updateDisplay();
  }

  function applyPercent() {
    if (hasError) return;
    const current = parseFloat(currentOperand);
    if (isNaN(current)) return;

    if (previousOperand !== '' && operation) {
      const prev = parseFloat(previousOperand);
      if (operation === '+' || operation === '-') {
        const percentVal = (prev * current) / 100;
        currentOperand = formatRawNumber(percentVal);
      } else if (operation === '×' || operation === '÷') {
        const percentVal = current / 100;
        currentOperand = formatRawNumber(percentVal);
      }
    } else {
      currentOperand = formatRawNumber(current / 100);
    }
    updateDisplay();
  }

  function chooseOperation(op) {
    if (hasError) clearAll();

    if (currentOperand === '' && previousOperand !== '') {
      operation = op;
      updateDisplay();
      return;
    }

    if (previousOperand !== '') {
      compute();
    }

    if (hasError) return;

    operation = op;
    previousOperand = currentOperand;
    currentOperand = '';
    shouldResetScreen = false;
    updateDisplay();
  }

  function compute() {
    if (operation === undefined || previousOperand === '' || currentOperand === '') {
      return;
    }

    const prev = parseFloat(previousOperand);
    const curr = parseFloat(currentOperand);

    if (isNaN(prev) || isNaN(curr)) return;

    let res;
    switch (operation) {
      case '+':
        res = prev + curr;
        break;
      case '-':
        res = prev - curr;
        break;
      case '×':
        res = prev * curr;
        break;
      case '÷':
        if (curr === 0) {
          setError('Cannot divide by 0');
          return;
        }
        res = prev / curr;
        break;
      default:
        return;
    }

    const formatted = formatResult(res);
    addHistory(`${formatDisplayNumber(previousOperand)} ${operation} ${formatDisplayNumber(currentOperand)}`, formatted);

    currentOperand = formatted.toString();
    previousOperand = '';
    operation = undefined;
    shouldResetScreen = true;
    updateDisplay();
  }

  function formatResult(num) {
    if (!isFinite(num)) {
      setError('Error');
      return 'Error';
    }
    const rounded = parseFloat(num.toPrecision(12));
    if (Math.abs(rounded) >= 1e12 || (Math.abs(rounded) > 0 && Math.abs(rounded) < 1e-7)) {
      return rounded.toExponential(6).replace(/\.?0+e/, 'e');
    }
    return rounded.toString();
  }

  function formatRawNumber(num) {
    if (!isFinite(num)) return '0';
    return parseFloat(num.toFixed(10)).toString();
  }

  function setError(msg) {
    hasError = true;
    currentOperand = msg;
    previousOperand = '';
    operation = undefined;
    shouldResetScreen = true;
    updateDisplay();
  }

  // Scientific Single-Operand Functions
  function applyScientificFunction(fn) {
    if (hasError) clearAll();
    const curr = parseFloat(currentOperand);
    if (isNaN(curr)) return;

    try {
      let result;
      let expr = `${fn}(${currentOperand})`;

      switch (fn) {
        case 'sin':
          result = sciEngine.toRadians(curr);
          result = Math.sin(result);
          break;
        case 'cos':
          result = sciEngine.toRadians(curr);
          result = Math.cos(result);
          break;
        case 'tan': {
          const rad = sciEngine.toRadians(curr);
          if (Math.abs(Math.cos(rad)) < 1e-15) throw new Error('Undefined');
          result = Math.tan(rad);
          break;
        }
        case 'asin':
          if (curr < -1 || curr > 1) throw new Error('Domain error');
          result = sciEngine.fromRadians(Math.asin(curr));
          break;
        case 'acos':
          if (curr < -1 || curr > 1) throw new Error('Domain error');
          result = sciEngine.fromRadians(Math.acos(curr));
          break;
        case 'atan':
          result = sciEngine.fromRadians(Math.atan(curr));
          break;
        case 'ln':
          if (curr <= 0) throw new Error('Domain error');
          result = Math.log(curr);
          break;
        case 'log':
          if (curr <= 0) throw new Error('Domain error');
          result = Math.log10(curr);
          break;
        case 'sqrt':
          if (curr < 0) throw new Error('Domain error');
          result = Math.sqrt(curr);
          expr = `√(${currentOperand})`;
          break;
        case 'sqr':
          result = curr * curr;
          expr = `(${currentOperand})²`;
          break;
        case 'fact':
          if (curr < 0 || !Number.isInteger(curr) || curr > 170) throw new Error('Invalid factorial');
          result = sciEngine.factorial(curr);
          expr = `${currentOperand}!`;
          break;
        case 'inv':
          if (curr === 0) throw new Error('Cannot divide by 0');
          result = 1 / curr;
          expr = `1/(${currentOperand})`;
          break;
        case 'abs':
          result = Math.abs(curr);
          expr = `|${currentOperand}|`;
          break;
        case 'exp':
          result = Math.exp(curr);
          expr = `e^(${currentOperand})`;
          break;
        case 'pi':
          currentOperand = Math.PI.toString();
          shouldResetScreen = true;
          updateDisplay();
          return;
        case 'e':
          currentOperand = Math.E.toString();
          shouldResetScreen = true;
          updateDisplay();
          return;
        case 'pow':
          // x^y sets up binary exponent
          previousOperand = currentOperand;
          operation = '^';
          currentOperand = '';
          shouldResetScreen = false;
          updateDisplay();
          return;
        case '(':
        case ')':
          // Append parenthesis into current input if needed
          return;
        default:
          return;
      }

      const formatted = formatResult(result);
      addHistory(expr, formatted);
      currentOperand = formatted.toString();
      shouldResetScreen = true;
      updateDisplay();
    } catch (err) {
      setError(err.message || 'Error');
    }
  }

  // Mode Switcher Handler
  function switchMode(mode) {
    currentMode = mode;
    tabButtons.forEach(btn => {
      const isActive = btn.dataset.mode === mode;
      btn.classList.toggle('active', isActive);
      btn.setAttribute('aria-selected', isActive ? 'true' : 'false');
    });

    if (mode === 'standard') {
      calcView.hidden = false;
      converterView.hidden = true;
      scientificGrid.hidden = true;
      appContainer.classList.remove('mode-scientific');
    } else if (mode === 'scientific') {
      calcView.hidden = false;
      converterView.hidden = true;
      scientificGrid.hidden = false;
      appContainer.classList.add('mode-scientific');
    } else if (mode === 'converter') {
      calcView.hidden = true;
      converterView.hidden = false;
      appContainer.classList.remove('mode-scientific');
      populateConverterCategories();
      recalculateConverter();
    }
  }

  // Converter Handlers
  function populateConverterCategories() {
    if (!converterCategory || converterCategory.options.length > 0) return;
    const cats = unitConverter.getCategories();
    cats.forEach(cat => {
      const opt = document.createElement('option');
      opt.value = cat.id;
      opt.textContent = cat.name;
      converterCategory.appendChild(opt);
    });

    updateConverterUnits();
  }

  function updateConverterUnits() {
    const catId = converterCategory.value;
    const units = unitConverter.getUnitsForCategory(catId);

    converterFromUnit.innerHTML = '';
    converterToUnit.innerHTML = '';

    units.forEach((u, i) => {
      const optFrom = document.createElement('option');
      optFrom.value = u.id;
      optFrom.textContent = u.name;
      converterFromUnit.appendChild(optFrom);

      const optTo = document.createElement('option');
      optTo.value = u.id;
      optTo.textContent = u.name;
      converterToUnit.appendChild(optTo);
    });

    // Select second unit for "To" if available
    if (units.length > 1) {
      converterToUnit.selectedIndex = 1;
    }

    recalculateConverter();
  }

  function recalculateConverter() {
    if (!converterFromValue || !converterCategory) return;
    const catId = converterCategory.value;
    const fromVal = parseFloat(converterFromValue.value);
    const fromUnit = converterFromUnit.value;
    const toUnit = converterToUnit.value;

    if (isNaN(fromVal)) {
      converterToValue.value = '—';
      converterFormula.textContent = '';
      return;
    }

    const converted = unitConverter.convert(fromVal, catId, fromUnit, toUnit);
    if (isNaN(converted)) {
      converterToValue.value = 'Error';
      converterFormula.textContent = '';
      return;
    }

    converterToValue.value = converted;

    // Show conversion rate
    const unit1 = unitConverter.convert(1, catId, fromUnit, toUnit);
    if (!isNaN(unit1)) {
      converterFormula.textContent = `1 ${fromUnit} = ${unit1} ${toUnit}`;
    } else {
      converterFormula.textContent = '';
    }
  }

  function swapConverterUnits() {
    const fromIdx = converterFromUnit.selectedIndex;
    const toIdx = converterToUnit.selectedIndex;
    converterFromUnit.selectedIndex = toIdx;
    converterToUnit.selectedIndex = fromIdx;
    recalculateConverter();
  }

  // Constants Modal
  function populateConstants(filter = '') {
    if (!constantsList || typeof SCIENTIFIC_CONSTANTS === 'undefined') return;
    constantsList.innerHTML = '';
    const q = filter.toLowerCase().trim();

    const filtered = SCIENTIFIC_CONSTANTS.filter(c => {
      return c.name.toLowerCase().includes(q) || c.symbol.toLowerCase().includes(q);
    });

    if (filtered.length === 0) {
      constantsList.innerHTML = '<li class="history-empty">No constants found</li>';
      return;
    }

    filtered.forEach(c => {
      const li = document.createElement('li');
      li.className = 'constant-item';
      li.tabIndex = 0;
      li.innerHTML = `
        <div class="constant-info">
          <span class="constant-symbol">${escapeHtml(c.symbol)}</span>
          <span class="constant-name">${escapeHtml(c.name)} ${c.unit ? `(${escapeHtml(c.unit)})` : ''}</span>
        </div>
        <span class="constant-val">${c.value}</span>
      `;

      li.addEventListener('click', () => {
        currentOperand = c.value.toString();
        shouldResetScreen = true;
        updateDisplay();
        constantsPanel.hidden = true;
        constantsToggleBtn.classList.remove('active');
      });

      constantsList.appendChild(li);
    });
  }

  // Angle Toggle (DEG / RAD)
  function toggleAngleMode() {
    const newMode = sciEngine.angleMode === 'DEG' ? 'RAD' : 'DEG';
    sciEngine.setAngleMode(newMode);
    angleIndicator.textContent = newMode;
    if (angleToggleBtn) angleToggleBtn.textContent = newMode;
  }

  // Event Listeners - Modes
  tabButtons.forEach(btn => {
    btn.addEventListener('click', () => {
      switchMode(btn.dataset.mode);
    });
  });

  // Event Listeners - Number & Operators
  document.querySelectorAll('[data-number]').forEach(button => {
    button.addEventListener('click', () => {
      appendNumber(button.dataset.number);
    });
  });

  document.querySelectorAll('[data-operator]').forEach(button => {
    button.addEventListener('click', () => {
      chooseOperation(button.dataset.operator);
    });
  });

  document.querySelectorAll('[data-action]').forEach(button => {
    button.addEventListener('click', () => {
      const action = button.dataset.action;
      switch (action) {
        case 'clear':
          clearAll();
          break;
        case 'backspace':
          deleteLast();
          break;
        case 'percent':
          applyPercent();
          break;
        case 'negate':
          negate();
          break;
        case 'equals':
          compute();
          break;
      }
    });
  });

  // Event Listeners - Scientific
  document.querySelectorAll('[data-sci]').forEach(button => {
    button.addEventListener('click', () => {
      applyScientificFunction(button.dataset.sci);
    });
  });

  if (angleToggleBtn) {
    angleToggleBtn.addEventListener('click', toggleAngleMode);
  }

  // Event Listeners - Memory
  document.querySelectorAll('[data-memory]').forEach(button => {
    button.addEventListener('click', () => {
      const memOp = button.dataset.memory;
      switch (memOp) {
        case 'mc':
          memory.clear();
          break;
        case 'mr':
          currentOperand = memory.recall().toString();
          shouldResetScreen = true;
          updateDisplay();
          break;
        case 'm+':
          memory.add(currentOperand);
          shouldResetScreen = true;
          break;
        case 'm-':
          memory.subtract(currentOperand);
          shouldResetScreen = true;
          break;
        case 'ms':
          memory.store(currentOperand);
          shouldResetScreen = true;
          break;
      }
    });
  });

  // Event Listeners - Converter
  if (converterCategory) {
    converterCategory.addEventListener('change', updateConverterUnits);
    converterFromUnit.addEventListener('change', recalculateConverter);
    converterToUnit.addEventListener('change', recalculateConverter);
    converterFromValue.addEventListener('input', recalculateConverter);
    converterSwapBtn.addEventListener('click', swapConverterUnits);
  }

  // Event Listeners - Panels
  historyToggleBtn.addEventListener('click', () => {
    constantsPanel.hidden = true;
    constantsToggleBtn.classList.remove('active');
    const isHidden = historyPanel.hidden;
    historyPanel.hidden = !isHidden;
    historyToggleBtn.classList.toggle('active', isHidden);
  });

  closeHistoryBtn.addEventListener('click', () => {
    historyPanel.hidden = true;
    historyToggleBtn.classList.remove('active');
  });

  clearHistoryBtn.addEventListener('click', clearHistory);
  exportHistoryBtn.addEventListener('click', exportHistory);

  constantsToggleBtn.addEventListener('click', () => {
    historyPanel.hidden = true;
    historyToggleBtn.classList.remove('active');
    const isHidden = constantsPanel.hidden;
    constantsPanel.hidden = !isHidden;
    constantsToggleBtn.classList.toggle('active', isHidden);
    if (isHidden) {
      constantsSearch.value = '';
      populateConstants();
    }
  });

  closeConstantsBtn.addEventListener('click', () => {
    constantsPanel.hidden = true;
    constantsToggleBtn.classList.remove('active');
  });

  constantsSearch.addEventListener('input', (e) => {
    populateConstants(e.target.value);
  });

  // Keyboard Support
  window.addEventListener('keydown', (e) => {
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA' || e.target.tagName === 'SELECT') {
      return;
    }

    let targetButton = null;

    if ((e.key >= '0' && e.key <= '9') || e.key === '.') {
      appendNumber(e.key);
      targetButton = document.querySelector(`[data-number="${e.key}"]`);
    } else if (e.key === '+') {
      chooseOperation('+');
      targetButton = document.querySelector('[data-operator="+"]');
    } else if (e.key === '-') {
      chooseOperation('-');
      targetButton = document.querySelector('[data-operator="-"]');
    } else if (e.key === '*' || e.key.toLowerCase() === 'x') {
      chooseOperation('×');
      targetButton = document.querySelector('[data-operator="×"]');
    } else if (e.key === '/') {
      e.preventDefault();
      chooseOperation('÷');
      targetButton = document.querySelector('[data-operator="÷"]');
    } else if (e.key === 'Enter' || e.key === '=') {
      e.preventDefault();
      compute();
      targetButton = document.querySelector('[data-action="equals"]');
    } else if (e.key === 'Backspace') {
      deleteLast();
      targetButton = document.querySelector('[data-action="backspace"]');
    } else if (e.key === 'Escape' || e.key.toLowerCase() === 'c') {
      clearAll();
      targetButton = document.querySelector('[data-action="clear"]');
    } else if (e.key === '%') {
      applyPercent();
      targetButton = document.querySelector('[data-action="percent"]');
    }

    if (targetButton) {
      targetButton.classList.add('keyboard-active');
      setTimeout(() => targetButton.classList.remove('keyboard-active'), 120);
    }
  });

  function escapeHtml(str) {
    return String(str)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }

  // Initialize
  loadHistory();
  clearAll();
  populateConverterCategories();
  populateConstants();
});
