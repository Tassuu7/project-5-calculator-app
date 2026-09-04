/**
 * Calculator Memory Module
 * Handles MC, MR, M+, M-, and MS with localStorage persistence.
 */

class MemoryManager {
  constructor(onMemoryChange) {
    this.storageKey = 'calc_memory_value';
    this.memoryValue = 0;
    this.onMemoryChange = onMemoryChange;
    this.load();
  }

  load() {
    try {
      const val = localStorage.getItem(this.storageKey);
      if (val !== null) {
        this.memoryValue = parseFloat(val) || 0;
      }
    } catch (e) {
      this.memoryValue = 0;
    }
    this.notify();
  }

  save() {
    try {
      localStorage.setItem(this.storageKey, this.memoryValue.toString());
    } catch (e) {}
    this.notify();
  }

  clear() {
    this.memoryValue = 0;
    this.save();
  }

  recall() {
    return this.memoryValue;
  }

  add(value) {
    const num = parseFloat(value);
    if (!isNaN(num)) {
      this.memoryValue += num;
      this.save();
    }
  }

  subtract(value) {
    const num = parseFloat(value);
    if (!isNaN(num)) {
      this.memoryValue -= num;
      this.save();
    }
  }

  store(value) {
    const num = parseFloat(value);
    if (!isNaN(num)) {
      this.memoryValue = num;
      this.save();
    }
  }

  hasMemory() {
    return this.memoryValue !== 0;
  }

  notify() {
    if (typeof this.onMemoryChange === 'function') {
      this.onMemoryChange(this.hasMemory(), this.memoryValue);
    }
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = MemoryManager;
} else {
  window.MemoryManager = MemoryManager;
}
