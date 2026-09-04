/**
 * Quantum Physics & Special Relativity Module
 * Computes Lorentz factor, relativistic kinetic energy, and de Broglie wavelengths.
 */

class QuantumRelativity {
  static SPEED_OF_LIGHT = 299792458; // m/s
  static PLANCK_CONSTANT = 6.62607015e-34; // J*s
  static ELECTRON_MASS = 9.1093837e-31; // kg

  static lorentzFactor(v) {
    if (v >= this.SPEED_OF_LIGHT) throw new Error('Velocity cannot reach or exceed speed of light');
    const beta = v / this.SPEED_OF_LIGHT;
    return 1 / Math.sqrt(1 - beta * beta);
  }

  static relativisticEnergy(mass, v) {
    const gamma = this.lorentzFactor(v);
    const totalEnergy = gamma * mass * Math.pow(this.SPEED_OF_LIGHT, 2);
    const restEnergy = mass * Math.pow(this.SPEED_OF_LIGHT, 2);
    return {
      gamma,
      restEnergy,
      totalEnergy,
      kineticEnergy: totalEnergy - restEnergy
    };
  }

  static deBroglieWavelength(mass, velocity) {
    const momentum = mass * velocity;
    if (momentum === 0) return Infinity;
    return this.PLANCK_CONSTANT / momentum;
  }
}

if (typeof module !== 'undefined') module.exports = QuantumRelativity;
