/**
 * Comprehensive Scientific & Mathematical Constants
 * Standard values from NIST CODATA & IUPAC
 */

const SCIENTIFIC_CONSTANTS = [
  { symbol: 'π', name: 'Pi', value: 3.141592653589793, unit: '' },
  { symbol: 'e', name: "Euler's Number", value: 2.718281828459045, unit: '' },
  { symbol: 'c', name: 'Speed of Light in Vacuum', value: 299792458, unit: 'm/s' },
  { symbol: 'G', name: 'Newtonian Constant of Gravitation', value: 6.67430e-11, unit: 'm³/(kg·s²)' },
  { symbol: 'h', name: "Planck's Constant", value: 6.62607015e-34, unit: 'J·s' },
  { symbol: 'ħ', name: 'Reduced Planck Constant', value: 1.054571817e-34, unit: 'J·s' },
  { symbol: 'e_charge', name: 'Elementary Charge (q)', value: 1.602176634e-19, unit: 'C' },
  { symbol: 'N_A', name: "Avogadro's Number", value: 6.02214076e23, unit: 'mol⁻¹' },
  { symbol: 'k_B', name: 'Boltzmann Constant', value: 1.380649e-23, unit: 'J/K' },
  { symbol: 'R', name: 'Molar Gas Constant', value: 8.314462618, unit: 'J/(mol·K)' },
  { symbol: 'm_e', name: 'Electron Mass', value: 9.1093837015e-31, unit: 'kg' },
  { symbol: 'm_p', name: 'Proton Mass', value: 1.67262192369e-27, unit: 'kg' },
  { symbol: 'm_n', name: 'Neutron Mass', value: 1.67492749804e-27, unit: 'kg' },
  { symbol: 'ε_0', name: 'Vacuum Electric Permittivity', value: 8.8541878128e-12, unit: 'F/m' },
  { symbol: 'μ_0', name: 'Vacuum Magnetic Permeability', value: 1.25663706212e-6, unit: 'N/A²' },
  { symbol: 'σ', name: 'Stefan-Boltzmann Constant', value: 5.670374419e-8, unit: 'W/(m²·K⁴)' },
  { symbol: 'α', name: 'Fine-Structure Constant', value: 0.0072973525693, unit: '' },
  { symbol: 'R_inf', name: 'Rydberg Constant', value: 10973731.568160, unit: 'm⁻¹' },
  { symbol: 'F', name: 'Faraday Constant', value: 96485.33212, unit: 'C/mol' },
  { symbol: 'g', name: 'Standard Gravity', value: 9.80665, unit: 'm/s²' },
  { symbol: 'atm', name: 'Standard Atmosphere', value: 101325, unit: 'Pa' },
  { symbol: 'φ', name: 'Golden Ratio', value: 1.618033988749895, unit: '' }
];

if (typeof module !== 'undefined' && module.exports) {
  module.exports = SCIENTIFIC_CONSTANTS;
} else {
  window.SCIENTIFIC_CONSTANTS = SCIENTIFIC_CONSTANTS;
}
