/**
 * Comprehensive Unit Conversion Engine
 * High-accuracy conversion across 14 scientific and everyday physical categories.
 */

const UNIT_CATEGORIES = {
  length: {
    name: 'Length',
    base: 'm',
    units: {
      m: { name: 'Meter (m)', factor: 1 },
      km: { name: 'Kilometer (km)', factor: 1000 },
      cm: { name: 'Centimeter (cm)', factor: 0.01 },
      mm: { name: 'Millimeter (mm)', factor: 0.001 },
      um: { name: 'Micrometer (µm)', factor: 1e-6 },
      nm: { name: 'Nanometer (nm)', factor: 1e-9 },
      pm: { name: 'Picometer (pm)', factor: 1e-12 },
      mi: { name: 'Mile (mi)', factor: 1609.344 },
      yd: { name: 'Yard (yd)', factor: 0.9144 },
      ft: { name: 'Foot (ft)', factor: 0.3048 },
      in: { name: 'Inch (in)', factor: 0.0254 },
      nmi: { name: 'Nautical Mile', factor: 1852 },
      ly: { name: 'Light Year (ly)', factor: 9.460730472e15 },
      au: { name: 'Astronomical Unit (AU)', factor: 1.495978707e11 },
      pc: { name: 'Parsec (pc)', factor: 3.085677581e16 },
      fathom: { name: 'Fathom', factor: 1.8288 },
      rod: { name: 'Rod', factor: 5.0292 },
      chain: { name: 'Chain', factor: 20.1168 },
      furlong: { name: 'Furlong', factor: 201.168 },
      angstrom: { name: 'Angstrom (Å)', factor: 1e-10 }
    }
  },

  mass: {
    name: 'Mass & Weight',
    base: 'kg',
    units: {
      kg: { name: 'Kilogram (kg)', factor: 1 },
      g: { name: 'Gram (g)', factor: 0.001 },
      mg: { name: 'Milligram (mg)', factor: 1e-6 },
      ug: { name: 'Microgram (µg)', factor: 1e-9 },
      tonne: { name: 'Metric Ton (t)', factor: 1000 },
      lb: { name: 'Pound (lb)', factor: 0.45359237 },
      oz: { name: 'Ounce (oz)', factor: 0.028349523125 },
      st: { name: 'Stone (st)', factor: 6.35029318 },
      carat: { name: 'Carat (ct)', factor: 0.0002 },
      grain: { name: 'Grain (gr)', factor: 0.00006479891 },
      short_ton: { name: 'Short Ton (US)', factor: 907.18474 },
      long_ton: { name: 'Long Ton (UK)', factor: 1016.0469088 },
      slug: { name: 'Slug', factor: 14.593903 },
      amu: { name: 'Atomic Mass Unit (u)', factor: 1.6605390666e-27 }
    }
  },

  temperature: {
    name: 'Temperature',
    base: 'celsius',
    custom: true,
    units: {
      c: { name: 'Celsius (°C)' },
      f: { name: 'Fahrenheit (°F)' },
      k: { name: 'Kelvin (K)' },
      r: { name: 'Rankine (°R)' }
    },
    convert: function (value, from, to) {
      if (from === to) return value;
      // Convert to Celsius first
      let c;
      switch (from) {
        case 'c': c = value; break;
        case 'f': c = (value - 32) * (5 / 9); break;
        case 'k': c = value - 273.15; break;
        case 'r': c = (value - 491.67) * (5 / 9); break;
        default: return NaN;
      }
      // Convert Celsius to destination
      switch (to) {
        case 'c': return c;
        case 'f': return (c * (9 / 5)) + 32;
        case 'k': return c + 273.15;
        case 'r': return (c + 273.15) * (9 / 5);
        default: return NaN;
      }
    }
  },

  area: {
    name: 'Area',
    base: 'sq_m',
    units: {
      sq_m: { name: 'Square Meter (m²)', factor: 1 },
      sq_km: { name: 'Square Kilometer (km²)', factor: 1e6 },
      sq_cm: { name: 'Square Centimeter (cm²)', factor: 1e-4 },
      sq_mm: { name: 'Square Millimeter (mm²)', factor: 1e-6 },
      sq_mi: { name: 'Square Mile (mi²)', factor: 2589988.110336 },
      sq_yd: { name: 'Square Yard (yd²)', factor: 0.83612736 },
      sq_ft: { name: 'Square Foot (ft²)', factor: 0.09290304 },
      sq_in: { name: 'Square Inch (in²)', factor: 0.00064516 },
      acre: { name: 'Acre (ac)', factor: 4046.8564224 },
      hectare: { name: 'Hectare (ha)', factor: 10000 },
      are: { name: 'Are (a)', factor: 100 },
      barn: { name: 'Barn (b)', factor: 1e-28 }
    }
  },

  volume: {
    name: 'Volume',
    base: 'l',
    units: {
      l: { name: 'Liter (L)', factor: 1 },
      ml: { name: 'Milliliter (mL)', factor: 0.001 },
      cu_m: { name: 'Cubic Meter (m³)', factor: 1000 },
      cu_cm: { name: 'Cubic Centimeter (cm³)', factor: 0.001 },
      cu_mm: { name: 'Cubic Millimeter (mm³)', factor: 1e-6 },
      gal_us: { name: 'Gallon (US)', factor: 3.785411784 },
      qt_us: { name: 'Quart (US)', factor: 0.946352946 },
      pt_us: { name: 'Pint (US)', factor: 0.473176473 },
      cup_us: { name: 'Cup (US)', factor: 0.2365882365 },
      fl_oz_us: { name: 'Fluid Ounce (US)', factor: 0.0295735295625 },
      tbsp_us: { name: 'Tablespoon (US)', factor: 0.01478676478125 },
      tsp_us: { name: 'Teaspoon (US)', factor: 0.00492892159375 },
      gal_uk: { name: 'Gallon (UK)', factor: 4.54609 },
      pt_uk: { name: 'Pint (UK)', factor: 0.56826125 },
      cu_ft: { name: 'Cubic Foot (ft³)', factor: 28.316846592 },
      cu_in: { name: 'Cubic Inch (in³)', factor: 0.016387064 },
      barrel_oil: { name: 'Oil Barrel (bbl)', factor: 158.987294928 }
    }
  },

  speed: {
    name: 'Speed',
    base: 'm_s',
    units: {
      m_s: { name: 'Meter/Second (m/s)', factor: 1 },
      km_h: { name: 'Kilometer/Hour (km/h)', factor: 1 / 3.6 },
      mph: { name: 'Mile/Hour (mph)', factor: 0.44704 },
      knot: { name: 'Knot (kn)', factor: 0.5144444444 },
      ft_s: { name: 'Foot/Second (ft/s)', factor: 0.3048 },
      mach: { name: 'Mach (at STP)', factor: 340.29 },
      c: { name: 'Speed of Light (c)', factor: 299792458 }
    }
  },

  time: {
    name: 'Time',
    base: 's',
    units: {
      s: { name: 'Second (s)', factor: 1 },
      ms: { name: 'Millisecond (ms)', factor: 0.001 },
      us: { name: 'Microsecond (µs)', factor: 1e-6 },
      ns: { name: 'Nanosecond (ns)', factor: 1e-9 },
      ps: { name: 'Picosecond (ps)', factor: 1e-12 },
      min: { name: 'Minute (min)', factor: 60 },
      h: { name: 'Hour (h)', factor: 3600 },
      d: { name: 'Day (d)', factor: 86400 },
      wk: { name: 'Week (wk)', factor: 604800 },
      mo: { name: 'Month (average 30.44d)', factor: 2629746 },
      yr: { name: 'Year (365.25d)', factor: 31557600 },
      decade: { name: 'Decade', factor: 315576000 },
      century: { name: 'Century', factor: 3155760000 }
    }
  },

  digital: {
    name: 'Digital Storage',
    base: 'b',
    units: {
      b: { name: 'Bit (b)', factor: 1 },
      B: { name: 'Byte (B)', factor: 8 },
      KB: { name: 'Kilobyte (KB)', factor: 8e3 },
      MB: { name: 'Megabyte (MB)', factor: 8e6 },
      GB: { name: 'Gigabyte (GB)', factor: 8e9 },
      TB: { name: 'Terabyte (TB)', factor: 8e12 },
      PB: { name: 'Petabyte (PB)', factor: 8e15 },
      KiB: { name: 'Kibibyte (KiB)', factor: 8 * 1024 },
      MiB: { name: 'Mebibyte (MiB)', factor: 8 * (1024 ** 2) },
      GiB: { name: 'Gibibyte (GiB)', factor: 8 * (1024 ** 3) },
      TiB: { name: 'Tebibyte (TiB)', factor: 8 * (1024 ** 4) },
      PiB: { name: 'Pebibyte (PiB)', factor: 8 * (1024 ** 5) }
    }
  },

  energy: {
    name: 'Energy',
    base: 'j',
    units: {
      j: { name: 'Joule (J)', factor: 1 },
      kj: { name: 'Kilojoule (kJ)', factor: 1000 },
      cal: { name: 'Calorie (cal)', factor: 4.184 },
      kcal: { name: 'Kilocalorie (kcal)', factor: 4184 },
      wh: { name: 'Watt-hour (Wh)', factor: 3600 },
      kwh: { name: 'Kilowatt-hour (kWh)', factor: 3.6e6 },
      ev: { name: 'Electronvolt (eV)', factor: 1.602176634e-19 },
      btu: { name: 'British Thermal Unit (BTU)', factor: 1055.056 },
      ft_lb: { name: 'Foot-pound (ft·lb)', factor: 1.355818 },
      erg: { name: 'Erg', factor: 1e-7 }
    }
  },

  pressure: {
    name: 'Pressure',
    base: 'pa',
    units: {
      pa: { name: 'Pascal (Pa)', factor: 1 },
      kpa: { name: 'Kilopascal (kPa)', factor: 1000 },
      mpa: { name: 'Megapascal (MPa)', factor: 1e6 },
      bar: { name: 'Bar', factor: 100000 },
      mbar: { name: 'Millibar (mbar)', factor: 100 },
      psi: { name: 'Pounds/Sq Inch (psi)', factor: 6894.757293 },
      atm: { name: 'Atmosphere (atm)', factor: 101325 },
      torr: { name: 'Torr (mmHg)', factor: 133.322368 },
      in_hg: { name: 'Inches of Mercury (inHg)', factor: 3386.389 }
    }
  },

  power: {
    name: 'Power',
    base: 'w',
    units: {
      w: { name: 'Watt (W)', factor: 1 },
      kw: { name: 'Kilowatt (kW)', factor: 1000 },
      mw: { name: 'Megawatt (MW)', factor: 1e6 },
      hp: { name: 'Horsepower (mechanical)', factor: 745.699872 },
      hp_metric: { name: 'Metric Horsepower (PS)', factor: 735.49875 },
      btu_h: { name: 'BTU/hour', factor: 0.293071 }
    }
  },

  frequency: {
    name: 'Frequency',
    base: 'hz',
    units: {
      hz: { name: 'Hertz (Hz)', factor: 1 },
      khz: { name: 'Kilohertz (kHz)', factor: 1000 },
      mhz: { name: 'Megahertz (MHz)', factor: 1e6 },
      ghz: { name: 'Gigahertz (GHz)', factor: 1e9 },
      rpm: { name: 'Revolutions/min (RPM)', factor: 1 / 60 },
      rad_s: { name: 'Radians/second (rad/s)', factor: 1 / (2 * Math.PI) }
    }
  }
};

class UnitConverter {
  constructor() {
    this.categories = UNIT_CATEGORIES;
  }

  getCategories() {
    return Object.keys(this.categories).map(key => ({
      id: key,
      name: this.categories[key].name
    }));
  }

  getUnitsForCategory(categoryId) {
    const cat = this.categories[categoryId];
    if (!cat) return [];
    return Object.keys(cat.units).map(unitKey => ({
      id: unitKey,
      name: cat.units[unitKey].name
    }));
  }

  convert(value, categoryId, fromUnit, toUnit) {
    const num = parseFloat(value);
    if (isNaN(num)) return NaN;
    if (fromUnit === toUnit) return num;

    const cat = this.categories[categoryId];
    if (!cat) return NaN;

    if (cat.custom && typeof cat.convert === 'function') {
      return cat.convert(num, fromUnit, toUnit);
    }

    const fromFactor = cat.units[fromUnit]?.factor;
    const toFactor = cat.units[toUnit]?.factor;

    if (fromFactor == null || toFactor == null) return NaN;

    // Convert from unit to base, then base to target
    const baseValue = num * fromFactor;
    const targetValue = baseValue / toFactor;

    // Clean precision
    return parseFloat(targetValue.toPrecision(12));
  }
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { UnitConverter, UNIT_CATEGORIES };
} else {
  window.UnitConverter = UnitConverter;
  window.UNIT_CATEGORIES = UNIT_CATEGORIES;
}
