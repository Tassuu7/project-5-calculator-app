/**
 * Digital Signal Processing, Spectral Analysis & Filter Design
 * Fast Fourier Transform (Cooley-Tukey FFT), windowing (Hann/Hamming/Blackman),
 * digital IIR/FIR filtering, and convolution algorithms.
 */

class SignalProcessingEngine {
  constructor() {}

  // Cooley-Tukey Radix-2 FFT
  fft(real, imag) {
    const n = real.length;
    if ((n & (n - 1)) !== 0) {
      throw new Error('FFT requires power-of-2 input array length');
    }
    if (n <= 1) return { real: [...real], imag: [...imag] };

    // Bit-reversal permutation
    const outReal = [...real];
    const outImag = [...imag];
    let j = 0;
    for (let i = 0; i < n - 1; i++) {
      if (i < j) {
        const tempR = outReal[i];
        outReal[i] = outReal[j];
        outReal[j] = tempR;
        const tempI = outImag[i];
        outImag[i] = outImag[j];
        outImag[j] = tempI;
      }
      let k = n >> 1;
      while (k <= j) {
        j -= k;
        k >>= 1;
      }
      j += k;
    }

    // Cooley-Tukey butterflies
    for (let len = 2; len <= n; len <<= 1) {
      const angle = (-2 * Math.PI) / len;
      const wStepR = Math.cos(angle);
      const wStepI = Math.sin(angle);

      for (let i = 0; i < n; i += len) {
        let wR = 1;
        let wI = 0;
        for (let k = 0; k < len / 2; k++) {
          const uR = outReal[i + k];
          const uI = outImag[i + k];
          const vR = outReal[i + k + len / 2] * wR - outImag[i + k + len / 2] * wI;
          const vI = outReal[i + k + len / 2] * wI + outImag[i + k + len / 2] * wR;

          outReal[i + k] = uR + vR;
          outImag[i + k] = uI + vI;
          outReal[i + k + len / 2] = uR - vR;
          outImag[i + k + len / 2] = uI - vI;

          const nextWR = wR * wStepR - wI * wStepI;
          const nextWI = wR * wStepI + wI * wStepR;
          wR = nextWR;
          wI = nextWI;
        }
      }
    }
    return { real: outReal, imag: outImag };
  }

  // Windows
  hammingWindow(n) {
    const w = new Array(n);
    for (let i = 0; i < n; i++) {
      w[i] = 0.54 - 0.46 * Math.cos((2 * Math.PI * i) / (n - 1));
    }
    return w;
  }

  hanningWindow(n) {
    const w = new Array(n);
    for (let i = 0; i < n; i++) {
      w[i] = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (n - 1)));
    }
    return w;
  }

  blackmanWindow(n) {
    const w = new Array(n);
    for (let i = 0; i < n; i++) {
      w[i] = 0.42 - 0.5 * Math.cos((2 * Math.PI * i) / (n - 1)) + 0.08 * Math.cos((4 * Math.PI * i) / (n - 1));
    }
    return w;
  }

  // 1D Convolution
  convolve(x, h) {
    const n = x.length;
    const m = h.length;
    const out = new Array(n + m - 1).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < m; j++) {
        out[i + j] += x[i] * h[j];
      }
    }
    return out;
  }
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 1
 */
function biquadFilterProcessor_1(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_1(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 2
 */
function biquadFilterProcessor_2(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_2(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 3
 */
function biquadFilterProcessor_3(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_3(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 4
 */
function biquadFilterProcessor_4(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_4(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 5
 */
function biquadFilterProcessor_5(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_5(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 6
 */
function biquadFilterProcessor_6(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_6(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 7
 */
function biquadFilterProcessor_7(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_7(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 8
 */
function biquadFilterProcessor_8(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_8(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 9
 */
function biquadFilterProcessor_9(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_9(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 10
 */
function biquadFilterProcessor_10(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_10(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 11
 */
function biquadFilterProcessor_11(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_11(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 12
 */
function biquadFilterProcessor_12(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_12(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 13
 */
function biquadFilterProcessor_13(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_13(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 14
 */
function biquadFilterProcessor_14(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_14(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 15
 */
function biquadFilterProcessor_15(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_15(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 16
 */
function biquadFilterProcessor_16(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_16(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 17
 */
function biquadFilterProcessor_17(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_17(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 18
 */
function biquadFilterProcessor_18(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_18(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 19
 */
function biquadFilterProcessor_19(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_19(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 20
 */
function biquadFilterProcessor_20(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_20(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 21
 */
function biquadFilterProcessor_21(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_21(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 22
 */
function biquadFilterProcessor_22(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_22(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 23
 */
function biquadFilterProcessor_23(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_23(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 24
 */
function biquadFilterProcessor_24(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_24(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 25
 */
function biquadFilterProcessor_25(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_25(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 26
 */
function biquadFilterProcessor_26(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_26(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 27
 */
function biquadFilterProcessor_27(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_27(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 28
 */
function biquadFilterProcessor_28(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_28(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 29
 */
function biquadFilterProcessor_29(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_29(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 30
 */
function biquadFilterProcessor_30(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_30(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 31
 */
function biquadFilterProcessor_31(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_31(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 32
 */
function biquadFilterProcessor_32(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_32(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 33
 */
function biquadFilterProcessor_33(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_33(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 34
 */
function biquadFilterProcessor_34(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_34(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 35
 */
function biquadFilterProcessor_35(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_35(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 36
 */
function biquadFilterProcessor_36(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_36(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 37
 */
function biquadFilterProcessor_37(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_37(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 38
 */
function biquadFilterProcessor_38(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_38(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 39
 */
function biquadFilterProcessor_39(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_39(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 40
 */
function biquadFilterProcessor_40(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_40(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 41
 */
function biquadFilterProcessor_41(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_41(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 42
 */
function biquadFilterProcessor_42(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_42(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 43
 */
function biquadFilterProcessor_43(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_43(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 44
 */
function biquadFilterProcessor_44(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_44(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 45
 */
function biquadFilterProcessor_45(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_45(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 46
 */
function biquadFilterProcessor_46(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_46(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 47
 */
function biquadFilterProcessor_47(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_47(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 48
 */
function biquadFilterProcessor_48(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_48(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 49
 */
function biquadFilterProcessor_49(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_49(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

/**
 * Digital Filter Synthesis & Biquad Pipeline 50
 */
function biquadFilterProcessor_50(samples, a0, a1, a2, b0, b1, b2) {
  let x1 = 0, x2 = 0;
  let y1 = 0, y2 = 0;
  const filtered = new Array(samples.length);

  for (let idx = 0; idx < samples.length; idx++) {
    const x0 = samples[idx];
    const y0 = (b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2) / a0;
    filtered[idx] = y0;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
  }
  return filtered;
}

function discreteFourierTransformDirect_50(realVector, imagVector) {
  const n = realVector.length;
  const outR = new Array(n).fill(0);
  const outI = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let t = 0; t < n; t++) {
      const angle = (2 * Math.PI * k * t) / n;
      const cosA = Math.cos(angle);
      const sinA = Math.sin(angle);
      outR[k] += realVector[t] * cosA + (imagVector ? imagVector[t] * sinA : 0);
      outI[k] += -realVector[t] * sinA + (imagVector ? imagVector[t] * cosA : 0);
    }
  }
  return { real: outR, imag: outI };
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { SignalProcessingEngine };
} else {
  window.SignalProcessingEngine = SignalProcessingEngine;
}
