// Checks of physis's provisional JavaScript stand-in (ADR 0003).
// Run from the repository root:  node js/tests/run-tests.mjs

import { gridPosition, bilinear, timeBracket, areaMean, gaussLegendre, SphericalHarmonics, FFT } from '../index.js';

let failures = 0, checks = 0;
function check(name, value, limit, unit = '') {
  checks++;
  const ok = Number.isFinite(value) && Math.abs(value) <= limit;
  if (!ok) failures++;
  console.log(`${ok ? 'ok  ' : 'FAIL'}  ${name}: ${Number.isFinite(value) ? value.toExponential(3) : value} ${unit} (limit ${limit} ${unit})`);
}
function assert(name, condition, detail = '') {
  checks++;
  if (!condition) failures++;
  console.log(`${condition ? 'ok  ' : 'FAIL'}  ${name}${detail ? `: ${detail}` : ''}`);
}
function seeded(seed) {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
}
const random = seeded(20261022);

// ---- grids -----------------------------------------------------------------
{
  // A linear function of longitude and latitude is reproduced exactly by bilinear interpolation.
  const grid = { lon0: -180, dlon: 2, nlon: 180, lat0: 89, dlat: 2, nlat: 90, periodic: true };
  const values = new Float32Array(grid.nlon * grid.nlat);
  for (let r = 0; r < grid.nlat; r++) for (let c = 0; c < grid.nlon; c++) values[r * grid.nlon + c] = 3 * (grid.lat0 - r * grid.dlat) + 0.5 * c;
  let worst = 0;
  for (let i = 0; i < 500; i++) {
    const lon = -180 + random() * 356, lat = -88 + random() * 176;
    worst = Math.max(worst, Math.abs(bilinear(values, grid, lon, lat) - (3 * lat + 0.5 * ((lon + 180) / 2))));
  }
  check('grid: bilinear interpolation of a linear field', worst, 1e-3);
  // Periodic: between the last column (178 E) and the first (180 W), halfway.
  const wrap = bilinear(values, grid, 179, 1);
  check('grid: interpolation across the antimeridian', wrap - (3 * 1 + 0.5 * (179 + 0) / 2), 1e-4);
  const pos = gridPosition(grid, -181, 0);
  check('grid: longitude wraps', pos[0] - 179.5, 1e-9);
  // A missing corner: weights renormalised while they carry half the weight.
  const g2 = { lon0: 0, dlon: 1, nlon: 2, lat0: 1, dlat: 1, nlat: 2, periodic: false };
  const v2 = new Float32Array([1, NaN, 3, 5]);
  check('grid: renormalised weights near a missing value', bilinear(v2, g2, 0.25, 0.25) - ((0.75 * 0.25 * 1 + 0.75 * 0.75 * 3 + 0.25 * 0.75 * 5) / (1 - 0.25 * 0.25)), 1e-6);
  assert('grid: no value next to a missing value', Number.isNaN(bilinear(v2, g2, 0.9, 0.9)));
  const times = [0, 3600e3, 7200e3];
  const b = timeBracket(times, 5400e3);
  assert('grid: time bracket', b.i0 === 1 && b.i1 === 2 && Math.abs(b.w - 0.5) < 1e-12);
  const c = timeBracket(times, 9000e3);
  assert('grid: time after the last is clamped, with its offset', c.i0 === 2 && c.offset === 1800e3);
  assert('grid: nearest time for means', timeBracket(times, 2000e3, true).i0 === 1);
  // Area mean of the field cos^0 = 1 is 1; of sin(lat) is 0.
  const ones = new Float32Array(values.length).fill(1);
  check('grid: area mean of a constant', areaMean(ones, grid) - 1, 1e-6);
}

// ---- spherical harmonics ----------------------------------------------------
{
  // A sphere of the Earth's size (Williamson et al. 1992), so that the tolerances
  // speak of realistic magnitudes; physis itself knows no particular radius.
  const RADIUS = 6371220;
  const { mu, w } = gaussLegendre(64);
  check('spectral: Gauss-Legendre weights sum to 2', w.reduce((s, x) => s + x, 0) - 2, 1e-13);
  check('spectral: quadrature of mu^2 is 2/3', w.reduce((s, x, i) => s + x * mu[i] ** 2, 0) - 2 / 3, 1e-13);
  const sh = new SphericalHarmonics(42);
  let worst = 0;
  for (let k = 0; k < sh.size; k++) {
    let s = 0;
    for (let j = 0; j < sh.nlat; j++) s += sh.w[j] * sh.P[j * sh.size + k] ** 2;
    worst = Math.max(worst, Math.abs(s - 1));
  }
  check('spectral: associated Legendre functions are normalised (T42)', worst, 1e-12);
  const rnd = seeded(42);
  const r = sh.zeros();
  for (let k = 0; k < sh.size; k++) { r.re[k] = rnd() - 0.5; r.im[k] = sh.m[k] === 0 ? 0 : rnd() - 0.5; }
  const back = sh.analyse(sh.synthesise(r));
  let rt = 0;
  for (let k = 0; k < sh.size; k++) rt = Math.max(rt, Math.abs(back.re[k] - r.re[k]), Math.abs(back.im[k] - r.im[k]));
  check('spectral: synthesis then analysis returns the coefficients', rt, 1e-12);
  // Solid-body rotation u = u0 cos(phi): psi = -a u0 mu, vorticity 2 u0 mu / a.
  const u0 = 20, psi = sh.zeros();
  psi.re[sh.index(0, 1)] = (-RADIUS * u0) / Math.sqrt(1.5);
  const { U, V } = sh.windFromStreamFunction(psi, RADIUS);
  let eu = 0, ev = 0;
  for (let j = 0; j < sh.nlat; j++) for (let k = 0; k < sh.nlon; k++) {
    eu = Math.max(eu, Math.abs(U[j * sh.nlon + k] - u0 * (1 - sh.mu[j] ** 2)));
    ev = Math.max(ev, Math.abs(V[j * sh.nlon + k]));
  }
  check('spectral: wind of a solid-body rotation', Math.max(eu, ev), 1e-12, 'm/s');
  const zg = sh.synthesise(sh.curl(U, V, RADIUS));
  let ez = 0;
  for (let j = 0; j < sh.nlat; j++) ez = Math.max(ez, Math.abs(zg[j * sh.nlon] - (2 * u0 * sh.mu[j]) / RADIUS));
  check('spectral: curl of that wind is 2 u0 sin(phi) / a', ez / ((2 * u0) / RADIUS), 1e-12);
  // The curl of the wind of a random stream function is its Laplacian; its divergence is zero.
  const ps = sh.zeros();
  for (let k = 0; k < sh.size; k++) if (sh.n[k] > 0 && sh.n[k] < 30) { ps.re[k] = 1e7 * (rnd() - 0.5); ps.im[k] = sh.m[k] === 0 ? 0 : 1e7 * (rnd() - 0.5); }
  const wind = sh.windFromStreamFunction(ps, RADIUS);
  const curl = sh.curl(wind.U, wind.V, RADIUS), div = sh.divergence(wind.U, wind.V, RADIUS);
  const lap = sh.scale({ re: Float64Array.from(ps.re), im: Float64Array.from(ps.im) }, (n) => sh.laplacianEigen(n, RADIUS));
  let ec = 0, ed = 0, scale = 0;
  for (let k = 0; k < sh.size; k++) {
    ec = Math.max(ec, Math.hypot(curl.re[k] - lap.re[k], curl.im[k] - lap.im[k]));
    ed = Math.max(ed, Math.hypot(div.re[k], div.im[k]));
    scale = Math.max(scale, Math.hypot(lap.re[k], lap.im[k]));
  }
  check('spectral: curl of the wind of psi is the Laplacian of psi', ec / scale, 1e-10);
  check('spectral: the wind of psi has no divergence', ed / scale, 1e-10);
  let refused = false;
  try { sh.laplacianEigen(3); } catch { refused = true; }
  assert('spectral: an operation without the radius of its sphere is refused', refused);
}

// ---- the fast Fourier transform ------------------------------------------------------
{
  const n = 64, fft = new FFT(n), rnd = seeded(3);
  const x = Float64Array.from({ length: n }, () => rnd() - 0.5), y = Float64Array.from({ length: n }, () => rnd() - 0.5);
  const re = Float64Array.from(x), im = Float64Array.from(y);
  fft.transform(re, im, -1);
  let err = 0;
  for (let m = 0; m < n; m++) {
    let a = 0, b = 0;
    for (let k = 0; k < n; k++) { const t = (-2 * Math.PI * m * k) / n; a += x[k] * Math.cos(t) - y[k] * Math.sin(t); b += x[k] * Math.sin(t) + y[k] * Math.cos(t); }
    err = Math.max(err, Math.abs(a - re[m]), Math.abs(b - im[m]));
  }
  check('fft: 64 points agree with the direct sum', err, 1e-12);
  fft.transform(re, im, 1);
  check('fft: forward then backward returns n times the data', Math.max(...x.map((v, k) => Math.abs(re[k] / n - v))), 1e-14);
}

console.log(`\n${checks - failures} of ${checks} checks passed.`);
process.exit(failures ? 1 : 0);
