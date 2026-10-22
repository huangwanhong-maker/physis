// Values on a regular longitude-latitude grid: where a place falls on the grid,
// interpolation in space and time, and the area-weighted mean.
//
// Provisional JavaScript stand-in for physis's discretisation module (physis
// ADR 0003), first written for the Earth Atlas (ADR 0014 of The Science
// Infrastructures). Interpolation as estimation may later move to the
// infrastructure's planned inference library.
//
// A grid is { lon0, dlon, nlon, lat0, dlat, nlat, periodic }: the centre of
// column c is at longitude lon0 + c dlon, the centre of row r at latitude
// lat0 - r dlat (rows run north to south). Longitudes and latitudes are in
// degrees, the chart in which gridded data are written; the sphere they lie on
// is any sphere. Values are row-major arrays with NaN for no value.
//
// Interpolation is bilinear between the four surrounding grid points. Where some
// of the four have no value (along a coast, for an ocean field), the weights of
// those that have one are renormalised, provided they carry at least half of the
// total weight; otherwise the place has no value.

// Fractional column and row of a place, or null outside the grid. A place within
// half a cell of the first or last row is clamped onto it.
export function gridPosition(grid, lon, lat) {
  let x = (lon - grid.lon0) / grid.dlon;
  if (grid.periodic) x = ((x % grid.nlon) + grid.nlon) % grid.nlon;
  else if (x < -0.5 || x > grid.nlon - 0.5) return null;
  let y = (grid.lat0 - lat) / grid.dlat;
  if (y < -0.5 || y > grid.nlat - 0.5) return null;
  return [x, Math.max(0, Math.min(grid.nlat - 1, y))];
}

export function bilinear(values, grid, lon, lat) {
  const pos = gridPosition(grid, lon, lat);
  if (!pos) return NaN;
  const [x, y] = pos;
  const nlon = grid.nlon;
  let c0 = Math.floor(x), r0 = Math.floor(y);
  const fx = x - c0, fy = y - r0;
  let c1 = c0 + 1;
  const r1 = Math.min(r0 + 1, grid.nlat - 1);
  if (c1 >= nlon) c1 = grid.periodic ? c1 - nlon : nlon - 1;
  if (c0 < 0) c0 = 0;
  const v00 = values[r0 * nlon + c0], v01 = values[r0 * nlon + c1];
  const v10 = values[r1 * nlon + c0], v11 = values[r1 * nlon + c1];
  const w00 = (1 - fx) * (1 - fy), w01 = fx * (1 - fy), w10 = (1 - fx) * fy, w11 = fx * fy;
  let sum = 0, weight = 0;
  if (v00 === v00) { sum += w00 * v00; weight += w00; }
  if (v01 === v01) { sum += w01 * v01; weight += w01; }
  if (v10 === v10) { sum += w10 * v10; weight += w10; }
  if (v11 === v11) { sum += w11 * v11; weight += w11; }
  return weight >= 0.5 ? sum / weight : NaN;
}

// The bracketing times of an instant in a sorted list of times (ms), and the
// weight of the later one. Outside the list the nearest end is taken, and
// `offset` says by how much (ms; negative before the first time).
export function timeBracket(times, ms, nearest = false) {
  const n = times.length;
  if (n === 0) return { i0: -1, i1: -1, w: 0, shown: null, offset: 0 };
  if (ms <= times[0]) return { i0: 0, i1: 0, w: 0, shown: times[0], offset: ms - times[0] };
  if (ms >= times[n - 1]) return { i0: n - 1, i1: n - 1, w: 0, shown: times[n - 1], offset: ms - times[n - 1] };
  let lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (times[mid] <= ms) lo = mid; else hi = mid;
  }
  const w = (ms - times[lo]) / (times[hi] - times[lo]);
  if (nearest) {
    const i = w < 0.5 ? lo : hi;
    return { i0: i, i1: i, w: 0, shown: times[i], offset: ms - times[i] };
  }
  return { i0: lo, i1: hi, w, shown: ms, offset: 0 };
}

// Area-weighted (cos latitude) mean of a gridded field, ignoring missing values.
export function areaMean(values, grid) {
  let sum = 0, weight = 0;
  for (let r = 0; r < grid.nlat; r++) {
    const w = Math.cos(((grid.lat0 - r * grid.dlat) * Math.PI) / 180);
    for (let c = 0; c < grid.nlon; c++) {
      const v = values[r * grid.nlon + c];
      if (v === v) { sum += w * v; weight += w; }
    }
  }
  return weight > 0 ? sum / weight : NaN;
}
