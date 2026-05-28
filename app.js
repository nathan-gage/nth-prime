const statusEl = document.querySelector("#wasmStatus");
const nInput = document.querySelector("#nInput");
const iterationsInput = document.querySelector("#iterationsInput");
const runOnceButton = document.querySelector("#runOnceButton");
const runBenchButton = document.querySelector("#runBenchButton");
const traceButton = document.querySelector("#traceButton");
const primeOutput = document.querySelector("#primeOutput");
const singleOutput = document.querySelector("#singleOutput");
const avgOutput = document.querySelector("#avgOutput");
const throughputOutput = document.querySelector("#throughputOutput");
const historyBody = document.querySelector("#historyBody");
const canvas = document.querySelector("#sieveCanvas");
const traceSummary = document.querySelector("#traceSummary");
const traceGrid = document.querySelector("#traceGrid");
const timerResolution = document.querySelector("#timerResolution");

let wasm = null;
let historyHasRows = false;
let primeCacheLimit = 0;
let primeCache = [];

const MAX_N = 216_289_611_853_439_384n;
const MAX_CANVAS_CENTER = 1_000_000_000_000;
const MAX_VISUAL_ROOT = 50_000;
const TIMER_RESOLUTION_MS = detectTimerResolution();

timerResolution.textContent = `about ${formatMs(TIMER_RESOLUTION_MS)} per visible timer tick`;

function detectTimerResolution() {
  let best = Number.POSITIVE_INFINITY;
  let previous = performance.now();
  for (let i = 0; i < 10_000; i += 1) {
    const current = performance.now();
    const delta = current - previous;
    if (delta > 0 && delta < best) {
      best = delta;
    }
    previous = current;
  }
  return Number.isFinite(best) ? best : 0.001;
}

function formatInteger(value) {
  return value.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
}

function formatSignedInteger(value) {
  if (value < 0n) {
    return `-${formatInteger(-value)}`;
  }
  return `+${formatInteger(value)}`;
}

function absBigInt(value) {
  return value < 0n ? -value : value;
}

function formatPositiveMs(value) {
  const ns = value * 1_000_000;
  if (value < 0.001) {
    return ns < 0.1 ? "<0.1 ns" : `${ns.toFixed(ns < 10 ? 1 : 0)} ns`;
  }
  if (value < 1) {
    const us = value * 1000;
    return `${us.toFixed(us < 10 ? 2 : 1)} us`;
  }
  return `${value.toFixed(value < 10 ? 3 : 2)} ms`;
}

function formatMs(value, respectTimer = false) {
  if (!Number.isFinite(value)) {
    return "-";
  }
  if (respectTimer && value < TIMER_RESOLUTION_MS) {
    return `<${formatPositiveMs(TIMER_RESOLUTION_MS)}`;
  }
  return formatPositiveMs(value);
}

function formatRate(callsPerSecond, lowerBound = false) {
  if (!Number.isFinite(callsPerSecond)) {
    return "-";
  }
  const prefix = lowerBound ? ">" : "";
  if (callsPerSecond >= 1_000_000) {
    return `${prefix}${(callsPerSecond / 1_000_000).toFixed(2)}M/s`;
  }
  if (callsPerSecond >= 1_000) {
    return `${prefix}${(callsPerSecond / 1_000).toFixed(1)}k/s`;
  }
  return `${prefix}${Math.round(callsPerSecond)}/s`;
}

function parseBigIntInput(input, label) {
  const normalized = input.value.replaceAll(",", "").trim();
  if (!/^\d+$/.test(normalized)) {
    throw new Error(`${label} must be a non-negative decimal integer.`);
  }
  const parsed = BigInt(normalized);
  if (parsed > MAX_N) {
    throw new Error(`${label} exceeds the supported 64-bit search boundary.`);
  }
  return parsed;
}

function parseIterations() {
  const normalized = iterationsInput.value.replaceAll(",", "").trim();
  if (!/^\d+$/.test(normalized)) {
    throw new Error("iterations must be a positive integer.");
  }
  const parsed = Number(normalized);
  if (!Number.isSafeInteger(parsed) || parsed < 1 || parsed > 1_000_000) {
    throw new Error("iterations must be between 1 and 1,000,000.");
  }
  return parsed;
}

function setStatus(message, kind = "") {
  statusEl.textContent = message;
  statusEl.className = `build-status ${kind}`.trim();
}

function setBusy(busy) {
  runOnceButton.disabled = busy || wasm === null;
  runBenchButton.disabled = busy || wasm === null;
  traceButton.disabled = busy || wasm === null;
  for (const button of document.querySelectorAll(".preset")) {
    button.disabled = busy;
  }
}

function nextFrame() {
  return new Promise((resolve) => requestAnimationFrame(() => resolve()));
}

function smallPrimes(limit) {
  if (limit <= primeCacheLimit) {
    return primeCache.filter((p) => p <= limit);
  }

  const composite = new Uint8Array(limit + 1);
  const primes = [];
  for (let i = 2; i <= limit; i += 1) {
    if (composite[i] === 0) {
      primes.push(i);
      if (i * i <= limit) {
        for (let j = i * i; j <= limit; j += i) {
          composite[j] = 1;
        }
      }
    }
  }
  primeCacheLimit = limit;
  primeCache = primes;
  return primes;
}

function estimateForCanvas(n) {
  let estimate = 101;
  if (wasm !== null && n <= 20_000_000_000n) {
    estimate = Number(wasm.estimate_nth_prime(n));
  } else {
    const nn = Number(n);
    if (Number.isFinite(nn) && nn >= 6) {
      const l = Math.log(nn);
      const ll = Math.log(l);
      estimate = Math.floor(nn * (l + ll - 1));
    }
  }
  if (!Number.isFinite(estimate) || estimate < 101) {
    return 101;
  }
  return Math.min(Math.floor(estimate), MAX_CANVAS_CENTER);
}

function drawSieveWindow() {
  const ctx = canvas.getContext("2d");
  const ratio = window.devicePixelRatio || 1;
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  canvas.width = Math.floor(width * ratio);
  canvas.height = Math.floor(height * ratio);
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);

  ctx.fillStyle = "#101820";
  ctx.fillRect(0, 0, width, height);

  let n = 1_000_000_000n;
  try {
    n = parseBigIntInput(nInput, "n");
  } catch {
    // Keep the previous valid-looking default visualization.
  }

  const center = estimateForCanvas(n);
  const columns = width < 620 ? 14 : 28;
  const rows = 7;
  const count = columns * rows;
  const cellGap = 4;
  const marginX = 18;
  const marginTop = 66;
  const marginBottom = 24;
  const cellWidth = (width - 2 * marginX - (columns - 1) * cellGap) / columns;
  const cellHeight = (height - marginTop - marginBottom - (rows - 1) * cellGap) / rows;
  const low = Math.max(3, Math.floor(center / 2) * 2 - count + 1);
  const firstOdd = low % 2 === 0 ? low + 1 : low;
  const high = firstOdd + 2 * (count - 1);
  const root = Math.floor(Math.sqrt(high));
  const visualRoot = Math.min(root, MAX_VISUAL_ROOT);
  const basePrimes = smallPrimes(visualRoot);
  const composite = new Uint8Array(count);

  for (const p of basePrimes) {
    if (p === 2) {
      continue;
    }
    let first = p * p;
    if (first < firstOdd) {
      first = Math.ceil(firstOdd / p) * p;
    }
    if (first % 2 === 0) {
      first += p;
    }
    for (let value = first; value <= high; value += 2 * p) {
      const index = (value - firstOdd) / 2;
      if (index >= 0 && index < count) {
        composite[index] = p;
      }
    }
  }

  ctx.fillStyle = "#dbe7ea";
  ctx.font = "600 15px ui-monospace, SFMono-Regular, Menlo, monospace";
  ctx.fillText(`odd-only segment near x = ${formatInteger(BigInt(Math.floor(center)))}`, 18, 26);
  ctx.fillStyle = "#8fa4ad";
  ctx.font = "12px ui-monospace, SFMono-Regular, Menlo, monospace";
  ctx.fillText(`window [${formatInteger(BigInt(firstOdd))}, ${formatInteger(BigInt(high))}], visual base primes <= ${formatInteger(BigInt(visualRoot))}`, 18, 47);

  for (let i = 0; i < count; i += 1) {
    const row = Math.floor(i / columns);
    const col = i % columns;
    const x = marginX + col * (cellWidth + cellGap);
    const y = marginTop + row * (cellHeight + cellGap);
    const value = firstOdd + 2 * i;
    const isComposite = composite[i] !== 0;

    ctx.fillStyle = isComposite ? "#263640" : "#00a896";
    ctx.fillRect(x, y, cellWidth, cellHeight);
    if (isComposite) {
      ctx.fillStyle = "rgba(255, 193, 92, 0.9)";
      ctx.fillRect(x, y + cellHeight - 4, cellWidth, 4);
    }

    if (cellWidth >= 30) {
      ctx.fillStyle = isComposite ? "#c7d2d8" : "#06251f";
      ctx.font = "10px ui-monospace, SFMono-Regular, Menlo, monospace";
      ctx.fillText(String(value % 1000).padStart(3, "0"), x + 5, y + Math.min(cellHeight - 7, 17));
    }
  }
}

function addHistory({ n, result, iterations, totalMs, avgMs }) {
  if (!historyHasRows) {
    historyBody.textContent = "";
    historyHasRows = true;
  }

  const row = document.createElement("tr");
  const values = [
    formatInteger(n),
    formatInteger(result),
    formatInteger(BigInt(iterations)),
    totalMs.toFixed(3),
    avgMs.toFixed(6),
  ];

  for (const value of values) {
    const cell = document.createElement("td");
    cell.textContent = value;
    row.append(cell);
  }

  historyBody.prepend(row);
  while (historyBody.rows.length > 8) {
    historyBody.deleteRow(historyBody.rows.length - 1);
  }
}

function makeTraceCell(label, value, detail = "") {
  const item = document.createElement("div");
  item.className = "trace-cell";

  const key = document.createElement("span");
  key.textContent = label;

  const strong = document.createElement("strong");
  strong.textContent = value;

  item.append(key, strong);
  if (detail !== "") {
    const note = document.createElement("small");
    note.textContent = detail;
    item.append(note);
  }
  return item;
}

function renderTrace(cells) {
  traceGrid.replaceChildren(...cells.map((cell) => makeTraceCell(...cell)));
}

function correctedCoordinate(estimate, delta) {
  if (estimate <= 2n) {
    return 2n;
  }
  const adjusted = Number(estimate) - Number(delta) * Math.log(Number(estimate));
  if (!Number.isFinite(adjusted) || adjusted < 2) {
    return 2n;
  }
  return BigInt(Math.floor(adjusted));
}

async function computeOnce() {
  if (wasm === null) {
    return;
  }

  setBusy(true);
  try {
    const n = parseBigIntInput(nInput, "n");
    setStatus("Computing isolated call...");
    await nextFrame();

    const start = performance.now();
    const result = wasm.nth_prime_isolated(n);
    const stop = performance.now();
    const elapsed = stop - start;

    primeOutput.textContent = formatInteger(result);
    singleOutput.textContent = formatMs(elapsed, true);
    setStatus("WASM ready", "ready");
    drawSieveWindow();
  } catch (error) {
    setStatus(error.message, "error");
  } finally {
    setBusy(false);
  }
}

async function benchmarkIsolated() {
  if (wasm === null) {
    return;
  }

  setBusy(true);
  try {
    const n = parseBigIntInput(nInput, "n");
    const iterations = parseIterations();

    setStatus("Benchmarking isolated calls...");
    await nextFrame();
    let result = 0n;
    const start = performance.now();
    for (let i = 0; i < iterations; i += 1) {
      result = wasm.nth_prime_isolated(n);
    }
    const stop = performance.now();

    const totalMs = stop - start;
    const avgMs = totalMs / iterations;
    const belowTimer = avgMs < TIMER_RESOLUTION_MS;
    const callsPerSecond = 1000 / Math.max(avgMs, TIMER_RESOLUTION_MS);

    primeOutput.textContent = formatInteger(result);
    singleOutput.textContent = formatMs(avgMs, true);
    avgOutput.textContent = formatMs(avgMs, true);
    throughputOutput.textContent = formatRate(callsPerSecond, belowTimer);
    addHistory({ n, result, iterations, totalMs, avgMs });
    setStatus("WASM ready", "ready");
    drawSieveWindow();
  } catch (error) {
    setStatus(error.message, "error");
  } finally {
    setBusy(false);
  }
}

async function traceCorrection() {
  if (wasm === null) {
    return;
  }

  setBusy(true);
  try {
    const n = parseBigIntInput(nInput, "n");
    if (n === 0n) {
      renderTrace([
        ["n", "0", "Boundary before the first prime"],
        ["p_n", "0", "Defined by this program as pi-inverse at rank zero"],
      ]);
      traceSummary.textContent = "Rank zero returns zero directly; no count or local sieve is needed.";
      setStatus("WASM ready", "ready");
      return;
    }

    setStatus("Tracing exact ranks...");
    await nextFrame();

    const estimate = wasm.estimate_nth_prime(n);
    const pi0 = wasm.prime_pi_exact(estimate);
    const delta0 = pi0 - n;
    const x1 = correctedCoordinate(estimate, delta0);
    const pi1 = x1 === estimate ? pi0 : wasm.prime_pi_exact(x1);
    const delta1 = pi1 - n;
    const accepted = absBigInt(delta1) < absBigInt(delta0);
    const coordinate = accepted ? x1 : estimate;
    const rank = accepted ? pi1 : pi0;
    const residual = absBigInt(rank - n);
    const direction = rank < n ? "forward" : "backward";
    const approximateGap = Math.log(Math.max(3, Number(coordinate)));
    const oddWindow = Math.max(1, Math.ceil((Number(residual) + 32) * approximateGap / 2));

    renderTrace([
      ["x0 estimate", formatInteger(estimate), "Cipolla-style analytic coordinate"],
      ["pi(x0)", formatInteger(pi0), "Exact isolated Lehmer count"],
      ["rank error", formatSignedInteger(delta0), "Positive means x is past the target rank"],
      ["x1 correction", formatInteger(x1), "x0 - (pi(x0) - n) log x0"],
      ["pi(x1)", formatInteger(pi1), "Second exact isolated count for the trace"],
      ["accepted", accepted ? "yes" : "no", "Only accepted if the absolute rank error shrinks"],
      ["local direction", direction, `${formatInteger(residual)} rank steps remain`],
      ["odd window estimate", formatInteger(BigInt(oddWindow)), "Approximate candidates touched by final segment"],
    ]);

    traceSummary.textContent =
      `The trace uses extra isolated exact counts for observability. The selected coordinate ` +
      `then sieves ${direction} over about ${formatInteger(residual)} residual prime ranks.`;
    setStatus("WASM ready", "ready");
    drawSieveWindow();
  } catch (error) {
    setStatus(error.message, "error");
  } finally {
    setBusy(false);
  }
}

async function loadWasm() {
  setBusy(true);
  try {
    const response = await fetch("nprime.wasm", { cache: "no-store" });
    if (!response.ok) {
      throw new Error(`could not fetch nprime.wasm (${response.status})`);
    }
    const bytes = await response.arrayBuffer();
    const module = await WebAssembly.instantiate(bytes, {});
    wasm = module.instance.exports;
    if (wasm.self_test_small() !== 1) {
      throw new Error("WASM self-test failed");
    }
    setStatus("WASM ready", "ready");
    drawSieveWindow();
  } catch (error) {
    setStatus(`WASM load failed: ${error.message}`, "error");
  } finally {
    setBusy(false);
  }
}

for (const button of document.querySelectorAll(".preset")) {
  button.addEventListener("click", () => {
    nInput.value = formatInteger(BigInt(button.dataset.n));
    drawSieveWindow();
  });
}

nInput.addEventListener("change", drawSieveWindow);
window.addEventListener("resize", drawSieveWindow);
runOnceButton.addEventListener("click", computeOnce);
runBenchButton.addEventListener("click", benchmarkIsolated);
traceButton.addEventListener("click", traceCorrection);

drawSieveWindow();
loadWasm();
