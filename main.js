/* ============================================================
   903 STRANGERS' LIVERS — main.js
   Data loading, charts, scroll behaviour.
   All heavy compute happens in scripts/analyze.py; this file
   just renders pre-computed JSON.
   ============================================================ */

import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

/* ---------- helpers ---------- */
const fmtInt = d3.format(",");
const fmtPct = (n) => `${(n * 100).toFixed(1)}%`;

// Cache-bust JSON fetches with the page-load timestamp so editing and reloading
// the analysis pipeline never serves stale tag distributions from the disk
// cache. (We learned this the hard way after re-running apply_group_tags.py.)
const CACHE_BUST = `?v=${Date.now()}`;
const fetchJSON = (path) =>
  fetch(path + CACHE_BUST).then((r) => {
    if (!r.ok) throw new Error(`${path}: ${r.status}`);
    return r.json();
  });

/** create an SVG with viewBox; returns {svg, width, height, g} */
function makeSvg(container, { width, height, margin = { top: 20, right: 20, bottom: 30, left: 40 } } = {}) {
  const rect = container.getBoundingClientRect();
  const w = width ?? Math.max(480, rect.width - 48);
  const h = height ?? Math.round(w * 0.5);
  const svg = d3.select(container)
    .append("svg")
    .attr("viewBox", `0 0 ${w} ${h}`)
    .attr("width", "100%")
    .attr("height", "auto");
  const g = svg.append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);
  return {
    svg, g, width: w, height: h,
    innerW: w - margin.left - margin.right,
    innerH: h - margin.top - margin.bottom,
    margin,
  };
}

function makeTooltip(className) {
  return d3.select(document.body)
    .append("div")
    .attr("class", className)
    .style("left", "0px")
    .style("top", "0px");
}

/* ============================================================
   1. HERO GRID — 903 points drifting into place
   ============================================================ */
function renderHeroGrid(container, samples) {
  const rect = container.getBoundingClientRect();
  const w = rect.width || window.innerWidth;
  const h = rect.height || window.innerHeight * 0.88;

  const svg = d3.select(container).append("svg")
    .attr("viewBox", `0 0 ${w} ${h}`)
    .attr("preserveAspectRatio", "xMidYMid slice");

  const n = samples.length;
  const cols = Math.ceil(Math.sqrt(n * (w / h)));
  const rows = Math.ceil(n / cols);
  const cellW = w / cols;
  const cellH = h / rows;
  const r = Math.min(cellW, cellH) * 0.18;

  const points = samples.map((gsm, i) => {
    const row = Math.floor(i / cols);
    const col = i % cols;
    return {
      gsm,
      tx: col * cellW + cellW / 2,
      ty: row * cellH + cellH / 2,
      sx: Math.random() * w,
      sy: Math.random() * h,
    };
  });

  svg.selectAll("circle")
    .data(points)
    .join("circle")
    .attr("r", r)
    .attr("fill", "#0C8DC3")
    .attr("fill-opacity", 0.18)
    .attr("cx", (d) => d.sx)
    .attr("cy", (d) => d.sy)
    .transition()
    .delay((_, i) => 200 + (i % 60) * 12)
    .duration(2400)
    .ease(d3.easeCubicInOut)
    .attr("cx", (d) => d.tx)
    .attr("cy", (d) => d.ty)
    .attr("fill-opacity", 0.32);
}

/* ============================================================
   2. FIG-GRID — a clean, static 903-sample grid
   ============================================================ */
function renderGrid(container, samples) {
  const svgDim = makeSvg(container, { width: 960, height: 260, margin: { top: 20, right: 20, bottom: 20, left: 20 } });
  const { g, innerW, innerH } = svgDim;

  const n = samples.length;
  const cols = 43; // sqrt(903) ≈ 30, but wider grid reads better
  const rows = Math.ceil(n / cols);
  const cellW = innerW / cols;
  const cellH = innerH / rows;
  const r = Math.min(cellW, cellH) * 0.36;

  g.selectAll("circle")
    .data(samples)
    .join("circle")
    .attr("cx", (_, i) => (i % cols) * cellW + cellW / 2)
    .attr("cy", (_, i) => Math.floor(i / cols) * cellH + cellH / 2)
    .attr("r", r)
    .attr("fill", "#0C8DC3")
    .attr("fill-opacity", 0.7)
    .append("title")
    .text((d) => d);

  // caption-ish overlay text
  svgDim.svg.append("text")
    .attr("x", svgDim.width - 20)
    .attr("y", svgDim.height - 10)
    .attr("text-anchor", "end")
    .attr("font-family", "var(--sans)")
    .attr("font-size", 11)
    .attr("fill", "rgba(0,33,44,0.5)")
    .text(`903 samples × 35,238 genes`);
}

/* ============================================================
   3. MATRIX HEATMAP — 30 genes × 50 samples slice
   ============================================================ */
function renderMatrix(container, matrix) {
  const genes = matrix.genes;
  const samples = matrix.samples;
  const raw = matrix.raw; // genes × samples

  const svgDim = makeSvg(container, {
    width: 960, height: 520,
    margin: { top: 20, right: 20, bottom: 30, left: 110 },
  });
  const { g, innerW, innerH } = svgDim;

  const x = d3.scaleBand().domain(d3.range(samples.length)).range([0, innerW]).padding(0.05);
  const y = d3.scaleBand().domain(d3.range(genes.length)).range([0, innerH]).padding(0.05);

  // flatten raw for max/scale
  let maxVal = 0;
  for (const row of raw) for (const v of row) if (v > maxVal) maxVal = v;
  const color = d3.scaleSequential()
    .domain([0, Math.log10(maxVal + 1)])
    .interpolator(d3.interpolate("#EAF3F7", "#035C81"));

  const cells = [];
  for (let i = 0; i < genes.length; i++) {
    for (let j = 0; j < samples.length; j++) {
      cells.push({ i, j, v: raw[i][j], gene: genes[i], sample: samples[j] });
    }
  }

  const tooltip = makeTooltip("matrix-tooltip");

  g.selectAll("rect")
    .data(cells)
    .join("rect")
    .attr("class", "matrix-cell")
    .attr("x", (d) => x(d.j))
    .attr("y", (d) => y(d.i))
    .attr("width", x.bandwidth())
    .attr("height", y.bandwidth())
    .attr("fill", (d) => (d.v === 0 ? "#F6FAFC" : color(Math.log10(d.v + 1))))
    .on("mousemove", function (event, d) {
      d3.select(this).attr("stroke", "#CF2FB3").attr("stroke-width", 1.5);
      tooltip.classed("is-visible", true)
        .html(`<strong>${d.gene}</strong> · ${d.sample}<br>raw count: ${fmtInt(d.v)}`)
        .style("left", `${event.pageX + 12}px`)
        .style("top", `${event.pageY - 38}px`);
    })
    .on("mouseleave", function () {
      d3.select(this).attr("stroke", "none");
      tooltip.classed("is-visible", false);
    });

  // gene labels (left)
  g.selectAll(".matrix-label")
    .data(genes)
    .join("text")
    .attr("class", "matrix-label")
    .attr("x", -8)
    .attr("y", (_, i) => y(i) + y.bandwidth() / 2 + 3)
    .attr("text-anchor", "end")
    .attr("font-family", "var(--mono)")
    .attr("font-size", 10)
    .attr("fill", "rgba(0,33,44,0.7)")
    .text((d) => d.length > 14 ? d.slice(0, 13) + "…" : d);
}

/* ============================================================
   4. TOP GENES BAR CHART
   ============================================================ */
function renderTopGenes(container, topGenes) {
  const genes = topGenes.genes;
  const svgDim = makeSvg(container, {
    width: 960, height: 720,
    margin: { top: 10, right: 140, bottom: 30, left: 90 },
  });
  const { g, innerW, innerH } = svgDim;

  const y = d3.scaleBand().domain(genes.map((d) => d.gene)).range([0, innerH]).padding(0.18);
  const x = d3.scaleLinear()
    .domain([0, d3.max(genes, (d) => d.median_logcpm) * 1.05])
    .range([0, innerW]);

  // axis
  const xAxis = d3.axisBottom(x).ticks(6).tickSize(-innerH);
  g.append("g")
    .attr("class", "top-genes-axis")
    .attr("transform", `translate(0, ${innerH})`)
    .call(xAxis)
    .call((g2) => g2.select(".domain").remove())
    .call((g2) => g2.selectAll(".tick line")
      .attr("stroke", "rgba(0,33,44,0.08)"))
    .call((g2) => g2.selectAll(".tick text")
      .attr("font-family", "var(--sans)")
      .attr("font-size", 11)
      .attr("fill", "rgba(0,33,44,0.5)"));

  // x axis label
  g.append("text")
    .attr("x", innerW / 2)
    .attr("y", innerH + 28)
    .attr("text-anchor", "middle")
    .attr("font-family", "var(--sans)")
    .attr("font-size", 11)
    .attr("fill", "rgba(0,33,44,0.6)")
    .text("median log₂(CPM + 1) across 903 samples");

  // IQR range whiskers
  g.selectAll(".iqr-line")
    .data(genes)
    .join("line")
    .attr("class", "iqr-line")
    .attr("x1", (d) => x(d.q25))
    .attr("x2", (d) => x(d.q75))
    .attr("y1", (d) => y(d.gene) + y.bandwidth() / 2)
    .attr("y2", (d) => y(d.gene) + y.bandwidth() / 2)
    .attr("stroke", "rgba(0,33,44,0.22)")
    .attr("stroke-width", 1.5);

  // bars
  const tooltip = makeTooltip("top-genes-tooltip");

  g.selectAll(".top-genes-bar")
    .data(genes)
    .join("rect")
    .attr("class", (d) => `top-genes-bar${d.is_liver_marker ? " is-marker" : ""}`)
    .attr("x", 0)
    .attr("y", (d) => y(d.gene))
    .attr("height", y.bandwidth())
    .attr("width", (d) => x(d.median_logcpm))
    .attr("rx", 1)
    .on("mousemove", function (event, d) {
      const parts = [
        `<div class="gene">${d.gene}</div>`,
        `<div class="stat">median log₂(CPM+1): <strong>${d.median_logcpm.toFixed(2)}</strong></div>`,
        `<div class="stat">IQR: ${d.q25.toFixed(2)} — ${d.q75.toFixed(2)}</div>`,
        `<div class="stat">min → max: ${d.min_logcpm.toFixed(2)} → ${d.max_logcpm.toFixed(2)}</div>`,
      ];
      if (d.annotation) parts.push(`<div class="stat" style="margin-top:6px;color:#47b883">${d.annotation}</div>`);
      tooltip.classed("is-visible", true)
        .html(parts.join(""))
        .style("left", `${event.pageX + 14}px`)
        .style("top", `${event.pageY - 16}px`);
    })
    .on("mouseleave", () => tooltip.classed("is-visible", false));

  // gene labels (right of bar)
  g.selectAll(".top-genes-label")
    .data(genes)
    .join("text")
    .attr("class", (d) => `top-genes-label${d.is_liver_marker ? " is-marker" : ""}`)
    .attr("x", (d) => x(d.median_logcpm) + 6)
    .attr("y", (d) => y(d.gene) + y.bandwidth() / 2 + 3)
    .text((d) => d.gene);

  // legend — positioned at bottom-right of the inner chart area so the sticky
  // site header can never occlude it when the figure is scrolled into view
  const legend = g.append("g")
    .attr("transform", `translate(${innerW - 170}, ${innerH - 48})`);
  legend.append("rect").attr("width", 12).attr("height", 12).attr("fill", "#34A270").attr("rx", 2);
  legend.append("text")
    .attr("x", 18).attr("y", 10)
    .attr("font-family", "var(--sans)").attr("font-size", 11).attr("fill", "rgba(0,33,44,0.7)")
    .text("liver marker");
  legend.append("rect").attr("width", 12).attr("height", 12).attr("fill", "#035C81").attr("rx", 2).attr("y", 20);
  legend.append("text")
    .attr("x", 18).attr("y", 30)
    .attr("font-family", "var(--sans)").attr("font-size", 11).attr("fill", "rgba(0,33,44,0.7)")
    .text("other highly-expressed");
}

/* ============================================================
   5. NORMALIZATION — raw / CPM / log-CPM slider
   ============================================================ */
function renderNormalization(chartContainer, norm) {
  const genes = norm.genes;
  const small = norm.small_sample;
  const large = norm.large_sample;

  // derive CPM on the fly for the middle state
  const withCPM = genes.map((g) => ({
    gene: g.gene,
    raw_small: g.raw_small,
    raw_large: g.raw_large,
    cpm_small: (g.raw_small / small.library_size) * 1e6,
    cpm_large: (g.raw_large / large.library_size) * 1e6,
    logcpm_small: g.logcpm_small,
    logcpm_large: g.logcpm_large,
  }));

  // layout
  const svgDim = makeSvg(chartContainer, {
    width: 960, height: 420,
    margin: { top: 10, right: 20, bottom: 60, left: 110 },
  });
  const { g, innerW, innerH } = svgDim;

  const y = d3.scaleBand()
    .domain(withCPM.map((d) => d.gene))
    .range([0, innerH])
    .padding(0.2);

  // three possible x scales depending on mode
  const scales = {
    raw: d3.scaleLog().domain([1, d3.max(withCPM, (d) => Math.max(d.raw_small + 1, d.raw_large + 1))]).range([0, innerW]).nice(),
    cpm: d3.scaleLog().domain([1, d3.max(withCPM, (d) => Math.max(d.cpm_small + 1, d.cpm_large + 1))]).range([0, innerW]).nice(),
    log: d3.scaleLinear().domain([0, 20]).range([0, innerW]).nice(),
  };
  const valuesFor = {
    raw: (d) => [d.raw_small + 1, d.raw_large + 1],
    cpm: (d) => [d.cpm_small + 1, d.cpm_large + 1],
    log: (d) => [d.logcpm_small, d.logcpm_large],
  };

  const labelFor = {
    raw: "raw read count (log scale)",
    cpm: "counts per million (log scale)",
    log: "log₂(CPM + 1)",
  };

  // gene labels
  g.selectAll(".gene-label")
    .data(withCPM)
    .join("text")
    .attr("class", "norm-gene-label")
    .attr("x", -8)
    .attr("y", (d) => y(d.gene) + y.bandwidth() / 2 + 3)
    .attr("text-anchor", "end")
    .text((d) => d.gene);

  // dumbbell geometry: line + two dots per gene
  const lines = g.selectAll(".dumbbell-line")
    .data(withCPM)
    .join("line")
    .attr("class", "dumbbell-line")
    .attr("y1", (d) => y(d.gene) + y.bandwidth() / 2)
    .attr("y2", (d) => y(d.gene) + y.bandwidth() / 2)
    .attr("stroke", "rgba(0,33,44,0.2)")
    .attr("stroke-width", 1.5);

  const dotSmall = g.selectAll(".dot-small")
    .data(withCPM)
    .join("circle")
    .attr("class", "dot-small")
    .attr("cy", (d) => y(d.gene) + y.bandwidth() / 2)
    .attr("r", 5)
    .attr("fill", "#B7B7B7")
    .attr("stroke", "#fff")
    .attr("stroke-width", 1);

  const dotLarge = g.selectAll(".dot-large")
    .data(withCPM)
    .join("circle")
    .attr("class", "dot-large")
    .attr("cy", (d) => y(d.gene) + y.bandwidth() / 2)
    .attr("r", 5)
    .attr("fill", "#035C81")
    .attr("stroke", "#fff")
    .attr("stroke-width", 1);

  // axis group
  const axisG = g.append("g")
    .attr("class", "norm-axis")
    .attr("transform", `translate(0, ${innerH})`);

  const axisLabel = g.append("text")
    .attr("x", innerW / 2)
    .attr("y", innerH + 44)
    .attr("text-anchor", "middle")
    .attr("font-family", "var(--sans)")
    .attr("font-size", 12)
    .attr("fill", "rgba(0,33,44,0.6)");

  // legend
  const legendWrap = d3.select(chartContainer).insert("div", "svg")
    .attr("class", "norm-legend")
    .html(`
      <span><span class="norm-legend-swatch" style="background:#B7B7B7"></span>${small.gsm} — ${fmtInt(small.library_size)} reads</span>
      <span><span class="norm-legend-swatch" style="background:#035C81"></span>${large.gsm} — ${fmtInt(large.library_size)} reads</span>
    `);

  function update(mode) {
    const x = scales[mode];
    const ticks = mode === "log" ? 6 : 5;
    const axis = mode === "log"
      ? d3.axisBottom(x).ticks(6)
      : d3.axisBottom(x).ticks(ticks, "~s");
    axisG.transition().duration(500).call(axis)
      .call((sel) => sel.selectAll(".tick text")
        .attr("font-family", "var(--sans)")
        .attr("font-size", 11)
        .attr("fill", "rgba(0,33,44,0.55)"))
      .call((sel) => sel.selectAll(".tick line").attr("stroke", "rgba(0,33,44,0.1)"))
      .call((sel) => sel.select(".domain").attr("stroke", "rgba(0,33,44,0.2)"));

    axisLabel.text(labelFor[mode]);

    lines.transition().duration(600).ease(d3.easeCubicInOut)
      .attr("x1", (d) => x(valuesFor[mode](d)[0]))
      .attr("x2", (d) => x(valuesFor[mode](d)[1]));
    dotSmall.transition().duration(600).ease(d3.easeCubicInOut)
      .attr("cx", (d) => x(valuesFor[mode](d)[0]));
    dotLarge.transition().duration(600).ease(d3.easeCubicInOut)
      .attr("cx", (d) => x(valuesFor[mode](d)[1]));
  }

  update("raw");

  // wire up slider
  const slider = document.getElementById("norm-slider");
  const steps = ["raw", "cpm", "log"];
  const stepNodes = document.querySelectorAll(".norm-steps span");
  slider.addEventListener("input", () => {
    const mode = steps[+slider.value];
    update(mode);
    stepNodes.forEach((n, i) => n.classList.toggle("is-active", i === +slider.value));
  });
}

/* ============================================================
   6. PCA scatter with brushing + selection panel
   ============================================================ */
function renderPCA(container, pca, loadings, samplesMeta) {
  const samples = pca.samples;
  const pc1Var = pca.variance_explained[0];
  const pc2Var = pca.variance_explained[1];

  // outer layout: scatter + side panel
  const root = d3.select(container).append("div").attr("class", "pca-layout");
  const scatterWrap = root.append("div");
  const side = root.append("div").attr("class", "pca-side");
  side.html(`
    <h4>Selection</h4>
    <div class="count" data-pca-count>—</div>
    <div class="count-label" data-pca-count-label>Drag to brush a region</div>
    <div data-pca-info style="margin-top:18px"></div>
  `);

  // scatter SVG
  const svgDim = makeSvg(scatterWrap.node(), {
    width: 760, height: 560,
    margin: { top: 24, right: 24, bottom: 50, left: 60 },
  });
  const { svg, g, innerW, innerH } = svgDim;

  const ext1 = d3.extent(samples, (d) => d.pc[0]);
  const ext2 = d3.extent(samples, (d) => d.pc[1]);
  const pad1 = (ext1[1] - ext1[0]) * 0.05;
  const pad2 = (ext2[1] - ext2[0]) * 0.05;

  const x = d3.scaleLinear().domain([ext1[0] - pad1, ext1[1] + pad1]).range([0, innerW]);
  const y = d3.scaleLinear().domain([ext2[0] - pad2, ext2[1] + pad2]).range([innerH, 0]);

  // axes
  g.append("g")
    .attr("class", "pca-axis")
    .attr("transform", `translate(0, ${innerH})`)
    .call(d3.axisBottom(x).ticks(6));
  g.append("g")
    .attr("class", "pca-axis")
    .call(d3.axisLeft(y).ticks(6));

  g.append("text")
    .attr("class", "pca-axis-title")
    .attr("x", innerW / 2)
    .attr("y", innerH + 38)
    .attr("text-anchor", "middle")
    .text(`PC1 — ${fmtPct(pc1Var)} of total variance`);
  g.append("text")
    .attr("class", "pca-axis-title")
    .attr("transform", `translate(-42, ${innerH / 2}) rotate(-90)`)
    .attr("text-anchor", "middle")
    .text(`PC2 — ${fmtPct(pc2Var)} of total variance`);

  // dots
  const dots = g.append("g").attr("class", "pca-dots")
    .selectAll("circle")
    .data(samples)
    .join("circle")
    .attr("class", "pca-dot")
    .attr("cx", (d) => x(d.pc[0]))
    .attr("cy", (d) => y(d.pc[1]))
    .attr("r", 4.2);

  // Native browser <title> tooltip per dot — works alongside the d3.brush
  // overlay (which would otherwise block hover-based tooltips on the dots).
  // Less rich than the Act II hover, but it doesn't fight the brush.
  if (samplesMeta) {
    dots.append("title").text((d) => {
      const meta = samplesMeta[d.gsm] || {};
      const lines = [
        `${d.gsm}    [${meta.disease_tag || "unknown"}]`,
        meta.gse && meta.study_title
          ? `${meta.gse} — ${(meta.study_title || "").slice(0, 100)}`
          : null,
        meta.title ? `Sample: ${(meta.title || "").slice(0, 80)}` : null,
        meta.summary ? `"${(meta.summary || "").slice(0, 120)}"` : null,
      ].filter(Boolean);
      return lines.join("\n");
    });
  } else {
    dots.append("title").text((d) => d.gsm);
  }

  // brush
  const brush = d3.brush()
    .extent([[0, 0], [innerW, innerH]])
    .on("start", () => g.selectAll(".pca-dot").classed("is-faded", false).classed("is-selected", false))
    .on("brush end", ({ selection }) => {
      if (!selection) {
        side.select("[data-pca-count]").text("—");
        side.select("[data-pca-count-label]").text("Drag to brush a region");
        side.select("[data-pca-info]").html("");
        dots.classed("is-selected", false).classed("is-faded", false);
        return;
      }
      const [[x0, y0], [x1, y1]] = selection;
      const selected = [];
      dots.each(function (d) {
        const cx = x(d.pc[0]);
        const cy = y(d.pc[1]);
        const inside = cx >= x0 && cx <= x1 && cy >= y0 && cy <= y1;
        d._selected = inside;
        if (inside) selected.push(d);
      });
      dots.classed("is-selected", (d) => d._selected)
          .classed("is-faded", (d) => !d._selected);

      const n = selected.length;
      side.select("[data-pca-count]").text(n);
      side.select("[data-pca-count-label]")
        .text(`${n === 1 ? "sample" : "samples"} in selection`);

      if (n === 0) {
        side.select("[data-pca-info]").html(`<p class="pca-side-empty">No samples in that region.</p>`);
        return;
      }

      const meanPC1 = d3.mean(selected, (d) => d.pc[0]);
      const meanPC2 = d3.mean(selected, (d) => d.pc[1]);
      const mag = Math.sqrt(meanPC1 * meanPC1 + meanPC2 * meanPC2) || 1;
      const u1 = meanPC1 / mag;
      const u2 = meanPC2 / mag;

      // rank loadings by alignment with selection direction
      const allLoads = []
        .concat(loadings.PC1.positive.map((l) => ({ gene: l.gene, v1: l.loading, v2: 0 })))
        .concat(loadings.PC1.negative.map((l) => ({ gene: l.gene, v1: l.loading, v2: 0 })))
        .concat(loadings.PC2.positive.map((l) => ({ gene: l.gene, v1: 0, v2: l.loading })))
        .concat(loadings.PC2.negative.map((l) => ({ gene: l.gene, v1: 0, v2: l.loading })));

      const merged = new Map();
      allLoads.forEach((l) => {
        const cur = merged.get(l.gene) || { gene: l.gene, v1: 0, v2: 0 };
        cur.v1 += l.v1; cur.v2 += l.v2;
        merged.set(l.gene, cur);
      });
      const scored = [...merged.values()]
        .map((l) => ({ gene: l.gene, score: l.v1 * u1 + l.v2 * u2 }))
        .filter((l) => Math.abs(l.score) > 0)
        .sort((a, b) => b.score - a.score);

      const topK = scored.slice(0, 6);
      side.select("[data-pca-info]").html(`
        <h4 style="margin-top:6px;margin-bottom:6px">Centroid (PC1, PC2)</h4>
        <div style="font-family:var(--mono);font-size:12px;color:rgba(0,33,44,0.7);line-height:1.6">
          PC1 = ${meanPC1.toFixed(2)}<br>PC2 = ${meanPC2.toFixed(2)}
        </div>
        <h4 style="margin-top:18px;margin-bottom:6px">Genes pulling in this direction</h4>
        <ol>
          ${topK.map((t) => `<li><span>${t.gene}</span><span>${t.score >= 0 ? "+" : ""}${t.score.toFixed(3)}</span></li>`).join("")}
        </ol>
      `);
    });

  g.append("g").attr("class", "pca-brush").call(brush);

  return {
    recolor(fn) {
      dots.transition().duration(600).attr("fill", fn);
    },
    resetColor() {
      dots.transition().duration(600).attr("fill", null);
    },
    dots,
  };
}

/* ============================================================
   7. LOADINGS PANEL
   ============================================================ */
function renderLoadings(container, loadings, pcaRef) {
  const root = d3.select(container).append("div").attr("class", "loadings-grid");

  const makeCol = (pcName, data) => {
    const col = root.append("div").attr("class", "loadings-col");
    col.append("h4").html(`Defining <strong>${pcName}</strong> — top positive ↑ / top negative ↓`);

    const all = [...data.positive.map((d) => ({ ...d, sign: "pos" })),
                 ...data.negative.map((d) => ({ ...d, sign: "neg" }))];
    const maxAbs = d3.max(all, (d) => Math.abs(d.loading));

    all.forEach((d) => {
      const row = col.append("div")
        .attr("class", "loading-row")
        .attr("data-gene", d.gene)
        .attr("data-pc", pcName);
      row.append("div").attr("class", "loading-gene").text(d.gene);
      const bar = row.append("div").attr("class", "loading-bar");
      bar.append("div")
        .attr("class", `loading-bar__fill${d.sign === "neg" ? " is-neg" : ""}`)
        .style("width", `${(Math.abs(d.loading) / maxAbs) * 100}%`);
    });
  };

  makeCol("PC1", loadings.PC1);
  makeCol("PC2", loadings.PC2);

  // click handlers just visual for now: highlight
  root.selectAll(".loading-row").on("click", function () {
    const alreadyActive = this.classList.contains("is-selected");
    root.selectAll(".loading-row").classed("is-selected", false);
    if (alreadyActive) return;
    this.classList.add("is-selected");
  });
}

/* ============================================================
   8. META BREAKDOWN (Act II)
   ============================================================ */
function renderMetaBreakdown(diseaseEl, gseEl, stats) {
  function renderBars(container, entries, limit, title) {
    const top = entries.slice(0, limit);
    const total = entries.reduce((s, [, v]) => s + v, 0);
    const maxV = d3.max(top, ([, v]) => v);

    const root = d3.select(container).attr("class", "meta-breakdown");
    root.append("h4").text(title);
    top.forEach(([label, count]) => {
      const row = root.append("div").attr("class", "meta-bar-row");
      row.append("div").attr("class", "meta-bar-label").text(label);
      const track = row.append("div").attr("class", "meta-bar-track");
      track.append("div")
        .attr("class", "meta-bar-fill")
        .style("width", `${(count / maxV) * 100}%`);
      row.append("div").attr("class", "meta-bar-count").text(count);
    });
    if (entries.length > limit) {
      root.append("div")
        .attr("class", "meta-bar-row")
        .style("color", "rgba(0,33,44,0.45)")
        .style("font-style", "italic")
        .html(`<div class="meta-bar-label">+ ${entries.length - limit} more</div><div></div><div></div>`);
    }
  }

  const diseaseEntries = Object.entries(stats.disease_counts).sort((a, b) => b[1] - a[1]);
  const gseEntries = Object.entries(stats.gse_counts)
    .filter(([k]) => k && k !== "unknown")
    .sort((a, b) => b[1] - a[1]);

  renderBars(diseaseEl, diseaseEntries, 8, "Disease tag (inferred)");
  renderBars(gseEl, gseEntries, 12, "Samples per GEO study");
}

/* ============================================================
   9. PROJECTION RECOLOR (Act II) — used for both PCA and UMAP
   ============================================================ */

// Shared palette / coloring helpers, defined once at module scope so the PCA
// and UMAP scatters render with identical colours.
const DISEASE_PALETTE = {
  // primary liver cancers
  "HCC":                                  "#CF2FB3",  // magenta
  "cholangiocarcinoma":                   "#C94E4E",  // red
  "hepatoblastoma":                       "#A02050",  // dark crimson
  // metastasis
  "liver metastasis (non-liver primary)": "#7B1FA2",  // deep purple
  // benign disease
  "NAFLD/NASH":                           "#8C4E9F",  // purple
  "viral hepatitis":                      "#E8803D",  // orange
  "fibrosis/cirrhosis":                   "#D4A13A",  // gold
  "cholestasis":                          "#B66B2C",  // burnt sienna
  // tissue
  "normal adult liver":                   "#34A270",  // emerald
  "normal fetal liver":                   "#7BB662",  // olive green
  "non-tumor adjacent liver":             "#A8D58F",  // pale green
  "other liver":                          "#46A292",  // sea green
  // in-vitro / model systems
  "primary hepatocytes (in vitro)":       "#046865",  // dark teal
  "hepatocyte cell line":                 "#1B998B",  // teal
  "hepatic organoid":                     "#F59E0B",  // amber
  "iPSC/ESC-derived hepatocyte-like":     "#3B82F6",  // bright blue
  // sorted cell types
  "hepatic stellate cell":                "#6A4C93",  // deep purple
  "liver-resident immune cell":           "#5B7CD8",  // periwinkle
  "liver-derived hematopoietic":          "#88B4DC",  // pale blue
  // catch-alls
  "non-liver":                            "#9AA5B1",  // desaturated gray-blue
  "unknown":                              "#B7B7B7",  // gray
};

const GSE_PALETTE = [
  "#CF2FB3", "#F59E0B", "#34A270", "#3B82F6", "#8C4E9F",
  "#E8803D", "#1B998B", "#C94E4E", "#6A4C93", "#D4A13A",
  "#0C8DC3", "#046865",
];

function buildGseColorFn(samples, samplesMeta) {
  const counts = {};
  for (const s of samples) {
    const gse = (samplesMeta[s.gsm] || {}).gse;
    if (gse) counts[gse] = (counts[gse] || 0) + 1;
  }
  const topGSEs = Object.entries(counts)
    .filter(([, n]) => n >= 10)
    .sort((a, b) => b[1] - a[1])
    .map(([g]) => g);
  const scale = d3.scaleOrdinal().domain(topGSEs).range(GSE_PALETTE);
  const fn = (d) => {
    const gse = (samplesMeta[d.gsm] || {}).gse;
    if (!gse) return "#D8DEE4";
    return counts[gse] >= 10 ? scale(gse) : "#D8DEE4";
  };
  fn._counts = counts;
  fn._top = topGSEs;
  return fn;
}

// Warn once per unknown tag so mismatches between samples_meta.json and the
// DISEASE_PALETTE surface as a console error instead of a silent gray cluster.
const _warnedTags = new Set();
function diseaseColorFn(samplesMeta) {
  return (d) => {
    const tag = (samplesMeta[d.gsm] || {}).disease_tag || "unknown";
    if (!(tag in DISEASE_PALETTE) && !_warnedTags.has(tag)) {
      _warnedTags.add(tag);
      console.warn(
        `[disease palette] tag "${tag}" not in DISEASE_PALETTE — falling back to gray. ` +
        `Either samples_meta.json is stale or main.js needs a new color for this tag.`
      );
    }
    return DISEASE_PALETTE[tag] || "#B7B7B7";
  };
}

/**
 * Render a 2D projection (PCA or UMAP) with three recolor modes.
 *
 * @param container         The figure DOM element to render into.
 * @param samples           Array of sample records, each with { gsm, ... }.
 * @param accessor          (sample) -> [x, y] in projection coordinates.
 * @param samplesMeta       Per-sample GEO metadata, keyed by GSM.
 * @param controlsSelector  CSS selector for the recolor button group.
 * @param axisLabels        { x, y } strings for axis titles.
 */
function renderProjectionRecolor({
  container,
  chartEl,
  samples,
  accessor,
  samplesMeta,
  controlsSelector,
  axisLabels,
}) {
  const svgDim = makeSvg(chartEl, {
    width: 960, height: 560,
    margin: { top: 24, right: 24, bottom: 50, left: 60 },
  });
  const { svg, g, innerW, innerH, margin } = svgDim;

  // axis domain = full data extent (with a hair of padding)
  const ext1 = d3.extent(samples, (d) => accessor(d)[0]);
  const ext2 = d3.extent(samples, (d) => accessor(d)[1]);
  const pad1 = (ext1[1] - ext1[0]) * 0.05;
  const pad2 = (ext2[1] - ext2[0]) * 0.05;
  const x = d3.scaleLinear().domain([ext1[0] - pad1, ext1[1] + pad1]).range([0, innerW]);
  const y = d3.scaleLinear().domain([ext2[0] - pad2, ext2[1] + pad2]).range([innerH, 0]);

  // axis groups (kept around so the zoom handler can rewrite them)
  const gxAxis = g.append("g").attr("class", "pca-axis").attr("transform", `translate(0, ${innerH})`).call(d3.axisBottom(x).ticks(6));
  const gyAxis = g.append("g").attr("class", "pca-axis").call(d3.axisLeft(y).ticks(6));

  g.append("text").attr("class", "pca-axis-title")
    .attr("x", innerW / 2).attr("y", innerH + 38).attr("text-anchor", "middle")
    .text(axisLabels.x);
  g.append("text").attr("class", "pca-axis-title")
    .attr("transform", `translate(-42, ${innerH / 2}) rotate(-90)`)
    .attr("text-anchor", "middle")
    .text(axisLabels.y);

  // clip-path so panned dots don't draw over the axis labels
  const clipId = `clip-${Math.random().toString(36).slice(2)}`;
  svg.append("defs").append("clipPath")
    .attr("id", clipId)
    .append("rect")
    .attr("width", innerW)
    .attr("height", innerH);

  const dotsLayer = g.append("g").attr("clip-path", `url(#${clipId})`);
  const dots = dotsLayer.selectAll("circle")
    .data(samples)
    .join("circle")
    .attr("class", "pca-dot")
    .attr("cx", (d) => x(accessor(d)[0]))
    .attr("cy", (d) => y(accessor(d)[1]))
    .attr("r", 4.2);

  // ---- hover halo (drawn above dots, in its own clipped layer) -------------
  const haloLayer = g.append("g").attr("clip-path", `url(#${clipId})`);
  const halo = haloLayer.append("circle")
    .attr("class", "projection-halo")
    .attr("r", 9)
    .attr("fill", "none")
    .attr("stroke", "#00212C")
    .attr("stroke-width", 2)
    .attr("pointer-events", "none")
    .style("opacity", 0);

  // ---- pan + zoom -----------------------------------------------------------
  // Computes the d3.zoomIdentity transform that fits a percentile-clipped
  // bounding box into the chart, with 8% padding.
  //
  // With UMAP parameters tuned to produce a roughly square, compact embedding,
  // the default [0, 100] fits all data. If a dataset grows back an outlier
  // cluster, callers can pass e.g. [5, 95] to clip it out of the initial view.
  function fitPercentileTransform(p = [0, 100], padFrac = 0.05) {
    const xs = samples.map((d) => accessor(d)[0]).sort(d3.ascending);
    const ys = samples.map((d) => accessor(d)[1]).sort(d3.ascending);
    const x0 = d3.quantile(xs, p[0] / 100);
    const x1 = d3.quantile(xs, p[1] / 100);
    const y0 = d3.quantile(ys, p[0] / 100);
    const y1 = d3.quantile(ys, p[1] / 100);
    const px0 = x(x0), px1 = x(x1);
    const py0 = y(y0), py1 = y(y1); // y is inverted
    const boxW = Math.abs(px1 - px0);
    const boxH = Math.abs(py1 - py0);
    const k = Math.max(0.5, Math.min(20, (1 - padFrac) * Math.min(innerW / boxW, innerH / boxH)));
    const cx = (px0 + px1) / 2;
    const cy = (py0 + py1) / 2;
    const tx = innerW / 2 - k * cx;
    const ty = innerH / 2 - k * cy;
    return d3.zoomIdentity.translate(tx, ty).scale(k);
  }

  let currentTransform = d3.zoomIdentity;

  // Quadtree built in screen space, used for nearest-neighbour hit-testing
  // on hover. Rebuilt every time the transform changes so the lookup is
  // always current with the on-screen positions of the dots.
  let quadtree = null;
  function rebuildQuadtree() {
    quadtree = d3.quadtree()
      .x((d) => currentTransform.applyX(x(accessor(d)[0])))
      .y((d) => currentTransform.applyY(y(accessor(d)[1])))
      .addAll(samples);
  }

  function applyTransform(t) {
    currentTransform = t;
    const newX = t.rescaleX(x);
    const newY = t.rescaleY(y);
    gxAxis.call(d3.axisBottom(newX).ticks(6))
      .call((sel) => sel.selectAll("text")
        .attr("font-family", "var(--sans)").attr("font-size", 11).attr("fill", "rgba(0,33,44,0.55)"))
      .call((sel) => sel.selectAll(".tick line").attr("stroke", "rgba(0,33,44,0.1)"))
      .call((sel) => sel.select(".domain").attr("stroke", "rgba(0,33,44,0.2)"));
    gyAxis.call(d3.axisLeft(newY).ticks(6))
      .call((sel) => sel.selectAll("text")
        .attr("font-family", "var(--sans)").attr("font-size", 11).attr("fill", "rgba(0,33,44,0.55)"))
      .call((sel) => sel.selectAll(".tick line").attr("stroke", "rgba(0,33,44,0.1)"))
      .call((sel) => sel.select(".domain").attr("stroke", "rgba(0,33,44,0.2)"));
    dots.attr("cx", (d) => newX(accessor(d)[0]))
        .attr("cy", (d) => newY(accessor(d)[1]));
    rebuildQuadtree();
  }

  // ---- tooltip --------------------------------------------------------------
  // Declared BEFORE the zoom handlers because the zoom "start" handler
  // calls hideTooltip(), which references `tooltip`. If tooltip is in TDZ
  // when the initial fitPercentileTransform fires its zoom event, the start
  // handler throws and the entire zoom-init chain (including
  // rebuildQuadtree) silently fails — leaving the chart with no hover data.
  const tooltip = makeTooltip("projection-tooltip");

  let isZooming = false;
  function hideTooltip() {
    tooltip.classed("is-visible", false);
    halo.style("opacity", 0);
  }

  const zoom = d3.zoom()
    .scaleExtent([0.5, 20])
    .on("start", () => { isZooming = true; hideTooltip(); })
    .on("zoom", (event) => applyTransform(event.transform))
    .on("end", () => { isZooming = false; });

  function escHTML(s) {
    return String(s == null ? "" : s).replace(/[&<>"']/g, (c) => ({
      "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
    }[c]));
  }

  function tooltipHTML(sample) {
    const meta = samplesMeta[sample.gsm] || {};
    const tag = meta.disease_tag || "unknown";
    const color = DISEASE_PALETTE[tag] || "#B7B7B7";
    const studyTitle = (meta.study_title || "").slice(0, 110);
    const sampleTitle = (meta.title || "").slice(0, 80);
    const sampleSummary = (meta.summary || "").slice(0, 130);
    const rationale = (meta.llm_rationale || "").slice(0, 180);
    const conf = meta.llm_confidence;
    const isLiver = meta.llm_is_liver_focused;
    return `
      <div class="projection-tooltip__row">
        <span class="projection-tooltip__chip" style="background:${color}">${escHTML(tag)}</span>
        <span class="projection-tooltip__gsm">${escHTML(sample.gsm)}</span>
      </div>
      ${meta.gse ? `<div class="projection-tooltip__study">${escHTML(meta.gse)}${studyTitle ? " — " + escHTML(studyTitle) : ""}</div>` : ""}
      ${sampleTitle ? `<div class="projection-tooltip__title">${escHTML(sampleTitle)}</div>` : ""}
      ${sampleSummary ? `<div class="projection-tooltip__summary">&ldquo;${escHTML(sampleSummary)}&rdquo;</div>` : ""}
      ${rationale ? `<div class="projection-tooltip__rationale">${escHTML(rationale)}${conf ? ` <em>(${escHTML(conf)} conf.)</em>` : ""}</div>` : ""}
      ${isLiver === false ? `<div class="projection-tooltip__warn">⚠ flagged not liver-focused</div>` : ""}
    `;
  }

  // Invisible rect on top of dots to receive pan/zoom gestures AND
  // mousemove for nearest-neighbour hover hit-testing.
  const catcher = g.append("rect")
    .attr("class", "zoom-catcher")
    .attr("width", innerW)
    .attr("height", innerH)
    .attr("fill", "transparent")
    .style("cursor", "grab")
    .call(zoom)
    .on("dblclick.zoom", null)
    .on("dblclick", () => {
      g.select(".zoom-catcher").transition().duration(600).call(zoom.transform, d3.zoomIdentity);
    });

  catcher
    .on("mousemove", function (event) {
      if (isZooming || !quadtree) return;
      const [mx, my] = d3.pointer(event, this);
      const found = quadtree.find(mx, my, 14);
      if (found) {
        const px = currentTransform.applyX(x(accessor(found)[0]));
        const py = currentTransform.applyY(y(accessor(found)[1]));
        halo.attr("cx", px).attr("cy", py).style("opacity", 1);
        tooltip
          .classed("is-visible", true)
          .html(tooltipHTML(found))
          .style("left", `${event.pageX + 14}px`)
          .style("top", `${event.pageY - 10}px`);
        this.style.cursor = "pointer";
      } else {
        hideTooltip();
        this.style.cursor = "grab";
      }
    })
    .on("mouseleave", hideTooltip);

  // Apply the initial fit-to-densest-region zoom.
  // Wrapped in a microtask so the call stack unwinds before .call(zoom.transform)
  // fires its synchronous zoom event handler.
  queueMicrotask(() => {
    g.select(".zoom-catcher").call(zoom.transform, fitPercentileTransform());
  });

  // small instructional caption inside the chart, top-left corner
  svg.append("text")
    .attr("class", "zoom-hint")
    .attr("x", margin.left + 8)
    .attr("y", margin.top + 14)
    .attr("font-family", "var(--sans)")
    .attr("font-size", 10.5)
    .attr("fill", "rgba(0,33,44,0.45)")
    .text("scroll to zoom · drag to pan · double-click to fit all");

  const legend = svgDim.svg.append("g")
    .attr("class", "recolor-legend")
    .attr("transform", `translate(${svgDim.width - 220}, 30)`);

  function clearLegend() { legend.selectAll("*").remove(); }

  function drawLegend(entries) {
    clearLegend();
    entries.forEach(([label, color], i) => {
      const row = legend.append("g").attr("transform", `translate(0, ${i * 18})`);
      row.append("rect").attr("width", 12).attr("height", 12).attr("fill", color).attr("rx", 2);
      row.append("text")
        .attr("x", 18).attr("y", 10)
        .attr("font-family", "var(--sans)").attr("font-size", 11).attr("fill", "rgba(0,33,44,0.7)")
        .text(label.length > 24 ? label.slice(0, 23) + "…" : label);
    });
  }

  function apply(mode) {
    if (mode === "none") {
      dots.transition().duration(500).style("fill", null).style("fill-opacity", null);
      clearLegend();
      return;
    }
    if (mode === "disease") {
      const fn = diseaseColorFn(samplesMeta);
      dots.transition().duration(700).style("fill", fn).style("fill-opacity", 0.88);
      const counts = {};
      for (const s of samples) {
        const t = (samplesMeta[s.gsm] || {}).disease_tag || "unknown";
        counts[t] = (counts[t] || 0) + 1;
      }
      // Cap legend to top 12 categories to keep it readable; bucket the rest
      const sorted = Object.entries(counts).sort((a, b) => b[1] - a[1]);
      const top = sorted.slice(0, 12);
      const rest = sorted.slice(12);
      const entries = top.map(([k, v]) => [`${k} (${v})`, DISEASE_PALETTE[k] || "#B7B7B7"]);
      if (rest.length) {
        const restTotal = rest.reduce((s, [, n]) => s + n, 0);
        entries.push([`+ ${rest.length} smaller categories (${restTotal})`, "#D8DEE4"]);
      }
      drawLegend(entries);
      return;
    }
    if (mode === "gse") {
      const fn = buildGseColorFn(samples, samplesMeta);
      dots.transition().duration(700).style("fill", fn).style("fill-opacity", 0.88);
      const counts = fn._counts;
      const top = Object.entries(counts)
        .filter(([, n]) => n >= 10)
        .sort((a, b) => b[1] - a[1]);
      const entries = top.map(([gse, count]) => {
        const fakeSample = samples.find((s) => (samplesMeta[s.gsm] || {}).gse === gse);
        return [`${gse} (${count})`, fn(fakeSample)];
      });
      const smallTotal = Object.entries(counts)
        .filter(([, n]) => n < 10)
        .reduce((s, [, n]) => s + n, 0);
      if (smallTotal) entries.push([`studies < 10 (${smallTotal})`, "#D8DEE4"]);
      drawLegend(entries);
      return;
    }
  }

  // button wiring — scoped to this projection's controls only
  const buttons = container.querySelectorAll(`${controlsSelector} button`);
  buttons.forEach((btn) => {
    btn.addEventListener("click", () => {
      buttons.forEach((b) => b.classList.remove("is-active"));
      btn.classList.add("is-active");
      apply(btn.dataset.recolor);
    });
  });
}

function renderPCArecolor(container, pca, samplesMeta) {
  renderProjectionRecolor({
    container,
    chartEl: document.getElementById("fig-pca-recolor-chart"),
    samples: pca.samples,
    accessor: (d) => d.pc,
    samplesMeta,
    controlsSelector: ".recolor-controls",
    axisLabels: {
      x: `PC1 — ${fmtPct(pca.variance_explained[0])}`,
      y: `PC2 — ${fmtPct(pca.variance_explained[1])}`,
    },
  });
}

function renderUMAPrecolor(container, umap, samplesMeta) {
  renderProjectionRecolor({
    container,
    chartEl: document.getElementById("fig-umap-chart"),
    samples: umap.samples,
    accessor: (d) => d.u,
    samplesMeta,
    controlsSelector: ".recolor-controls--umap",
    axisLabels: {
      x: "UMAP-1",
      y: "UMAP-2",
    },
  });
}

/* ============================================================
   10. STAT SLOTS & scroll observer
   ============================================================ */
function fillStats(meta, pca) {
  const setStat = (key, value) => {
    document.querySelectorAll(`[data-stat="${key}"]`).forEach((el) => { el.textContent = value; });
  };
  setStat("total-cells", fmtInt(meta.total_cells));
  setStat("zero-frac", fmtPct(meta.zero_fraction));
  setStat("pc1-var", fmtPct(pca.variance_explained[0]));
  setStat("pc2-var", fmtPct(pca.variance_explained[1]));
}

function fillMetaStats(stats) {
  document.querySelectorAll(`[data-stat="n-gse"]`).forEach((el) => {
    el.textContent = `${stats.n_gse}`;
  });
}

function initScrollObserver() {
  const io = new IntersectionObserver(
    (entries) => {
      entries.forEach((e) => {
        if (e.isIntersecting) {
          e.target.classList.add("is-in-view");
          io.unobserve(e.target);
        }
      });
    },
    { threshold: 0.15, rootMargin: "0px 0px -10% 0px" }
  );
  document.querySelectorAll("[data-fig], .prose, .paragraph, .section-title, .colophon").forEach((el) => io.observe(el));
}

/* ============================================================
   MAIN
   ============================================================ */
async function main() {
  // ?visible or ?reveal query param bypasses scroll-triggered hiding
  // (useful for full-page screenshots and printing).
  const params = new URLSearchParams(location.search);
  const revealAll = params.has("visible") || params.has("reveal");
  if (revealAll) document.body.classList.add("reveal-all");

  // Try to load metadata in the same batch as the rest. If it's not available
  // (older snapshot, fetch script not yet run), the projection charts will
  // still render but without per-sample tooltips and recoloring.
  let samplesMeta = null;
  let metaStats = null;
  const dataPromises = [
    fetchJSON("data/meta.json"),
    fetchJSON("data/samples.json"),
    fetchJSON("data/top_genes.json"),
    fetchJSON("data/normalization_demo.json"),
    fetchJSON("data/pca.json"),
    fetchJSON("data/umap.json"),
    fetchJSON("data/loadings.json"),
    fetchJSON("data/matrix_sample.json"),
  ];
  const metaPromises = [
    fetchJSON("data/meta_stats.json").then((d) => { metaStats = d; }).catch(() => {}),
    fetchJSON("data/samples_meta.json").then((d) => { samplesMeta = d; }).catch(() => {}),
  ];
  const [meta, samplesDoc, topGenes, norm, pca, umap, loadings, matrixSample] = await Promise.all(dataPromises);
  await Promise.all(metaPromises);

  fillStats(meta, pca);

  // Act I
  renderHeroGrid(document.getElementById("hero-grid"), samplesDoc.samples);
  renderGrid(document.getElementById("fig-grid"), samplesDoc.samples);
  renderMatrix(document.getElementById("fig-matrix"), matrixSample);
  renderTopGenes(document.getElementById("fig-top-genes"), topGenes);
  renderNormalization(document.getElementById("fig-norm-chart"), norm);

  const pcaRef = renderPCA(document.getElementById("fig-pca"), pca, loadings, samplesMeta);
  renderLoadings(document.getElementById("fig-loadings"), loadings, pcaRef);

  // Act II — render only if metadata loaded successfully
  if (samplesMeta && metaStats) {
    fillMetaStats(metaStats);
    renderMetaBreakdown(
      document.getElementById("meta-disease"),
      document.getElementById("meta-gse"),
      metaStats,
    );
    renderPCArecolor(document.getElementById("fig-pca-recolor"), pca, samplesMeta);
    renderUMAPrecolor(document.getElementById("fig-umap"), umap, samplesMeta);
  } else {
    console.warn("Metadata not available — Act II charts disabled");
    const warn = (text) => `<p class="pca-side-empty" style="padding:24px;text-align:center">${text}</p>`;
    document.getElementById("fig-metabreakdown").innerHTML = warn(
      "Metadata fetch is still running (or was skipped). Run <code>scripts/fetch_metadata.py</code> and reload."
    );
    document.getElementById("fig-pca-recolor-chart").innerHTML = warn(
      "Metadata-based recolouring will appear here once the fetch finishes."
    );
    const umapChart = document.getElementById("fig-umap-chart");
    if (umapChart) umapChart.innerHTML = warn("UMAP recolouring requires metadata.");
  }

  if (!revealAll) initScrollObserver();
  // Signal to automated consumers (screenshot tooling) that rendering is done.
  document.body.dataset.ready = "true";
}

main().catch((err) => {
  console.error(err);
  document.body.insertAdjacentHTML("beforeend",
    `<pre style="padding:2em;background:#fee;color:#900;font-family:monospace">${err.stack || err}</pre>`
  );
});
