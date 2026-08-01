import { mkdir, writeFile } from "node:fs/promises";
import { execFile } from "node:child_process";
import path from "node:path";
import { promisify } from "node:util";

const accessions = [
  "GSE201348",
  "GSE161277",
  "GSE132465",
  "GSE41657",
  "GSE100179",
  "GSE8671",
  "GSE99573",
  "GSE226997",
];

const outDir = path.resolve("metadata", "source_cache");
await mkdir(outDir, { recursive: true });

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const execFileAsync = promisify(execFile);

async function fetchText(url, attempts = 4) {
  let lastError;
  for (let attempt = 1; attempt <= attempts; attempt += 1) {
    try {
      const { stdout } = await execFileAsync("curl.exe", [
        "-sS",
        "-L",
        "--max-time",
        "30",
        "--retry",
        "2",
        "--retry-delay",
        "1",
        "-A",
        "CRC-carcinogenesis-metadata-audit/1.0",
        url,
      ], {
        encoding: "utf8",
        maxBuffer: 8 * 1024 * 1024,
      });
      return stdout;
    } catch (error) {
      lastError = error;
      await sleep(attempt * 750);
    }
  }
  throw lastError;
}

function parseSeriesSamples(html) {
  const pattern =
    /acc=(GSM\d+)[^>]*>\1<\/a><\/td>\s*<td[^>]*>([^<]+)<\/td>/g;
  return [...html.matchAll(pattern)].map((match) => ({
    sample_id: match[1],
    title: match[2].trim(),
  }));
}

function parseSoft(text) {
  const result = { characteristics: {}, raw_lines: [] };
  for (const line of text.split(/\r?\n/)) {
    if (!line.startsWith("!Sample_") && !line.startsWith("^SAMPLE")) continue;
    result.raw_lines.push(line);
    const pair = line.match(/^!Sample_([^=]+?)\s*=\s*(.*)$/);
    if (!pair) {
      if (line.startsWith("^SAMPLE")) {
        result.sample_id = line.split("=").slice(1).join("=").trim();
      }
      continue;
    }
    const key = pair[1].trim();
    const value = pair[2].trim();
    if (key.startsWith("characteristics")) {
      const characteristic = value.match(/^([^:]+):\s*(.*)$/);
      if (characteristic) {
        result.characteristics[characteristic[1].trim()] =
          characteristic[2].trim();
      } else {
        result.characteristics[`unparsed_${Object.keys(result.characteristics).length + 1}`] =
          value;
      }
    } else if (result[key] === undefined) {
      result[key] = value;
    } else if (Array.isArray(result[key])) {
      result[key].push(value);
    } else {
      result[key] = [result[key], value];
    }
  }
  return result;
}

async function mapConcurrent(items, concurrency, worker) {
  const output = new Array(items.length);
  let cursor = 0;
  async function run() {
    while (cursor < items.length) {
      const index = cursor;
      cursor += 1;
      output[index] = await worker(items[index], index);
    }
  }
  await Promise.all(Array.from({ length: concurrency }, run));
  return output;
}

for (const accession of accessions) {
  const seriesUrl = `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${accession}`;
  const html = await fetchText(seriesUrl);
  const listed = parseSeriesSamples(html);
  if (listed.length === 0) {
    throw new Error(`${accession}: no GSM records parsed`);
  }
  console.log(`${accession}: collecting ${listed.length} GSM metadata records`);
  const samples = await mapConcurrent(listed, 10, async (sample) => {
    const url =
      `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${sample.sample_id}` +
      "&targ=self&form=text&view=quick";
    const soft = parseSoft(await fetchText(url));
    return { ...sample, source_url: url, ...soft };
  });
  await writeFile(
    path.join(outDir, `${accession}_samples.json`),
    `${JSON.stringify({ accession, series_url: seriesUrl, samples }, null, 2)}\n`,
    "utf8",
  );
}
