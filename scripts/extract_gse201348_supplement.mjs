import { mkdir, writeFile } from "node:fs/promises";
import path from "node:path";
import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const inputPath = process.argv[2];
if (!inputPath) {
  throw new Error("Usage: node extract_gse201348_supplement.mjs <supplement.xlsx>");
}

const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(inputPath));
const inspected = await workbook.inspect({
  kind: "table",
  maxChars: 200_000,
  tableMaxRows: 120,
  tableMaxCols: 30,
  tableMaxCellChars: 500,
});

const tables = {};
for (const line of inspected.ndjson.split(/\r?\n/)) {
  if (!line.trim()) continue;
  const record = JSON.parse(line);
  if (record.kind === "table") tables[record.sheet] = record.values;
}

const outDir = path.resolve("metadata", "source_cache");
await mkdir(outDir, { recursive: true });
await writeFile(
  path.join(outDir, "GSE201348_supplement_tables.json"),
  `${JSON.stringify(
    {
      source_file: path.basename(inputPath),
      source_url:
        "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01088-x/MediaObjects/41588_2022_1088_MOESM3_ESM.xlsx",
      tables,
    },
    null,
    2,
  )}\n`,
  "utf8",
);
