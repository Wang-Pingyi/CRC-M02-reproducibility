import { execFile } from "node:child_process";
import { readFile, writeFile } from "node:fs/promises";
import path from "node:path";
import { promisify } from "node:util";

const execFileAsync = promisify(execFile);
const cacheDir = path.resolve("metadata", "source_cache");
const rows = [];

const cell = (value) =>
  String(value ?? "NA").replaceAll("\t", " ").replaceAll("\n", " ");
const fileNameFromUrl = (url) => {
  const parsed = new URL(url);
  const explicit = parsed.searchParams.get("file");
  if (explicit) return explicit;
  return path.basename(parsed.pathname);
};

function add({
  accession,
  fileName,
  url,
  dataClass,
  priority = "core",
  version = "GEO_current_public_record",
}) {
  rows.push({
    accession,
    relative_path: `${accession}/${fileName}`,
    official_url: url,
    data_class: dataClass,
    priority,
    source_version: version,
  });
}

async function curlText(url) {
  const { stdout } = await execFileAsync(
    "curl.exe",
    [
      "-sS",
      "-L",
      "--fail",
      "--retry",
      "4",
      "--retry-all-errors",
      "--connect-timeout",
      "30",
      "--max-time",
      "120",
      "-A",
      "CRC-carcinogenesis-stage4B/1.0",
      url,
    ],
    { encoding: "utf8", maxBuffer: 8 * 1024 * 1024 },
  );
  return stdout;
}

async function addGeoSampleSupplementary(accession) {
  const source = JSON.parse(
    await readFile(path.join(cacheDir, `${accession}_samples.json`), "utf8"),
  );
  for (const sample of source.samples) {
    const page = await curlText(
      `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${sample.sample_id}`,
    );
    const links = [
      ...page.matchAll(/href="([^"]*\/geo\/download\/\?[^"]+)"/g),
    ].map((match) => match[1].replaceAll("&amp;", "&"));
    for (const relative of [...new Set(links)]) {
      const url = new URL(relative, "https://www.ncbi.nlm.nih.gov").toString();
      const fileName = fileNameFromUrl(url);
      const dataClass = /\.(mtx|tsv)\.gz$/i.test(fileName)
        ? "processed_count_or_feature_matrix"
        : accession === "GSE226997"
          ? "processed_spatial_bundle"
          : "processed_supplementary";
      add({
        accession,
        fileName,
        url,
        dataClass,
        priority: accession === "GSE226997" ? "optional_spatial" : "core",
      });
    }
  }
}

for (const accession of [
  "GSE201348",
  "GSE161277",
  "GSE132465",
  "GSE41657",
  "GSE100179",
  "GSE8671",
  "GSE99573",
  "GSE226997",
]) {
  add({
    accession,
    fileName: `${accession}_series_metadata.soft.txt`,
    url:
      `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${accession}` +
      "&targ=self&form=text&view=quick",
    dataClass: "series_metadata",
  });
}

await addGeoSampleSupplementary("GSE201348");
await addGeoSampleSupplementary("GSE161277");
await addGeoSampleSupplementary("GSE226997");

for (const fileName of [
  "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz",
  "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz",
]) {
  add({
    accession: "GSE132465",
    fileName,
    url:
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132465&format=file&file=" +
      encodeURIComponent(fileName),
    dataClass: fileName.includes("annotation")
      ? "cell_metadata"
      : "processed_raw_UMI_count_matrix",
  });
}

for (const accession of ["GSE41657", "GSE100179", "GSE8671", "GSE99573"]) {
  add({
    accession,
    fileName: `${accession}_RAW.tar`,
    url: `https://www.ncbi.nlm.nih.gov/geo/download/?acc=${accession}&format=file`,
    dataClass:
      accession === "GSE41657" ? "raw_feature_extraction_TXT" : "raw_CEL_archive",
  });
}

const gse100179Processed =
  "GSE100179_2017-07-21-annotated-RMA-SKETCH.RMA-GENE-FULL-Group1.txt.gz";
add({
  accession: "GSE100179",
  fileName: gse100179Processed,
  url:
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE100179&format=file&file=" +
    encodeURIComponent(gse100179Processed),
  dataClass: "processed_expression_matrix",
});

const tcga = JSON.parse(
  await readFile(path.join(cacheDir, "TCGA-COAD_rnaseq_files.json"), "utf8"),
);
for (const file of tcga.data.hits) {
  add({
    accession: "TCGA-COAD",
    fileName: file.file_name,
    url: `https://api.gdc.cancer.gov/data/${file.file_id}`,
    dataClass: "STAR_gene_counts",
    version: `GDC_file_uuid:${file.file_id}`,
  });
}

rows.sort((a, b) =>
  `${a.priority}\t${a.accession}\t${a.relative_path}`.localeCompare(
    `${b.priority}\t${b.accession}\t${b.relative_path}`,
  ),
);
const duplicates = rows.filter(
  (row, index) =>
    rows.findIndex((candidate) => candidate.relative_path === row.relative_path) !==
    index,
);
if (duplicates.length) {
  throw new Error(
    `Duplicate relative paths: ${duplicates.map((row) => row.relative_path).join(", ")}`,
  );
}

const columns = [
  "accession",
  "relative_path",
  "official_url",
  "data_class",
  "priority",
  "source_version",
];
await writeFile(
  path.resolve("metadata", "download_plan.tsv"),
  [
    columns.join("\t"),
    ...rows.map((row) => columns.map((column) => cell(row[column])).join("\t")),
    "",
  ].join("\n"),
  "utf8",
);

console.log(
  JSON.stringify(
    {
      files: rows.length,
      core: rows.filter((row) => row.priority === "core").length,
      optional_spatial: rows.filter((row) => row.priority === "optional_spatial")
        .length,
      fastq_files: rows.filter((row) => /fastq|\.fq/i.test(row.relative_path)).length,
    },
    null,
    2,
  ),
);
