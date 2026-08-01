import { readFile, readdir, stat } from "node:fs/promises";
import path from "node:path";

const required = [
  "accession",
  "sample_id",
  "donor_id",
  "donor_id_status",
  "tissue",
  "condition",
  "histology",
  "colon_or_rectum",
  "tumor_location",
  "sporadic_or_FAP",
  "treatment_status",
  "platform",
  "assay_type",
  "data_format",
  "raw_counts_available",
  "clinical_data_available",
  "download_url",
  "inclusion",
  "exclusion_reason",
];

function parseTsv(text) {
  const lines = text.trimEnd().split(/\r?\n/);
  const header = lines[0].split("\t");
  const rows = lines.slice(1).map((line, index) => {
    const cells = line.split("\t");
    if (cells.length !== header.length) {
      throw new Error(`TSV width mismatch at data row ${index + 1}`);
    }
    return Object.fromEntries(header.map((column, i) => [column, cells[i]]));
  });
  return { header, rows };
}

const manifest = parseTsv(
  await readFile(path.resolve("metadata", "dataset_manifest.tsv"), "utf8"),
);
for (const column of required) {
  if (!manifest.header.includes(column)) throw new Error(`Missing required column: ${column}`);
}
for (const [index, row] of manifest.rows.entries()) {
  for (const column of manifest.header) {
    if (row[column] === "") throw new Error(`Blank value row ${index + 1}: ${column}`);
  }
}
if (new Set(manifest.rows.map((row) => row.sample_id)).size !== manifest.rows.length) {
  throw new Error("sample_id values must be globally unique");
}
if (manifest.rows.some((row) => row.donor_id === "NA")) {
  throw new Error("Every biological record must have a donor traceability key");
}
const unresolvedDonors = manifest.rows.filter(
  (row) => row.donor_id_status === "unresolved_sample_scoped_placeholder",
);
if (
  unresolvedDonors.length !== 148 ||
  unresolvedDonors.some(
    (row) => !["GSE41657", "GSE100179"].includes(row.accession),
  )
) {
  throw new Error(`Unexpected unresolved donor mapping: ${unresolvedDonors.length}`);
}

const expectedCounts = {
  GSE201348: 72,
  GSE161277: 13,
  GSE132465: 33,
  GSE41657: 88,
  GSE100179: 60,
  GSE8671: 64,
  GSE99573: 338,
  "TCGA-COAD": 514,
  GSE226997: 4,
};
for (const [accession, expected] of Object.entries(expectedCounts)) {
  const observed = manifest.rows.filter((row) => row.accession === accession).length;
  if (observed !== expected) {
    throw new Error(`${accession}: expected ${expected}, observed ${observed}`);
  }
}

const gse201348 = manifest.rows.filter((row) => row.accession === "GSE201348");
const biological201348 = new Map(
  gse201348.map((row) => [row.biological_sample_id, row]),
);
if (biological201348.size !== 70) throw new Error("GSE201348 must contain 70 biological tissues");
const conditionCounts = {};
for (const row of biological201348.values()) {
  conditionCounts[row.condition] = (conditionCounts[row.condition] ?? 0) + 1;
}
const expected201348 = {
  normal_mucosa: 8,
  unaffected_mucosa: 15,
  polyp_or_adenoma: 42,
  cancer: 5,
};
if (
  Object.entries(expected201348).some(
    ([condition, expected]) => conditionCounts[condition] !== expected,
  )
) {
  throw new Error(`GSE201348 biological counts differ: ${JSON.stringify(conditionCounts)}`);
}
const duplicate201348 = [...new Set(
  gse201348
    .filter(
      (row, index, all) =>
        all.findIndex((candidate) => candidate.biological_sample_id === row.biological_sample_id) !==
        index,
    )
    .map((row) => row.biological_sample_id),
)].sort();
if (duplicate201348.join(",") !== "A002-C-010,A002-C-121") {
  throw new Error(`Unexpected GSE201348 technical replicates: ${duplicate201348.join(",")}`);
}

const gse99573 = manifest.rows.filter((row) => row.accession === "GSE99573");
const splitCounts = Object.fromEntries(
  ["training", "testing", "not_used"].map((split) => [
    split,
    gse99573.filter((row) => row.validation_split === split).length,
  ]),
);
if (splitCounts.training !== 265 || splitCounts.testing !== 65 || splitCounts.not_used !== 8) {
  throw new Error(`GSE99573 split mismatch: ${JSON.stringify(splitCounts)}`);
}
if (
  gse99573
    .filter((row) => row.validation_split === "testing")
    .some((row) => row.inclusion !== "reserve_locked_test_set")
) {
  throw new Error("Every GSE99573 test sample must be locked");
}

const gse161277 = manifest.rows.filter((row) => row.accession === "GSE161277");
const donorCounts161277 = Object.fromEntries(
  ["Patient0", "Patient1", "Patient2", "Patient3"].map((donor) => [
    donor,
    gse161277.filter((row) => row.donor_id === donor).length,
  ]),
);
if (JSON.stringify(donorCounts161277) !== JSON.stringify({
  Patient0: 1,
  Patient1: 3,
  Patient2: 4,
  Patient3: 5,
})) {
  throw new Error(`GSE161277 donor matching mismatch: ${JSON.stringify(donorCounts161277)}`);
}

const prohibited = /\.(fastq|fq|bam|cram|h5|h5ad|rds|rdata|cel|mtx)(\.gz)?$/i;
async function walk(directory) {
  const results = [];
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    if (entry.name === ".git" || entry.name === "node_modules") continue;
    const full = path.join(directory, entry.name);
    if (entry.isDirectory()) results.push(...(await walk(full)));
    else results.push(full);
  }
  return results;
}
for (const file of await walk(path.resolve("."))) {
  if (prohibited.test(file)) throw new Error(`Prohibited large-data type present: ${file}`);
  if ((await stat(file)).size > 10 * 1024 * 1024) {
    throw new Error(`Unexpected file larger than 10 MB: ${file}`);
  }
}

console.log(
  JSON.stringify(
    {
      status: "PASS",
      manifest_rows: manifest.rows.length,
      accession_counts: expectedCounts,
      GSE201348_biological_counts: expected201348,
      GSE161277_donor_counts: donorCounts161277,
      GSE99573_splits: splitCounts,
      prohibited_large_data_files: 0,
      unique_sample_ids: manifest.rows.length,
      donor_keys_present: manifest.rows.length,
      unresolved_sample_scoped_donor_keys: unresolvedDonors.length,
    },
    null,
    2,
  ),
);
