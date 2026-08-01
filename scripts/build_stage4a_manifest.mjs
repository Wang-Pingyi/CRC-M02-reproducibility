import { readFile, writeFile } from "node:fs/promises";
import path from "node:path";

const cacheDir = path.resolve("metadata", "source_cache");
const outputDir = path.resolve("metadata");
const columns = [
  "accession",
  "sample_id",
  "donor_id",
  "donor_id_status",
  "biological_sample_id",
  "sample_role",
  "technical_replicate_of",
  "paired_group",
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
  "validation_split",
  "download_url",
  "inclusion",
  "exclusion_reason",
  "source_evidence",
  "metadata_notes",
];

const NA = "NA";
const rows = [];
const clean = (value) =>
  value === undefined || value === null || value === "" ? NA : String(value);
const tsvCell = (value) => clean(value).replaceAll("\t", " ").replaceAll("\n", " ");
const regionClass = (location) => {
  if (!location || location === NA) return NA;
  return /rectum/i.test(location) ? "rectum" : "colon";
};
const sourceUrl = (gsm) =>
  `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm}`;

function add(row) {
  if (row.donor_id_status === undefined) {
    row.donor_id_status = "verified_public_identifier";
  }
  rows.push(Object.fromEntries(columns.map((column) => [column, clean(row[column])])));
}

async function geo(accession) {
  return JSON.parse(
    await readFile(path.join(cacheDir, `${accession}_samples.json`), "utf8"),
  );
}

const gse201348 = await geo("GSE201348");
const supplement201348 = JSON.parse(
  await readFile(path.join(cacheDir, "GSE201348_supplement_tables.json"), "utf8"),
);
const s2 = supplement201348.tables.TableS2;
const s2Header = s2[1];
const s2Rows = s2.slice(2).filter((record) => record[0]);
const s2Map = new Map(
  s2Rows.map((record) => [
    record[0],
    Object.fromEntries(s2Header.map((header, index) => [header, record[index]])),
  ]),
);
const crcDonors = {
  "CRC-1-8810": "CRC1_8810",
  "CRC-2-15564": "CRC2_15564",
  "CRC-3-11773": "CRC3_11773",
  "CRC-4-8456": "CRC4_8456",
};
const crcSampleAliases = {
  CRC1_8810: "CRC-1-8810",
  CRC2_15564: "CRC-2-15564",
  CRC3_11773: "CRC-3-11773",
  CRC4_8456: "CRC-4-8456",
};

for (const sample of gse201348.samples) {
  const biological = sample.title.split(",")[0].trim();
  const supplementSample = crcSampleAliases[biological] ?? biological;
  const pathology = s2Map.get(supplementSample);
  if (!pathology) throw new Error(`GSE201348 supplement match missing: ${biological}`);
  const stage = sample.characteristics["disease stage"];
  const isReplicate = /Replicate/i.test(sample.title);
  const donor = crcDonors[supplementSample] ?? pathology.Donor;
  const histology =
    pathology.PolypType ??
    pathology.Pathologist0_PolypType ??
    (stage === "CRC" ? "Adenocarcinoma" : stage);
  add({
    accession: "GSE201348",
    sample_id: sample.sample_id,
    donor_id: donor,
    biological_sample_id: biological,
    sample_role: isReplicate ? "technical_replicate_library" : "biological_tissue_library",
    technical_replicate_of: isReplicate ? biological : NA,
    paired_group: donor,
    tissue: "colorectal mucosa",
    condition:
      stage === "CRC"
        ? "cancer"
        : stage === "Polyp"
          ? "polyp_or_adenoma"
          : stage === "Unaffected"
            ? "unaffected_mucosa"
            : "normal_mucosa",
    histology,
    colon_or_rectum: regionClass(pathology.Location),
    tumor_location: pathology.Location,
    sporadic_or_FAP:
      sample.characteristics["familial adenomatous_polyposis"] === "Y"
        ? "FAP"
        : stage === "CRC"
          ? "sporadic"
          : "non-FAP",
    treatment_status: NA,
    platform: "GPL24676; Illumina NovaSeq 6000",
    assay_type: "single-nucleus RNA-seq; 10x Genomics",
    data_format: "processed MTX/TSV counts; SRA raw reads",
    raw_counts_available: "yes",
    clinical_data_available: "yes_limited",
    validation_split: "discovery",
    download_url: sourceUrl(sample.sample_id),
    inclusion: isReplicate
      ? "include_for_technical_replicate_resolution"
      : "include",
    exclusion_reason: isReplicate
      ? "technical replicate; not an independent biological replicate"
      : NA,
    source_evidence:
      `${sourceUrl(sample.sample_id)}; Nature Genetics supplementary Tables S1-S3/S5`,
    metadata_notes:
      biological === "A001-C-007"
        ? "Gross pathology is adenocarcinoma; supplementary microscopic fields are internally discordant and require pathology-aware review."
        : NA,
  });
}

const gse161277 = await geo("GSE161277");
for (const sample of gse161277.samples) {
  const donor = sample.title.match(/^Patient\d+/)?.[0] ?? NA;
  const stage = sample.characteristics["tumor stage"];
  const blood = stage === "blood";
  add({
    accession: "GSE161277",
    sample_id: sample.sample_id,
    donor_id: donor,
    biological_sample_id: sample.sample_id,
    sample_role: "biological_sample",
    technical_replicate_of: NA,
    paired_group: donor,
    tissue: blood ? "blood" : "colorectal tissue",
    condition: stage,
    histology:
      stage === "carcinoma"
        ? "carcinoma"
        : stage === "adenoma"
          ? "adenoma"
          : stage,
    colon_or_rectum: NA,
    tumor_location: NA,
    sporadic_or_FAP: NA,
    treatment_status: NA,
    platform: "GPL20795; Illumina HiSeq X Ten",
    assay_type: "single-cell RNA-seq; 10x Genomics",
    data_format: "processed MTX/TSV counts",
    raw_counts_available: "yes",
    clinical_data_available: "yes_limited",
    validation_split: "single_cell_validation",
    download_url: sourceUrl(sample.sample_id),
    inclusion: blood ? "exclude_from_tissue_sequence" : "include",
    exclusion_reason: blood ? "blood sample is outside the tissue progression endpoint" : NA,
    source_evidence: sourceUrl(sample.sample_id),
    metadata_notes:
      donor === "Patient3" && /^Patient3 adenoma/.test(sample.title)
        ? "One of two distinct adenoma tissues from Patient3; not a technical replicate."
        : NA,
  });
}

const gse132465 = await geo("GSE132465");
for (const sample of gse132465.samples) {
  const c = sample.characteristics;
  const tumor = c["tissue type"] === "Colorectal cancer";
  add({
    accession: "GSE132465",
    sample_id: sample.sample_id,
    donor_id: c.patient_id,
    biological_sample_id: sample.sample_id,
    sample_role: "biological_tissue",
    technical_replicate_of: NA,
    paired_group: c.patient_id,
    tissue: tumor ? "colorectal cancer tissue" : "normal colorectal mucosa",
    condition: tumor ? "cancer" : "normal_mucosa",
    histology: tumor ? c.pathologic : "normal_mucosa",
    colon_or_rectum: tumor ? regionClass(c.region) : NA,
    tumor_location: tumor ? c.region : NA,
    sporadic_or_FAP: NA,
    treatment_status: NA,
    platform: "GPL20301; Illumina HiSeq 4000",
    assay_type: "single-cell RNA-seq; 10x Chromium 3' v2",
    data_format: "processed count matrix and cell annotation",
    raw_counts_available: "yes",
    clinical_data_available: "yes_limited",
    validation_split: "cell_specificity_validation",
    download_url: sourceUrl(sample.sample_id),
    inclusion: "include",
    exclusion_reason: NA,
    source_evidence: sourceUrl(sample.sample_id),
    metadata_notes:
      tumor
        ? NA
        : "Normal-mucosa record has no sample-level anatomic region in GEO; paired donor tumor location must not be copied to the normal sample.",
  });
}

const gse41657 = await geo("GSE41657");
for (const sample of gse41657.samples) {
  const grade = sample.characteristics["pathologic grade"];
  const condition = grade.includes("Normal")
    ? "normal_mucosa"
    : grade.includes("dysplasia")
      ? "adenoma"
      : "cancer";
  add({
    accession: "GSE41657",
    sample_id: sample.sample_id,
    donor_id: `UNRESOLVED_GSE41657_${sample.sample_id}`,
    donor_id_status: "unresolved_sample_scoped_placeholder",
    biological_sample_id: sample.sample_id,
    sample_role: "biological_tissue",
    technical_replicate_of: NA,
    paired_group: NA,
    tissue: "colorectal tissue",
    condition,
    histology: grade,
    colon_or_rectum: NA,
    tumor_location: NA,
    sporadic_or_FAP: "sporadic",
    treatment_status: NA,
    platform: "GPL6480; Agilent-014850 4x44K G4112F",
    assay_type: "expression microarray",
    data_format: "raw feature-extraction TXT; processed sample table",
    raw_counts_available: "not_applicable_microarray",
    clinical_data_available: "yes_limited",
    validation_split: "bulk_validation",
    download_url: sourceUrl(sample.sample_id),
    inclusion: "include",
    exclusion_reason: NA,
    source_evidence: sourceUrl(sample.sample_id),
    metadata_notes:
      "Patient identifier and colon-versus-rectum site are not reported in GEO. donor_id is a sample-scoped traceability placeholder, not evidence that each array came from a distinct patient.",
  });
}

const gse100179 = await geo("GSE100179");
for (const sample of gse100179.samples) {
  const tissue = sample.characteristics.tissue;
  const condition = /healthy/i.test(tissue)
    ? "normal_mucosa"
    : /adenoma/i.test(tissue)
      ? "adenoma"
      : "cancer";
  add({
    accession: "GSE100179",
    sample_id: sample.sample_id,
    donor_id: `UNRESOLVED_GSE100179_${sample.sample_id}`,
    donor_id_status: "unresolved_sample_scoped_placeholder",
    biological_sample_id: sample.sample_id,
    sample_role: "biological_biopsy",
    technical_replicate_of: NA,
    paired_group: NA,
    tissue,
    condition,
    histology: tissue,
    colon_or_rectum: /healthy colon/i.test(tissue) ? "colon" : NA,
    tumor_location: NA,
    sporadic_or_FAP: NA,
    treatment_status: NA,
    platform: "GPL17586; Affymetrix Human Transcriptome Array 2.0",
    assay_type: "expression microarray",
    data_format: "raw CEL; processed matrix",
    raw_counts_available: "not_applicable_microarray",
    clinical_data_available: "yes_limited",
    validation_split: "bulk_validation",
    download_url: sourceUrl(sample.sample_id),
    inclusion: "include",
    exclusion_reason: NA,
    source_evidence: sourceUrl(sample.sample_id),
    metadata_notes:
      "No donor identifier or matching key is reported in GEO; numeric title suffixes were not interpreted as patient matching. donor_id is a sample-scoped traceability placeholder and must not be used to claim independent patients.",
  });
}

const gse8671 = await geo("GSE8671");
for (const sample of gse8671.samples) {
  const c = sample.characteristics;
  const donor = `GSE8671_${c["Patient ID"]}`;
  const location = c.Location;
  add({
    accession: "GSE8671",
    sample_id: sample.sample_id,
    donor_id: donor,
    biological_sample_id: sample.sample_id,
    sample_role: "biological_tissue",
    technical_replicate_of: NA,
    paired_group: donor,
    tissue: c.Tissue,
    condition: c.Tissue === "normal" ? "normal_mucosa" : "adenoma",
    histology: c.Tissue === "normal" ? "normal_mucosa" : "sporadic_adenoma",
    colon_or_rectum: regionClass(location),
    tumor_location: location,
    sporadic_or_FAP: "sporadic",
    treatment_status: NA,
    platform: "GPL570; Affymetrix Human Genome U133 Plus 2.0",
    assay_type: "expression microarray",
    data_format: "raw CEL; processed matrix",
    raw_counts_available: "not_applicable_microarray",
    clinical_data_available: "yes_limited",
    validation_split: "bulk_validation",
    download_url: sourceUrl(sample.sample_id),
    inclusion: "include",
    exclusion_reason: NA,
    source_evidence: sourceUrl(sample.sample_id),
    metadata_notes: "Normal and adenoma tissues are paired by GEO Patient ID.",
  });
}

const gse99573 = await geo("GSE99573");
for (const sample of gse99573.samples) {
  const disease = sample.characteristics["disease status"];
  const split = sample.characteristics.set;
  const description = Array.isArray(sample.description)
    ? sample.description[0]
    : sample.description;
  const donor = description?.split(/\s+/)[0] ?? NA;
  const notUsed = split === "Not Used";
  add({
    accession: "GSE99573",
    sample_id: sample.sample_id,
    donor_id: donor,
    biological_sample_id: sample.sample_id,
    sample_role: "biological_stool_sample",
    technical_replicate_of: NA,
    paired_group: donor,
    tissue: "stool",
    condition:
      disease === "Normal"
        ? "normal"
        : disease === "Adenoma"
          ? "adenoma"
          : disease === "Benign"
            ? "benign_not_used"
            : "cancer",
    histology: disease,
    colon_or_rectum: NA,
    tumor_location: NA,
    sporadic_or_FAP: NA,
    treatment_status: NA,
    platform: "GPL17586; Affymetrix Human Transcriptome Array 2.0",
    assay_type: "stool eukaryotic RNA expression microarray",
    data_format: "raw CEL; processed matrix",
    raw_counts_available: "not_applicable_microarray",
    clinical_data_available: "yes",
    validation_split: split.toLowerCase().replaceAll(" ", "_"),
    download_url: sourceUrl(sample.sample_id),
    inclusion: notUsed
      ? "exclude"
      : split === "Testing"
        ? "reserve_locked_test_set"
        : "include_training_only",
    exclusion_reason: notUsed
      ? "GEO set field is Not Used; outside the published 330-patient training/testing cohort"
      : NA,
    source_evidence: sourceUrl(sample.sample_id),
    metadata_notes:
      split === "Testing"
        ? "Must remain untouched until the stool model is locked."
        : NA,
  });
}

const gse226997 = await geo("GSE226997");
for (const sample of gse226997.samples) {
  const donor = sample.title.match(/patient\s+\d+/i)?.[0].replace(/\s+/g, "_") ?? NA;
  add({
    accession: "GSE226997",
    sample_id: sample.sample_id,
    donor_id: donor,
    biological_sample_id: sample.sample_id,
    sample_role: "spatial_tissue_section",
    technical_replicate_of: NA,
    paired_group: donor,
    tissue: sample.characteristics.tissue,
    condition: "cancer",
    histology: "colorectal_cancer",
    colon_or_rectum: NA,
    tumor_location: NA,
    sporadic_or_FAP: NA,
    treatment_status: NA,
    platform: "GPL30173; Illumina NextSeq 2000",
    assay_type: "10x Genomics Visium spatial transcriptomics",
    data_format: "Space Ranger processed matrices/images; raw reads available",
    raw_counts_available: "yes",
    clinical_data_available: "yes_limited",
    validation_split: "optional_spatial_validation",
    download_url: sourceUrl(sample.sample_id),
    inclusion: "optional_include",
    exclusion_reason: NA,
    source_evidence: sourceUrl(sample.sample_id),
    metadata_notes: "GEO does not identify colon versus rectum or the precise tumor location.",
  });
}

const tcga = JSON.parse(
  await readFile(path.join(cacheDir, "TCGA-COAD_rnaseq_files.json"), "utf8"),
);
const tcgaSamples = new Map();
for (const file of tcga.data.hits) {
  for (const caseRecord of file.cases) {
    for (const sample of caseRecord.samples) {
      const key = sample.submitter_id;
      if (!tcgaSamples.has(key)) {
        tcgaSamples.set(key, {
          donor: caseRecord.submitter_id,
          sample,
          files: [],
        });
      }
      tcgaSamples.get(key).files.push({
        file_id: file.file_id,
        file_name: file.file_name,
      });
    }
  }
}
for (const [sampleId, record] of [...tcgaSamples].sort()) {
  const sampleType = record.sample.sample_type;
  const include = ["Primary Tumor", "Solid Tissue Normal"].includes(sampleType);
  add({
    accession: "TCGA-COAD",
    sample_id: sampleId,
    donor_id: record.donor,
    biological_sample_id: sampleId,
    sample_role: "bulk_RNA_tissue",
    technical_replicate_of: NA,
    paired_group: record.donor,
    tissue: sampleType,
    condition:
      sampleType === "Solid Tissue Normal"
        ? "normal_mucosa"
        : sampleType === "Primary Tumor"
          ? "primary_cancer"
          : sampleType.toLowerCase().replaceAll(" ", "_"),
    histology: sampleType === "Solid Tissue Normal" ? "normal" : "colon_adenocarcinoma",
    colon_or_rectum: "colon",
    tumor_location: NA,
    sporadic_or_FAP: NA,
    treatment_status: NA,
    platform: "Illumina RNA-seq; GDC harmonized",
    assay_type: "bulk RNA-seq",
    data_format: "STAR augmented gene counts TSV",
    raw_counts_available: "yes",
    clinical_data_available: "yes",
    validation_split: "auxiliary_clinical_validation",
    download_url: record.files
      .map((file) => `https://api.gdc.cancer.gov/data/${file.file_id}`)
      .join("|"),
    inclusion: include ? "include_auxiliary" : "exclude_from_primary_endpoint",
    exclusion_reason: include
      ? NA
      : `${sampleType} is outside the prespecified primary-tumor/normal auxiliary comparison`,
    source_evidence:
      "https://portal.gdc.cancer.gov/projects/TCGA-COAD; GDC files API metadata query",
    metadata_notes:
      record.files.length > 1
        ? `${record.files.length} STAR-count files map to this biological sample barcode; resolve file-level replicate/version choice before analysis.`
        : "Sample-level anatomic subsite and treatment fields require a later clinical-table join.",
  });
}

rows.sort((a, b) =>
  `${a.accession}\t${a.sample_id}`.localeCompare(`${b.accession}\t${b.sample_id}`),
);
const manifest = [
  columns.join("\t"),
  ...rows.map((row) => columns.map((column) => tsvCell(row[column])).join("\t")),
].join("\n");
await writeFile(path.join(outputDir, "dataset_manifest.tsv"), `${manifest}\n`, "utf8");

const exclusionColumns = [
  "accession",
  "sample_id",
  "donor_id",
  "inclusion",
  "exclusion_reason",
  "source_evidence",
];
const excluded = rows.filter(
  (row) =>
    row.inclusion === "exclude" ||
    row.inclusion === "exclude_from_tissue_sequence" ||
    row.inclusion === "exclude_from_primary_endpoint",
);
const exclusionLog = [
  exclusionColumns.join("\t"),
  ...excluded.map((row) =>
    exclusionColumns.map((column) => tsvCell(row[column])).join("\t"),
  ),
].join("\n");
await writeFile(path.join(outputDir, "exclusion_log.tsv"), `${exclusionLog}\n`, "utf8");

const counts = Object.fromEntries(
  [...new Set(rows.map((row) => row.accession))].map((accession) => [
    accession,
    rows.filter((row) => row.accession === accession).length,
  ]),
);
console.log(JSON.stringify({ rows: rows.length, excluded: excluded.length, counts }, null, 2));
