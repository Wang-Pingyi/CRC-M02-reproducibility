# Stage 10FG-CLOSE decision

Decision: **PASS**

## Authoritative branch closure

- Stage 10F is locked as `SKIPPED_NOT_ESTIMABLE`.
- Stage 10G is locked as `SKIPPED_GATE_NOT_MET`.
- Stage 10E-DESC is locked as `DESCRIPTIVE_ONLY`.
- Stage 10H is `AUTHORIZED` but was not started by this stage.

## Stage 10F rationale

Only one of four patients passed the frozen ROI registration gate. This is
below the prespecified minimum of three independent patients, so the
patient-level primary spatial estimand is not estimable. This is a spatial ROI
registration sample-size limitation, not evidence that M02 is negative.

No Stage 10F numeric result table was created. No effect estimate, confidence
interval, P value or FDR exists for Stage 10F.

## Stage 10G rationale

There is no eligible Stage 10F patient-level spatial evidence. The case4
result is a single-patient descriptive localization and cannot open mechanism
analysis. Stage 10G is therefore skipped, and no Stage 10G numeric or
mechanistic result table was created.

## Claim boundary

Case4 must not be called a spatial validation cohort. It cannot support
spatial validation, population inference, mechanism, generalization or
biomarker claims. Stage 10H must preserve Stage 10F as not estimable and Stage
10G as gate not met; authorization of Stage 10H does not repair or replace
either missing branch.

## Historical-file integrity

The pre-existing skip markers
`results/stage10f/STAGE10F_SKIPPED.md` and
`results/stage10g/STAGE10G_SKIPPED.md` were already frozen as read-only inputs
to Stage 10E-DESC. They remain byte-for-byte unchanged. This decision and the
branch-closure table are the authoritative formal closure record.
