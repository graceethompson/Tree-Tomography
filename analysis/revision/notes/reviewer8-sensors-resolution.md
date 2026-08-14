# Reviewer 1, comment #8 — sensor counts, paths, and SoT resolution

Working notes for the response letter. Analysis: `scripts/sensor_counts.R`;
data: `data/sot_image_annotations.csv`; outputs: `output/CJFR-sensor-counts.{png,csv}`.

## What was actually deployed

Sensor counts were not recorded in the field data, but every PiCUS tomogram
prints the numbered sensor positions, so they were recovered from the images
and verified three ways:

1. **Visual reads** of all 61 main-survey scans and 24 validation scans
   (numbered labels, ambiguous digits re-read at magnification).
2. **Pixel detection**: the olive-green sensor-ring markers have a color
   signature (B < 45) distinct from the legend green and the decay-class green
   (B ≈ 86); clustering those pixels reproduces the visual count on 59/61
   scans. The two disagreements (351, 583) are pixel *over*-counts from
   decay-zone boundary pixels on heavily damaged trees — confirmed 7 by eye.
   No scan has a pixel count below the visual read (the direction a missed
   sensor would show).
3. **Plausibility**: 8 sensors occur only on the seven largest trees
   (EMS *Q. rubra*, DBH 35–43 cm); all others got 7. This matches the PiCUS
   circumference-based deployment recommendation, which was followed
   throughout (software defaults).

Result: **7 sensors (21 sensor-pair paths) for 50 of 57 trees; 8 sensors
(28 paths) for the seven largest EMS Q. rubra.** Validation hemlocks: 7–8
per scan. The reviewer's arithmetic assumed a minimum of 8 sensors / 28
paths; the true numbers are lower, so the resolution concern was, if
anything, understated — state this plainly rather than minimizing it.

## Key derived quantities (all 57 trees, in CJFR-sensor-counts.csv)

- Sensor spacing (circumference / count): mean 10.6 cm BGS vs 12.4 cm EMS
  (Welch p = 0.013); range 7.6–16.8 cm.
- DBH: BGS 23.7 ± 4.7 cm vs EMS 29.0 ± 8.3 cm (p = 0.006).
- Spacing vs percent damaged: Spearman ρ = 0.28, p = 0.038.
- percent_damaged ~ spacing + site (Type I): once spacing enters first, the
  site term attenuates to F = 3.14, p = 0.082.

## Position for the response

Concede, with quantification, rather than defend precision the method lacks:

- **Deployment was standard practice**: sensor count followed the
  manufacturer's circumference-based recommendation, as in the applied SoT
  literature; absolute spacing was therefore roughly constant (~11 cm)
  across the sample.
- **The classification does not hinge on sub-resolution detail**: SoT damage
  is effectively bimodal (0% or ≥6%, empty gap between — `boundaries.R`), so
  every detected defect is well above any plausible reconstruction floor,
  and the structural (>1%) classification is identical for any threshold in
  1–5% except 2 trees.
- **The confound is real and is now reported**: spacing covaries with site
  and weakly with the damage estimate; the site percent-damage contrast
  attenuates when spacing is controlled. This is consistent with (and
  further motivates) reframing the site contrast as exploratory rather than
  a headline finding.
- **Precision honesty**: report site-mean percent damage rounded to whole
  percent (not 1.93/7.07), and present percent-damage magnitude comparisons
  (species severity, site means) as exploratory given 21–28 paths per
  section and the resolution-severity confound.

## Related: PiCUS's own area tabulation (repro #2 tie-in)

Damaged trees' tomogram headers print both quantities separately, e.g. tree
351: "Solid wood: 48 % · Damaged: 26 %" (with 26% green implied as the
remainder), and tree 583: "Solid wood: 68 % · Damaged: 19 %". This is
PiCUS's own confirmation that its "Damaged" class excludes the green
intermediate zone — matching `percent_damaged` in Tree_ID_info.csv exactly
and settling how green was tabulated (see `green_zone.R`).
