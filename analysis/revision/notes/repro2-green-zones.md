# Green (intermediate-velocity) zones — what was recorded, and why it was conservative

Working notes for the response letter (reproducibility comment #2; also
supports #9 and #4). Analysis: `scripts/green_zone.R`; verification data:
`data/sot_image_annotations.csv`; figure: `output/CJFR-green-zone.png`.

## What was actually recorded (settled, three ways)

- `percent_solid_wood` in Tree_ID_info.csv equals the "Solid wood: X %"
  printed by PiCUS on the tomogram header (56/57 exact; one 1-point read
  discrepancy on tree 434). Solid wood = brown classes only.
- `percent_damaged` equals PiCUS's "Damaged" tabulation = violet + blue
  classes only. Damaged trees' headers print it directly (tree 351:
  "Solid wood: 48 % · Damaged: 26 %"; tree 583: "68 % · 19 %"), and
  solid + damaged < 100 for most damaged trees — the remainder is green
  (median 3%, up to 26% of the section).
- Independent pixel re-measurement of all 57 tomograms reproduces the
  recorded values only under the green-excluded definition
  (violet+blue vs percent_damaged: r = 0.993; green vs the remainder:
  r = 0.98).

**Manuscript correction required**: section 2.3.1 describes percent damage
as "non-brown colors combined." That wording is wrong — the recorded
quantity excluded green. State the actual definition: percent damage =
PiCUS decayed (violet) + cavity (blue) classes; intermediate-velocity green
counted as neither solid nor damaged.

## Why exclusion was the conservative choice (and helps #9)

Green is the intermediate-velocity class — the one most susceptible to
reconstruction artifact at 21–28 sensor-pair paths. Counting only the
unambiguous decay classes means every recorded damage call is well above
the plausible artifact floor:

- As recorded: smallest nonzero damage = 3%; **zero** trees in (0,1%];
  only two trees below 6%. The damage distribution is effectively bimodal
  (0% or ≥6% for 13 of 15 damaged trees) with an empty gap where artifacts
  would accumulate, which is why the >1% structural threshold is stable
  (identical classifications for any cut in 1–5% except 2 trees).
- Had green been counted as damage: 7 trees would fall in (0,1%] and 15
  more in (1,5%] — a near-threshold continuum, and 21 of 57 structural
  calls would flip. The crispness the classification relies on exists
  BECAUSE green was excluded.

Framing for the letter: "Percent damage counted only the unambiguous
decay/cavity velocity classes; intermediate-velocity (green) regions were
excluded to be conservative, since intermediate velocities are the most
susceptible to reconstruction artifact at these sensor counts. As a result
every recorded damage value is ≥3% of the section, and the structural
classification is insensitive to the threshold (1–5%)."

## Sensitivity (reviewer asked for it explicitly)

- Site percent-damage contrast is not sensitive to the green treatment:
  one-way site ANOVA p = 0.027 (recorded) vs p = 0.006 (green as damage);
  prevalence direction EMS > BGS under both.
- Green-as-damage flips 21/57 structural calls — reported in
  `green_zone.R` with the full list; all flips are 0%-damage trees with
  green shares of 2–26%.

## Tie-in to #4 (wetwood)

Green zones in otherwise sound wetland trees plausibly reflect elevated
moisture (wetwood) rather than degradation — the same ambiguity the
manuscript already concedes for the ERT incipient class. Excluding green
from the SoT damage measure keeps the structural axis free of that
ambiguity: the moisture story lives entirely on the ERT axis, where the
revision now frames it as incipient-or-wetwood.
