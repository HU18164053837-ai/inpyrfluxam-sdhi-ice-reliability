# Data dictionary and evidence units

Concentrations with `ug_L` or model `CI95_lower/upper` and `PI95_lower/upper` are in μg/L; `mg_L` is mg/L. Blank values mean unavailable or not evaluable, never zero. Boolean fields use True/False. ASCII field names are software identifiers. Compound names and scientific species names identify records; they do not imply independent experiments.

## Units and validation tiers

There are 8 compounds, 30 species–compound endpoint summaries and 90 dependent surrogate–target directions. The primary external-candidate tables contain 52 directions and 12 endpoints (34/7 inpyrfluxam, 18/5 isopyrazam). The separate conditional internal files contain 38 directions and 18 endpoints across 6 training-overlap compounds. A missing public-workbook match does not establish complete trial-level independence. Repeated regulatory citations or opposite directions do not add experiments.

## Saved 90-direction table

`all90_predictions_recomputed.csv` is the previously saved result table; its name does not mean it was recomputed during manuscript preparation. `compound`, `surrogate`, and `target` are direction keys. `validation` records the original tier. `input_ug_L`, `observed_ug_L`, and `predicted_ug_L` are surrogate observation, target observation and saved model point estimate. The model uses log₁₀ concentrations in μg/L.

`intercept`, `slope`, `R2`, `MSE`, `n_training`, `x_min`, `x_max`, `x_mean`, `Sxx` are saved regression parameters; the x-range and mean are log₁₀ values. Training n is not the number of validation trials. `leverage` is 1/n plus squared input departure from the training mean divided by Sxx. `extrapolation_distance_log10` records distance outside the numerical range; `in_domain` is numerical membership, not mechanistic applicability. `guidance_pass` is the saved screening flag.

`signed_log10_error` is log₁₀(predicted/observed), positive for underestimated continuous toxicity. `absolute_log10_error` is its absolute value. `fold_error` is max(predicted/observed, observed/predicted); `LC50_overestimated` is predicted > observed. `CI95_*` is the regression mean-response interval; `PI95_*` is the supplementary classical individual-response interval. `CI_fold_width` is upper/lower. `CI_covers` and `PI_covers` test whether the observed endpoint lies inside those different model intervals. Neither interval propagates bioassay uncertainty; the PI is not an independently calibrated EPA protective interval.

## Classification directions

`GB_hazard_classification_results.csv` has 52 external-candidate directions; the `GB_internal_` counterpart has 38 internal directions. Observed and predicted LC50 fields are unchanged saved concentrations. Category ranks run from 1 (very highly toxic) to 5 (practically nontoxic). In μg/L: rank 1 <100; rank 2 100–1000 inclusive; rank 3 >1000–10000; rank 4 >10000–100000; rank 5 >100000. `signed_class_shift` is predicted rank minus observed rank; positive means a less severe predicted description. `absolute_class_shift` counts category differences, and `point_class_mismatch` indicates unequal ranks. Category hazard under/overestimation differs from any continuous signed error.

`observed_` and `predicted_` distance fields identify the nearest threshold, absolute log₁₀ distance, fold distance, within-twofold flag and exact-boundary flag. Twofold is only a diagnostic convention. `target_CI_incompatible` is absent when no experimental CI is available; otherwise it asks whether the fixed predicted category lies outside the target interval's possible categories. `hazard_underestimated_for_all_target_CI_classes` asks whether the predicted rank exceeds every category permitted by that interval. These are not probabilities or complete uncertainty propagation.

`uncertainty_status` distinguishes point mismatch with target-CI incompatibility, compatibility or missing CI. Point matches have separate labels and do not prove agreement throughout an experimental interval.

## Experimental endpoints

`GB_endpoint_uncertainty_assessment.csv` contains 12 external endpoint summaries; internal endpoints are separate. `CI95_lower_ug_L/upper_ug_L` are actual reported target bioassay limits. `CI_status` is CI_SINGLE_CLASS, CI_MULTIPLE_CLASSES or CI_MISSING. Semicolon-separated possible classes form a set, not class probabilities. Touching a cutoff and spanning classes are distinguished by the fixed boundary rule. Missing inpyrfluxam CIs were not fabricated.

`specific_quality_limitation`, flags and notes document concentration deviations or replicate variability; they do not automatically invalidate a study. Study/report IDs, source URLs, concentration basis, exposure regime, material, purity, GLP and regulator reliability preserve reported information only. Missing age and water conditions remain unknown. Alternative regulatory point readings represent the same underlying trial and are not extra experiments or CIs.

## Paired baseline

`GB_baseline_comparison.csv` compares methods on identical external-candidate direction keys. `identity_prediction_ug_L` equals `surrogate_observed_LC50_ug_L`. The `webice_` and `identity_` fields report the same error, classification and CI-compatibility definitions under each method. `delta_*_webice_minus_identity` subtracts baseline from model: a negative error difference favors Web-ICE. A zero category difference can coexist with a continuous improvement. The baseline is not a validated regulatory model.

## Summaries and selection

Target summaries have one compound–scenario–target row, with within-target surrogate counts and mean absolute errors. Surrogate summaries group by the input species. Compound-direction summaries weight available directions; compound-equal-target summaries average target errors/rates equally among covered targets. `n_directions` is not independent sample size. Rates, ranges and counts of improved/tied targets describe different estimands. Zero CI-incompatible counts with missing CIs do not establish zero true disagreement; read the CI-missing counts alongside them.

Scenarios are ALL_AVAILABLE, GUIDANCE_AND_DOMAIN, STRICT_QUALITY (also MSE≤0.5 and n≥10), TARGET_AWAY_FACTOR2, and inpyrfluxam-specific bluegill omission checks. Both methods use the same subset. Strategy A minimizes MSE, B minimizes leverage, C integrates eligible log predictions, and D uses available trout or zebrafish surrogates. PRIMARY denotes the saved selection comparison; sensitivity scenario names retain the prior definitions. Selection did not receive target errors, but earlier outcomes were known, so the strategy comparison is exploratory.

## Structural and source tables

Pocket-alignment tables are sequence mappings with orthology limitations, not experimental fish structures. Redocking RMSD and contact recovery describe reference-ligand recovery. Proxy scores are saved human-scaffold controls, not binding free energies or four complete fish receptors. Fish model QC records rejected models. No new structural calculations were run.

Source facts remain attributed to their regulatory reports and verified articles. Full third-party reports, articles and training workbooks are not redistributed. Source-access history can differ from current evidence; the current Xu full-text comparison supersedes earlier unavailable-text notes. The data archive is local with no assigned public repository DOI.
