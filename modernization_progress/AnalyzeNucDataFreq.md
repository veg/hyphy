# Modernization Report: AnalyzeNucDataFreq.bf

## Baseline
- **File**: `res/TemplateBatchFiles/AnalyzeNucDataFreq.bf`
- **Original Behavior**: Used interactive menus to select standard models while uniquely exposing the equilibrium base frequencies (`estPiA`, `estPiC`, `estPiG`, `estPiT`) as global parameters. This allows nucleotide frequencies to be jointly optimized alongside standard substitution rates.
- **Dependencies**: Included a minor typo in the prompt string (prompting for "nuleotide data or aminoacid file" instead of strictly nucleotide).

## Revisions
- **(a) KeywordArgument**: Added chronological `KeywordArgument` parsing (`alignment`, `model`, `tree`, `display-options`, `output`). Corrected the typo in the prompt to match its actual logic. Re-aligned the model selection prompt to exact string "Select a standard model" to ensure correct keyword-to-input mapping.
- **(b) Modern Conventions**: 
  - Restructured to `libv3` namespace mapping using `afd.` prefixes to prevent collisions.
  - Used `EMBED_FREQUENCY_DEPENDENCE = 1` globally while preserving the estimation of base frequencies as standard globals within the analysis setup.
- **(c) MD Tables Output**: Added multiple `io.FormatTableRow` blocks:
  - Model Fit Summary (LogL, AIC, etc.)
  - Global Parameter Estimates
  - Branch Lengths
  - **Estimated Base Frequencies**: A new table cleanly displaying the Estimated vs Observed values for A, C, G, T.
- **(d) JSON Output**: Populates a `.AFD.json` output mapping the model fit results and custom frequencies correctly into the `terms.json.fits` nested structure.
- **(e) Analysis Banner**: Formatted standard banner.

## Documentation
- Documented testing outputs ensuring that the Estimated vs Observed Frequency markdown table outputs cleanly and accurately against the HKY85 standard models.
