# Modernization Report: Selection Analyses (FitModel, MEME, FEL)

## Baseline
- **Files**: 
  - `res/TemplateBatchFiles/SelectionAnalyses/FitModel.bf`
  - `res/TemplateBatchFiles/SelectionAnalyses/MEME.bf`
  - `res/TemplateBatchFiles/SelectionAnalyses/FEL.bf`
- **Original Behavior**: `FitModel.bf` was originally a standalone script in a separate repo that needed porting. `MEME.bf` and `FEL.bf` were mostly using `libv3` but contained bugs regarding stdout formatting, log-likelihood swapping, and regex matching for `--limit-to-sites`.

## Revisions
- **(a) KeywordArgument**: All scripts fully support non-interactive arguments. `FitModel.bf` dynamically handles Likelihood Ratio Test overrides by generating `KeywordArgument` prompts for global parameter bounds checking dynamically.
- **(b) Modern Conventions**: 
  - `FitModel.bf` ported fully into standard library layout. Fixed module directory paths to resolve against `"modules/"` correctly. Added `files.lst` registry mapping.
  - `MEME.bf` & `FEL.bf`: Updated `selection.io.sitelist_matches_pattern` regex (`"(^|\\,)" + idx + "(\\,|$)"`) inside `io_functions.ibf` so filtering matches indices exactly instead of partially prefix matching.
- **(c) MD Tables Output**: 
  - `FitModel.bf` prints log likelihoods and model parameters as MD tables.
  - `MEME.bf`: Added missing `"LRT MEME vs FEL"` table header column and swapped the previously incorrect column alignments. Synchronized optimization output weights to the global namespace to correct stale variables formatting.
- **(d) JSON Output**: All files accurately spool their components using standard `libv3` `io.SpoolJSON` functionality.

## Documentation
- The `FitModel.bf` script updates the `terms.io.info` banner to document its dynamic LRT feature (`--parameter_name value`).
- Testing routines mapped to successfully execute using standard `tests/hbltests/libv3/runner.bf`.
