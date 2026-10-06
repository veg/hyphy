# Modernization Report: AnalyzeDiNucData.bf

## Baseline
- **File**: `res/TemplateBatchFiles/AnalyzeDiNucData.bf`
- **Original Behavior**: Executed a simple interactive prompt for reading di-nucleotide alignments and fitting models. The model parameterization used `CreateFilter (ds, 2)`.

## Revisions
- **(a) KeywordArgument**: Adopted chronological keyword arguments matching `AnalyzeNucProtData.bf` (`alignment`, `model`, `tree`, `display-options`, `output`).
- **(b) Modern Conventions**: 
  - Restructured the file to match `libv3` standards using a specific variable prefix (`adn.`).
  - Added safe path resolution handling via `adn.resolve_path`.
- **(c) MD Tables Output**: Outputs model fit summary (AIC, AICc, LogL, BIC), Global Parameters, and Branch Lengths in standard Markdown tables.
- **(d) JSON Output**: Constructs and saves standard `libv3` compliant JSON data mapping to `.ADN.json` or custom output files.
- **(e) Analysis Banner**: Added `io.DisplayAnalysisBanner` declaring support for Di-nucleotide alignment data fitting.

## Documentation
- The modernized workflow was documented and validated using a MUSE95 model on a generated nucleotide string alignment grouped in pairs.
