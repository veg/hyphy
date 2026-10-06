# Modernization Report: AnalyzeNucProtData.bf

## Baseline
- **File**: `res/TemplateBatchFiles/AnalyzeNucProtData.bf`
- **Original Behavior**: Executed a simple linear prompt sequence using `SetDialogPrompt` and `inputOptions`. Formatted standard output manually and did not spool structured JSON results. Standard outputs were just basic likelihood values.
- **Dependencies**: Inherited old legacy patterns that leaked global namespaces and polluted `stdinRedirect`.

## Revisions
- **(a) KeywordArgument**: Converted all interactive prompts to `KeywordArgument` in their exact chronological order (`alignment`, `model`, `tree`, `display-options`, `output`) to support fully non-interactive CLI workflows.
- **(b) Modern Conventions**: 
  - Integrated `libv3` loaders and standardized table formats.
  - Implemented a path resolution helper `ad.resolve_path` to handle relative inputs.
  - Placed execution blocks in an isolated run sequence (using `dummy_redirect = "isolate"`) to prevent parent standard input redirects from propagating incorrectly to child script calls via `ExecuteAFile`.
- **(c) MD Tables Output**: Utilized `io.FormatTableRow` to print markdown tables for:
  - Model Fit Summary (Log Likelihood, Parameters, AIC, AIC-c, BIC)
  - Global Parameter Estimates
  - Branch Lengths
- **(d) JSON Output**: Built a structured JSON dictionary mapping `libv3` terms (e.g. `terms.json.analysis`, `terms.json.input`, `terms.json.fits`, `terms.json.trees`) and spooled it using `io.SpoolJSON` to the requested path.
- **(e) Analysis Banner**: Embedded an informative analysis banner containing version, authors, and requirements.

## Documentation
- Examples and walkthroughs were generated and appended to the main documentation tracking.
- Test suites (`TwoSequenceTest.bf`) were adjusted to pass given the new table and prompt output format.
