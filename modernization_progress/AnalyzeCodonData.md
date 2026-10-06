# Modernization Report: AnalyzeCodonData.bf

## Baseline
- **File**: `res/TemplateBatchFiles/AnalyzeCodonData.bf`
- **Original Behavior**: An extensive set of prompts mapping input configurations sequentially to load standard codon templates, generate trees, and fit sequences. Relied heavily on global standard inputs and older JSON formatting structures.
- **Dependencies**: Inherited multiple codon `.mdl` and `.ibf` templates (e.g. `MEC`, `GY94`, `MG94`, etc.).

## Revisions
- **(a) KeywordArgument**: Added standard options (`code`, `alignment`, `use-tree`, `tree`, `model`, `model-options`, `branch-lengths`, `output`). To accommodate `libv3` dynamic evaluation, arguments were rearranged precisely matching the chronological prompt order during execution.
- **(b) Modern Conventions**: 
  - Restructured to use `namespace acd { ... }` isolation.
  - Allowed skipping branch selections by defaulting `selected_branches = {}`.
  - Updated 23 supporting `.mdl` and `.ibf` codon templates inside `TemplateModels` to process `KeywordArgument` parsing without interactive breakages.
- **(c) MD Tables Output**: Outputs standard `libv3` markdown tables for:
  - Model Fit Summary (Calculates and displays `AIC-c` as well using codon size constraints)
  - Global Parameter Estimates
  - Branch Lengths
- **(d) JSON Output**: Converted old hardcoded JSON literal keys to map via the `terms.*` namespace. Saves to `.ACD.json` by default.
- **(e) Analysis Banner**: Displays the standard banner natively.

## Documentation
- Thorough documentation outlining the `AnalyzeCodonData.bf` refactor and usage of Neighbor-Joining tree fallback commands was added to `help/lib/README.md`.
