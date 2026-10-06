Dear @S-Lucas,

Thank you for reporting this issue in such detail! Your description, along with the debugging steps you tried, made it straightforward to trace the root cause.

### The Root Cause

The behavior you observed is due to a coercion and comparison bug in the underlying HyPhy Batch Language (HBL) file `CleanStopCodons.bf` (the script executed by the `hyphy cln` analysis shortcut).

Specifically, the script prompts the user for duplicate/gap filtering options using `io.SelectAnOption`:
```hyphy
filteringOption = io.SelectAnOption ({"No/No" :  "Keep all sequences and sites",
                                      "No/Yes" :   "Keep all sequences, filter sites with nothing but gaps",
                                      "Yes/No" :  "Filter duplicate sequences but keep all sites",
                                      "Yes/Yes" :  "Filter duplicate sequences and sites with nothing but gaps",
                                      "Disallow stops" : "Filter duplicate sequences and all sequences that have stop codons"},"Filter duplicates/gaps?");
```

1. **String/Numeric Mismatch**: `io.SelectAnOption` returns the selected key from the options dictionary as a **string** (e.g., `"No/No"`, `"No/Yes"`). However, the rest of the script assumes `filteringOption` is a **numeric index** from `0` to `4` and performs operations like `filteringOption >= 2` (to check if duplicate filtering is enabled) and `filteringOption % 2` (to check if gap filtering is enabled).
2. **Implicit Coercion Behavior**: In HBL, comparing a string to a number (e.g., `filteringOption >= 2`) is performed by converting the number to a string and executing an alphabetical comparison (i.e. comparing the choice string to `"2"`). Since `"No/No"`, `"No/Yes"`, `"Yes/No"`, `"Yes/Yes"`, and `"Disallow stops"` all start with characters (`N`, `Y`, `D`) that alphabetically succeed `"2"`, the condition `filteringOption >= 2` always evaluated to `1` (True) across all menu selections, causing duplicate sequence filtering to run unconditionally.
3. **Modulo and Typo Bugs**: In HBL, the `%` operator when applied to a string left-hand operand acts as a case-insensitive string equality check (`x % arg`). Since the right-hand operand `2` is a number (not a string), the comparison `filteringOption % 2` always returned `0` (False), meaning gap filtering was never executed. Furthermore, there was a typo on line 157 (`filterinOption == 4` instead of `filteringOption == 4`), which meant the `"Disallow stops"` option failed to correctly filter out sequences containing stop codons.

This explains why all five option settings produced the exact same output alignment where duplicate sequences were removed, and sequences with stop codons (which were converted to gaps and consequently became identical to other sequences) were filtered as duplicates.

### The Fix

We have resolved these bugs in `CleanStopCodons.bf`. The fix will be pushed in the next version and is currently available on the `develop` branch. We addressed the issue by comparing `filteringOption` directly to specific string choices and setting semantic boolean variables for downstream logic:
```hyphy
cln.filter_duplicates = (filteringOption == "Yes/No" || filteringOption == "Yes/Yes" || filteringOption == "Disallow stops");
cln.filter_gaps       = (filteringOption == "No/Yes" || filteringOption == "Yes/Yes");
cln.disallow_stops    = (filteringOption == "Disallow stops");
```

We then replaced the numeric comparisons and the typo with these boolean flags:
- `cln.filter_duplicates` is used to decide whether duplicate sequences should be filtered.
- `cln.filter_gaps` is used to decide whether gappy sites should be filtered.
- `cln.disallow_stops` is used to decide whether sequences with stop codons should be discarded.

These changes ensure that all five options behave exactly as documented in the menu (e.g. `"No/No"` now correctly preserves all sequences and sites).

Best regards,
