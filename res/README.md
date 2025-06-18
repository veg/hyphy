# Adding New Official Analyses to HyPhy

This README describes the process for adding new official analyses to HyPhy. Following these steps will ensure your analyses are properly registered and accessible via the command-line interface (`hyphy --help`) and the interactive menu.

## Steps to Add a New Analysis

1. **Add Analysis Files**
   - Place your analysis files (`.bf` only) in the `res/TemplateBatchFiles` directory.
   - You may organize analyses into subdirectories within `res/TemplateBatchFiles` and include Markdown-formatted READMEs in the subdirectories if desired.

2. **Add a Help URL**
   - In the analysis description at the top of your `.bf` analysis file, add a help URL to the analysis header. This provides users with documentation or further information about the analysis.
   - Example:
     ```
     terms.io.help : "https://your.documentation.url/analysis-name"
     ```

3. **Register the Analysis in `files.lst`**
   - Edit the `files.lst` file (located in the same directory) to include your new analysis file(s).
   - This step is required for HyPhy to recognize and list the analysis.

## Result

After completing these steps, your new analysis will:
- Appear when running `hyphy --help`.
- Be available in the interactive menu when launching HyPhy without arguments or through the menu system.

---

For further questions, refer to the main HyPhy documentation or contact the development team.