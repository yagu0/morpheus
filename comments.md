There were 3 NOTEs:

1. "New submission / Package was archived on CRAN"
   - This version addresses all previous issues that led to archiving, at least as far as I can check (R CMD check --as-cran is OK). I forgot the previous issues.

2. "checking for future file timestamps ... NOTE: unable to verify current time"
   - This appears to be an environmental issue on my machine during check and is not a property of the package files themselves (I reset all timestamps with "find . -exec touch {} +" )

3. "checking compilation flags used ... NOTE: Compilation used the following non-portable flag(s): ..."
   - These flags are the default configuration of the R installation on my machine. They are not hardcoded in the package's src/Makevars.

====

The version only changes to 1.0.5 because I just added the pvalue() function to the file R/utils.R.
