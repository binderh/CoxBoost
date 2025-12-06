## R CMD check results

0 errors | 0 warnings | 0 note

* John Zobolas (mlr3 core team) takes over maintenance after discussion with original owner (Harald Binder)
* I updated parts of the code to more recent standards - no change in functionality
* Motivation: this is a great model for single-event survival data and competing risks, would be great to have it back to CRAN!

## CRAN Review comments

Thanks for the comments! Here is a summary of the fixes:

* Added two publications in the description field of the DESCRIPTION file. CITATION file added as well.
* Added the two suggested missing @return fields.
* As most examples are very "heavy" (i.e. execution time more than 1min), \dontrun{} was replaced with \donttest{}.
* The `trace` argument was added in resample.CoxBoost(), following the previous convention in the CoxBoost() and cv.CoxBoost() functions.
* The set.seed() calls in resample.CoxBoost() were removed.
