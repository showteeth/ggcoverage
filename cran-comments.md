# Resubmission

This package has been submitted previously (2023, v0.7.1) and was removed
from CRAN due to several issues, which have been addressed now.

## Issues

### Update for current submission:

- fixed DESCRIPTION abstract, added quotation marks for package names etc
- examples that are running >5 sec were wrapped in `/donttest{}` instead of `/dontrun{}`
- example in `ggcoverage.Rd` does not contain out-commented code as stated by the reviewer
- changing of graphical parameters using `par(mfrow=c(2,2))` were removed from vignette
- however the package does not contain any `par()` statement in any of the
  functions (supposedly `R/geom_base.R`) as stated by the reviewer

### Previously fixed issues

- We have substantially removed the size of test files for the examples,
reducing overall package size from ~30 Mb to only ~6 Mb. A further reduction
was not possible as the package contains examples for many different NGS file
types.

- Examples were reduced to run in shorter time or excluded from running at all.
Please note that the 5 sec threshold feels very arbitrary and heavily depends
on the test environment -- this package never exceeds the threshold when
tested locally on a 3 year old average laptop, or in 3 different github actions 
workflows, but regularly fails on the CRAN server.

## Test environments

### Local

- Ubuntu 22.04

### with Github Actions

- windows-latest (release)
- macOS-latest (release)
- ubuntu-latest (release)

## R CMD check results

0 errors | 0 warnings | 1 note

There was 1 NOTE:

```
installed size is  5.8Mb
  sub-directories of 1Mb or more:
    doc       1.6Mb
    extdata   3.0Mb
```

## Downstream dependencies

- There are currently no downstream dependencies for this package
