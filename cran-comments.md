# Resubmission

Update for current submission: We have substantially removed the size of test
files for the examples, reducing overall package size from ~30 Mb to only ~6 Mb.
A further reduction was not possible as the package contains examples for 
many different NGS file types.

This package has been submitted previously (2023, v0.7.1) and was removed from
CRAN due to several issues. In the mean time many functions were re-factored,
more than 10 dependencies were removed to make the package lighter, and other
problems regarding documentation and style were fixed. Current version 1.4.0
now builds fine on the tested platforms.

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
