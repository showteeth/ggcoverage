# ggcoverage 1.1.0
## Major changes
* Mark SNV with twill or star.
* Moved `FormatTrack` to `LoadTrackFile`, this can reduce load time and memory for big files.

-------------

# ggcoverage 1.0.0
## Major changes
* Supporting Hi-C visualization (`geom_tad` and `geom_tad2`).
* Supporting genome region interaction visualization (`geom_link`).
* Supporting CNV visualization (`geom_cnv`).

## Minor changes
* `geom_ideogram` supporting the highlight of centromere.

-------------

# ggcoverage 0.9.0
## Major changes
* Supporting gtf file from Ensembl

## Minor changes
* Fixed bugs in getting gene group.

-------------

# ggcoverage 0.8.1
## Minor changes
* Fixed bugs in `ggcoverage`.

-------------

# ggcoverage 0.8.0
## Major changes
* Changed the plot type of coverage when visualizing at single-nucleotide level (`geom_bar` instead of `geom_rect`)
* Changed the plot type of amino acids (`geom_rect` instead of `geom_tile`)
* Changed x axis range
* Added `GetPlotData` to obtain raw plot data
* Added a new vignette to customize the plot

## Minor changes
* Fixed bugs in `geom_base`.

-------------

# ggcoverage 0.7.2
## New features
* Added `GetConsensusPeak` to get consensus peak from replicates with MSPC.

## Minor changes
* Added `peak.df` for `geom_peak` to support dataframe as input.

-------------

# ggcoverage 0.7.0
## New features
* Added `geom_base` to annotate genome coverage with base and amino acids.

## Minor changes
* `LoadTrackFile` supported visualization at single-nucleotide level.
* Added `rect.color` for `geom_coverage` to control  rect border color.
* Fixed bug in track file format identification.

-------------

# ggcoverage 0.6.0
## New features
* Added `geom_gc` to annotate genome coverage with GC content. 

-------------

# ggcoverage 0.5.0
## Minor changes
* Fixed bug in `getIdeogram`.

-------------

# ggcoverage 0.4.0
## New features
* Added `geom_peak` to enhance its usage in ChIP-seq or ATAC-seq data.
* Changed Y axis theme.

-------------

# ggcoverage 0.3.0
## New features
* Added `geom_transcript`, `geom_ideogram`.

## Minor changes
* Fixed bug in `GetGeneGroup`.

-------------

# ggcoverage 0.2.0
## New features
* Added `geom_gene`, `geom_ideogram`.

-------------

# ggcoverage 0.1.0

## New features
* Added a `NEWS.md` file to track changes to the package.
* Added `ggcoverage`, `geom_coverage`, `LoadTrack` and `FormatInput`.

-------------
