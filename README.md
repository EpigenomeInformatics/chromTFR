# chromTFR

Quantification of transcription factor footprints in chromatin
accessibility data, using the annotation of
[methylTFR](https://github.com/EpigenomeInformatics/methylTFR).

Tn5 insertion densities replace methylation levels. For every binding
site of a motif the insertions are counted per base and stacked on the
site centre, and the deviation score is the insertion density of the
central window over the density of the outer flanks, corrected by the
same ratio on an expected profile. The expected profile comes either
from the GC bin composition of the binding sites, as in methylTFR, or
from a k-mer model of the Tn5 sequence preference.

## Installation

```r
remotes::install_github("EpigenomeInformatics/chromTFR")
```

Reading a dataset requires [ChrAccR](https://github.com/GreenleafLab/ChrAccR),
and the k-mer model requires a `BSgenome` package. Both are suggested
rather than required, so the package installs without them.

## Usage

```r
library(chromTFR)

dsa <- loadAccDataset("dsATAC_filtered")
peaks <- getAccRegions(dsa, regionType = ".peaks.cons", extend = 500)
ins <- getTn5Insertions(dsa, getAccSamples(dsa)[1], regions = peaks)

accDeviationScore(accProfile(ins, prepareTFBS(tf_bindsites[["CTCF"]])))
```

## Citation

```r
citation("chromTFR")
```

The vignette *A worked example with the ChrAccR example data* runs the
whole workflow on the public `ChrAccRex` dataset, from the `DsATAC`
object to the footprints.

See the vignettes for the full workflow, and
`system.file("scripts", package = "chromTFR")` for template scripts
covering the steps from a dataset to the figures.
