To install this package:

```
library(devtools)
install_github('jpwrobinson/funk')
library(funk)
```

To update this package run 'Check' in RStudio.

**Documentation is shit**.

Functions included are:

**Data cleaning**

* ```uniques``` generic function for counting unique objects (e.g. how many species?)
* ```scaler``` generic function for scaling and centering continuous covariates, and creating dummy variables for categorical covariates

**Data visualization**

* ```add_labels``` for annotating multipanel base plots (thanks [Sean Anderson](https://seananderson.ca/2013/10/21/panel-letters/))
* ```theme_sleek``` for a publication-quality ggplot theme (thanks [Sean Anderson](https://github.com/seananderson/ggsidekick/tree/master/R) and Cameron Freshwater)
* ```pairs2``` for a correlation matrix, showing Pearson correlations, histograms and scatterplots with LOESS smoothers

**Statistical modelling**

* ```dfa_fits``` examine model fit for a DFA object (built from [NWFSC](https://nwfsc-timeseries.github.io/AFTSLabbook/sec-dfa-plot-data.html) tutorial)
* ```fit_dfa``` fits and plots a DFA model (also from [NWFSC](https://nwfsc-timeseries.github.io/AFTSLabbook/sec-dfa-plot-data.html))
* ```gam_jacknife``` runs a jacknife sensitivity analysis for GAMs
* ```gam_predict_newdata``` predicts fitted GAM relationships for new data
* ```gam_predict``` performs model diagnostics and predicts fitted relationships for GAM objects
* ```gam_dev_exp``` evaluates relative variable importance for GAMs, following [Simon Wood](http://r.789695.n4.nabble.com/variance-explained-by-each-term-in-a-GAM-td836513.html) and [Litzow et al. 2014, *Glob. Change Biol*.](https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.12373)
* ```glm_gam_test``` compares model fits of GAM and GLM (by AIC)
* ```mmi_tvalue``` evaluates variable importance by multimodel inference, following [Cade (2015)](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/14-1639.1)
* ```vif_mer``` estimates variance inflation factors for lme4 models (thanks [Austin Frank](https://github.com/aufrank/R-hacks/blob/master/mer-utils.R))

