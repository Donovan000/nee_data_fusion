Step 1 - Extract and prepare data:
a) Run <extract_fluxnet.m> to read directly from the FluxNet daily csv files.
b) Run <collocate_data_types.m> to pull in RSIF data.
c) Run <biweekly_averaging.m> to convert from daily to twice-monthly values, so that all rows have an RSIF value.

Step 2 - Perform linear sensitivity analysis:
<main_lr_fluxnet.m> assesses sensitivity to different input columns. This is probably not useful for nonlinear regression.

Step 3 - Train LOO global regressions:
Train and test nonlinear regression with a leave-one-out evaluation strategy: there are 209 sites, and we train 209 different ANN regressions, each with 208 sites worth of training data. Each regression model is tested on the remaining site not used for training. 
Statistics are caluculated on a site-by-site basis, as well as on the whole LOO data record - i.e., all predicticted values from all sites.

