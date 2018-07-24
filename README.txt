This file needs to be updated after recent significant changes (7-10-18). File should provide instructions for users.

To Do:
1) Add RBM
2) LSTM for site-specific regressions
3) ANN, LSTM, RBM sensitivity
4) Add more remote sensing layers
5) Choose best remote sensing spatial and temporal resolution
6) Add atmospheric CO2 data - in situ and remote sensing
7) Add AmeriFlux data to training/validation set
8) Check for obviusly spurious data (e.g., plots at bottom of <extract_fluxnet.m>)
9) Models based on IGPB classification - try the following
	9a) Model Averaging based on Bayesian performance and/or parameter similarity
	9b) Algorithm-specific model combinations (e.g., as discussed in Hinton's Coursera)
	9c) Spatial Kriging
	9d) Parameter Kriging
	9e) Different models for different IGBP classifications
	9f) Bring in MODIS annual IGBP map
10) Add DEM map to in situ and remote sensing
11) Plot local stats in k-fold, loo, sensitivity
12) Add k-fold to the sensitivity regressions

