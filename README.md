# SpatCorrNPR: A MATLAB and Octave Toolbox for Nonparametric Regression with Correlated Errors

Local polynomial regression is a widely used smoothing method in nonparametric statistics, and the main tuning parameter of this method is the bandwidth, which is susceptible to correlation. Several techniques exist for bandwidth selection with uncorrelated errors and even fewer for the correlated errors scenario; however, the availability of these implementations is meager and, in some cases, nonexistent. To address this gap, we present a new MATLAB and Octave toolbox, called SpatCorrNPR, which is compatible with Windows and Linux, for local polynomial regression estimation with correlated errors. The SpatCorrNPR toolbox supports one- and two-dimensional regression problems, featuring several bandwidth selection methods designed to mitigate the effects of correlation. Additionally, SpatCorrNPR offers covariance estimation capabilities, including the construction of semivariograms and the generation of geoplot maps for spatial datasets. This toolbox can be accessed in three forms: a MATLAB graphical user interface (GUI), a MATLAB command line toolbox, and an Octave toolbox. The toolbox is written in MATLAB to take advantage of the efficient computing system, and the Octave counterpart increases its availability to the user. We also demonstrate the utility of this toolbox with an illustrative example, and further compare the proposed toolbox to other software implementations using some practical datasets.

You can find the article at:

Please cite the SpatCorrNPR toolbox as:

