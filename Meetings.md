# Summary of meetings with adivsors Pr. Gijbels & Pr. Weltje

## Meeting 6th December 10am to 12 am

The main issue right now is the number of outliers. It appears that they have been overestimated by the robust methods. Flagging as outliers approximately 40 out of 90 observations is not reasonable. Pr Gijbels stresses the importance of knowing whether there are 6 or 40 outliers. Pr Weltje suggests of colouring observations on plot "Robust distance vs Index" the observations which were imputed with the geometric mean.

The approach I've taken so far could seriously go wrong if it appears that there was an artificial generation of outliers by filling missing values with geometric mean and therefore artificially generating a lot of identical values. Departures from these identical values will be wrongly assumed to be outliers while they are only data variability

As a rule of thumb, Pr Weltje suggests that no more than 5 or 10 % of observations should be flagged as outliers.

Advice : once you've identified the outliers redo the analysis without outliers and project the outliers in the PC space estimated without outliers.

Advice : see if classical methods conducted on the bulk of the data (without outliers) and robust method (with outliers) yield the same results.

Pr Weltje's comment on missing values : either the elements sum to less than one and therefore 1 - sum(elements) = rest either this sum is more than one and then the rest is a missing value.

Plot the average of values recorded for each element against the % of labs which have discovered the value (ex : 10,80%) for CaO means that 80 % of labs have reported a value whose average is 10 unit.
We can look at log(x/1-x) which is a marginal distribution. This approach automatically takes care of skewness of distribution.

Professors have recalled that there should be two approaches to handle missing values : impute them by geometric mean for major elements and use the EM-algorithm for non missing values.

Biplot : replace iteratively missing values by a value smaller and smaller. 

## Meeting 20 December 10 am to 12 am

^
