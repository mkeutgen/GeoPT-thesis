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

### Pr Weltje algorithm

So far, I've looked at a subproblem : major elements.
Pr Weltje proposes a unifying algorithm where instead of looking at a subcomposition we look at the total composition (summing to one) seen as :

1. Mass fraction of majors elements
2. Mass fraction of trace elements 
3. Mass fraction of oxygen
4. Mass fraction of the rest

In this way, we get feasible compositions from a physical point of view. How to handle missing values ? Pr Weltje proposes to take the log(x/1-x) which is a mapping from the Simplex to the Real Space. Pr Gijbels underlines that doing so creates dependance between the labs since we replace missing values with means.

### Principal Component Analysis 

Both professors were happy with the fact that we need a high number of PC's to explain the variability in the dataset. This suggests that the assumption of a multivariate standard normal distribution for the standardized clr data holds, all we have is random multivariate error.

If we analyze rock bodies, we would expect a few principal component analysis would be able to summarize the rock samples but in the case where we are looking at i.i.d analyses of the same rock we need a high number of PC's.

### Clusters, row-wise & column-wise missing values
For each rock, count the number of missing values row-wise & column-wise.
If number of missing values per column is larger than cut-off (let's say 2/3), we say all the column is missing.

One can also count missing values in each rows and if too many missing values in a lab, also throw them out.

## Meeting 14 February 10 am to 12 am 
### Comments regarding my work so far
The most important remark is that there is an issue in the algorithm. What I thought to be multiple interactions is only one iteration according to Pr Weltje.
What I should do is : 
1. Impute the missing values for the whole dataframe by the blr mean of each column.
2. If sum is larger than one, scale back to 1.
3. Delete the imputed values and replace them by the new blr means of the scaled dataframe.
4. If sum is larger than one, scale back to 1.
5. Delete the imputed values and replace them by the new blr means of the scaled dataframe.
6. Repeat until convergence.

To highlight in the thesis : BLR are sesnsitive to the size. They're equal to the alr transformation in the psecific case that they sum to one.
BLR are a way to exploit the scale of each observations in the dataset by using the constraint X + OX + Tr + R = 1. First everybody tried to get rid of the

Pr Gijbels : the normal distribution assumption is still an assumption. We should check if this assumption is reasonable through diagnostic plots, 2D histograms, etc...

Pr Gijbels : is PCA in the imputed dataframe similar to PCA in the NAs full dataframe ? 

We actually expect to have less outliers in the imputed df than in the original df.   

Overall, both professors were satisfied with what I did so far, Pr Weltje said I was on the right track. 
## Practicalities
Intermediate Presentation will be in the week of 21st of March.
What I'm doing ? 
What I did ?
Where I'm going ?

Illustration example
Algorithm
Show results.

Next meeting : I can also meet individually with Pr. Gijbels or Weltje to discuss specifics. (algorithm with Weltje for example).

## Meeting 2nd of March 10 a.m
# Still to do :
The outlier detection function which I implemented has a major drawback, it is univariate. Pr Weltje suggests instead to do a principal component analysis on multivariate normal data and compare it with PCA on my datasets. Under the null hypothesis that the generating process at the origin of my data is multivariate normal, one would need a very high number of dimensions. If correlated datasets, a lesser number of Principal Components and hence of dimensions would be sufficient to reconstruct the data matrix. How does the reconstruction error increases with the size of the dataset ? 

IF I am able to answer to this question, I shall be able to reconstruct the data matrix and see which are the observations which are contributing the most to this reconstruction error. Reconstruction error is the sum of squared distances between actual data matrix and its projection in the R space where R is the number of meaningful dimensions which was found by PCA. => Goal of estimating the actual rank of the data matrices.

IF analytical methods are unbiased, the assumption of multivariate normal distribution would be satisfied. Apparently analysis is biased and therefore we must have some kind of methods for outlier detection to remove dubious results from some laboratories.

Then, once outlier detection is done based on this PCA-based approach, we move to the next step and we summarize each rock to its mean and its variance covariance matrix. We further summarize them in a few dimensions only, the dimensions which were considered as meaningful through the PCA. 

! These dimensions should be the same throughout different rocks !


We then test the null hypothesis that covariance-variance matrices are the same, that would mean the biases have been the same, there is some systematic error going on rounds and rounds.

In the end, we generate random variables in the euclidean space which are normally distributed and then we transform them through ilr.inv for example and see what it would mean for the chemist.





