# GMM
Gaussian Mixture Model based on Expectation Maximization and DBSCAN

A Gaussian Mixture Model (GMM) is a parametric probability density function represented as a weighted sum of Gaussian
component densities, each with its own set of parameters : a mean vector, variance and a mixture weight. GMM computes the probabilties of a data point belonging to each cluster; the point is not assigned to any partcular cluster.

The parameters and output probabilities are computed using the Expectation Maximization(EM) algorithm. The algorithm operates on the assumption that the number of clusters is fixed and known. A variety of clustering methods may be used to estimate the number of clusters and the initial parameters of each cluster's distribution. This project uses DBSCAN.

### Functions ###
The main functions to use in the project are:
- void initGMM(GaussianMixtureModel *gmm, int iterations,
             double deltaLogLikelihood);
             
- int ClusterEM(double **points, int data_length, int dim_length,
               enum metric m, double distance_threshold,
               int min_points, GaussianMixtureModel *gmm,
               double **membershipWeights);

- int ClusterEMStepE(double **points, int data_length, int dim_length, int clusterCount,
                    double **initial_mean, double *initial_variance,
                    double *initial_mixtureWeights, GaussianMixtureModel *gmm,
                   double **membershipWeights);

- int ClusterEMStepM(double **points, int data_length, int dim_length,
                    GaussianMixtureModel *gmm,
                   double **membershipWeights);

### Function Parameters ###
- points is m x n matrix, where m is the number of data points and n is the number of dimensions of each point.
data_length = m, dim_length = n.

- metric is the method used to compute the distance between two points (used in DBSCAN). The options are EUCLIDEAN,
MANHATTAN and CHEBYSHEV.

- distance_threshold is the maximum distance/dissimilarity between two neighbouring points (used in DBSCAN). Neighbouring points belong to the same cluster. This property is transitive.

- min_points is the minimum number of points that must be in a mutual neighbourhood for them to form a cluster (used in DBSCAN).

- GaussianMixtureModel is a structure that carries information about the model including maximum number of iterations to be carried out for the EM algorithm (maxIterations) and deltaLogLikelihood, which is the measure of convergence. The EM algorithm terminates when the number of iterations exceeds maxIterations or the relative change in likelihood logarithm is lesser than deltaLogLokelihood. The structure also carries information about cluster parameters, such as mean, variance and mixture weight.
membershipWeights is a m x K matrix that stores the membership weight of the ith point in the kth cluster.
clusterCount (K) is the number of clusters in the model. 

### Avoiding Initialization ###
Initialization through DBSCAN may be avoided by calling either ClusterEMStepE and ClusterEMStepM. 

ClusterEMStepE starts with the Expectation step and requires initial_mean (K x n matrix of K points), initial_variance (vector of length K) and initial_mixtureWeights (vector of length K) as input.

ClusterEMStepM starts with the Maximization step and requires membershipWeights (m x K) as input. The function modifies this matrix itself to produce the output.




