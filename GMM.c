#include "GMM.h"

#include <string.h>
#include <math.h>

void initGMM(GaussianMixtureModel *gmm, int maxIterations,
             double deltaLogLikelihood)
{
    gmm->maxIterations = maxIterations;
    gmm->deltaLogLikelihood = deltaLogLikelihood;
    gmm->previouslogLikelihood = 0;
}

int clusterData(Point *data, int data_length,
                 enum metric m, double distance_threshold, int min_points,
                 GaussianMixtureModel *gmm)
{
    int row = 0, k = 0, retVal = 0;
    gmm->clusterCount = DBSCAN(data, data_length, m,
                    distance_threshold, min_points);
    if(gmm->clusterCount > 0)
    {
        gmm->mean = (Point*)calloc(sizeof(Point), gmm->clusterCount);
        gmm->variance = (double*)calloc(sizeof(double), gmm->clusterCount);
        gmm->mixtureWeights = (double*)calloc(sizeof(double), gmm->clusterCount);
        int *clusterStrength = (int*)calloc(sizeof(int), gmm->clusterCount);
        double **membershipWeights = (double**)malloc(sizeof(double*) * data_length);
        for(row = 0;row < data_length;row++)
        {
            membershipWeights[row] = (double*)malloc(sizeof(double*) * gmm->clusterCount);
            memset(membershipWeights[row], 1, gmm->clusterCount);

            clusterStrength[data[row].cluster - 1] += 1;
        }

        for(k = 0;k < gmm->clusterCount;k++)
            gmm->mixtureWeights[k] = ((double)clusterStrength[k]) /
                                     ((double)data_length);

        retVal = calculateMeanAndVariance(data, data_length, gmm->clusterCount, membershipWeights,
                                          gmm->mean, gmm->variance);

        for(row = 0;row < data_length;row++)
            free(membershipWeights[row]);

        free(membershipWeights);
        free(clusterStrength);
    }
    else
        retVal = -1;

    return retVal;
}

double calculateGaussianProbability(Point x, Point mean, double variance)
{
    if(x.features.length != mean.features.length)
        return -1;

    double pi = acos(-1);
    int f = 0;
    double avg = 0, var = 0;
    int dim = x.features.length;

    for(f = 0;f < dim;f++)
    {
        avg = mean.features.vector[f];
        var += ((x.features.vector[f] - avg) * (x.features.vector[f] - avg));
    }

    return exp(-0.5 * (var/ variance))/sqrt(pow(2 * pi, dim) * variance);
}

int isConverged(Point *data, int data_length,
                 GaussianMixtureModel *gmm)
{
    int row = 0, k = 0;
    int dim = data[0].features.length;
    Point p;
    double likelihood = 0;
    double logLikelihood = 0;

    for(row = 0;row < data_length;row++)
    {
        p = data[row];

        if(p.features.length != dim)
            return -1;

        likelihood = 0;

        for(k = 0;k < gmm->clusterCount;k++)
            likelihood += gmm->mixtureWeights[k] * calculateGaussianProbability(p, gmm->mean[k], gmm->variance[k]);

        logLikelihood += log(likelihood);
    }

    if((logLikelihood - gmm->previouslogLikelihood) <= gmm->deltaLogLikelihood)
        return 1;
    else
        gmm->previouslogLikelihood = logLikelihood;

    return 0;
}

int expectation(Point *data, int data_length,
                  GaussianMixtureModel *gmm,
                  double **membershipWeights)
{
    int row = 0, k = 0, retVal = 0;
    double *likelihoods = (double*)malloc(sizeof(double) * gmm->clusterCount);
    double sumlikelihood = 0;
    double probability = 0;
    Point p;

    for(row = 0;row < data_length;row++)
    {
        sumlikelihood = 0;
        p = data[row];

        for(k = 0;k < gmm->clusterCount;k++)
        {
            probability = calculateGaussianProbability(p, gmm->mean[k], gmm->variance[k]);
            if(probability == -1)
            {
                retVal = -1;
                break;
            }
            likelihoods[k] = gmm->mixtureWeights[k] * probability;
            sumlikelihood += likelihoods[k];
        }

        if(retVal != -1)
        {
            for(k = 0;k < gmm->clusterCount;k++)
                membershipWeights[row][k] = likelihoods[k] / sumlikelihood;
        }
        else
            break;
    }

    free(likelihoods);
    return retVal;
}

int maximization(Point *data, int data_length,
                 double **membershipWeights,
                 GaussianMixtureModel *gmm)
 {
     int row = 0, k = 0;
     double clusterStrength = 0;

     for(k = 0;k < gmm->clusterCount;k++)
     {
         clusterStrength = 0;

         for(row = 0;row < data_length;row++)
             clusterStrength += membershipWeights[row][k];

         gmm->mixtureWeights[k] = clusterStrength / data_length;
     }

     return calculateMeanAndVariance(data, data_length, gmm->clusterCount,
                                     membershipWeights,
                                     gmm->mean, gmm->variance);

 }

int EMWithDBSCAN(Point *data, int data_length,
                 enum metric m, double distance_threshold,
                 int min_points, GaussianMixtureModel *gmm,
                 double **membershipWeights)
 {
     int iter = 0;

     if(clusterData(data, data_length, m, distance_threshold, min_points, gmm) != -1)
     {
         for(iter = 0;iter < gmm->maxIterations;iter++)
         {
            if(expectation(data, data_length, gmm, membershipWeights))
            {
                if(maximization(data, data_length, membershipWeights, gmm))
                {
                    if(isConverged(data, data_length, gmm))
                        break;
                    else
                        continue;
                }
            }

            return -1;
         }
     }

     return 0;
 }

int EMStepE(Point *data, int data_length, int clusterCount,
            Point *initial_mean, double *initial_variance,
            double *initial_mixtureWeights, GaussianMixtureModel *gmm,
            double **membershipWeights)
 {
     int iter = 0, k = 0;
     gmm->clusterCount = clusterCount;
     gmm->mean = (Point*)calloc(sizeof(Point), clusterCount);
     gmm->variance = (double*)calloc(sizeof(double), clusterCount);
     gmm->mixtureWeights = (double*)calloc(sizeof(double), clusterCount);
     for(k = 0;k < clusterCount;k++)
     {
         gmm->mean[k] = initial_mean[k];
         gmm->variance[k] = initial_variance[k];
         gmm->mixtureWeights[k] = initial_mixtureWeights[k];
     }

     for(iter = 0;iter < gmm->maxIterations;iter++)
     {
         if(expectation(data, data_length, gmm, membershipWeights))
         {
             if(maximization(data, data_length, membershipWeights, gmm))
             {
                 if(isConverged(data, data_length, gmm))
                     break;
                 else
                     continue;
             }
         }

         return -1;
     }

     return 0;
 }

 int EMStepM(Point *data, int data_length,
              GaussianMixtureModel *gmm,
              double **membershipWeights)
 {
     int iter = 0;

     for(iter = 0;iter < gmm->maxIterations;iter++)
     {
        if(maximization(data, data_length, membershipWeights, gmm))
        {
            if(expectation(data, data_length, gmm, membershipWeights))
            {
                if(isConverged(data, data_length, gmm))
                    break;
                else
                    continue;
            }
        }

        return -1;
     }

     return 0;
 }

 int ClusterEM(double **points, int data_length, int dim_length,
                enum metric m, double distance_threshold,
                int min_points, GaussianMixtureModel *gmm,
               double **membershipWeights)
 {
     int i = 0;
     Point *data = (Point*)malloc(sizeof(Point) * data_length);
     for(i = 0;i < data_length;i++)
         initialize_point(&(data[i]), points[i], dim_length);

     return EMWithDBSCAN(data, data_length, m, distance_threshold,
                         min_points, gmm, membershipWeights);
 }

 int ClusterEMStepE(double **points, int data_length, int dim_length,
                     int clusterCount, double **initial_mean, double *initial_variance,
                     double *initial_mixtureWeights, GaussianMixtureModel *gmm,
                    double **membershipWeights)
 {
     int i = 0;
     Point *data = (Point*)malloc(sizeof(Point) * data_length);
     Point *initial_center = (Point*)malloc(sizeof(Point) * data_length);
     for(i = 0;i < data_length;i++)
         initialize_point(&(data[i]), points[i], dim_length);
     for(i = 0;i < clusterCount;i++)
         initialize_point(&(initial_center[i]), initial_mean[i], dim_length);


     return EMStepE(data, data_length, clusterCount, initial_center,
             initial_variance, initial_mixtureWeights, gmm, membershipWeights);
 }

 int ClusterEMStepM(double **points, int data_length, int dim_length,
                     GaussianMixtureModel *gmm,
                     double **membershipWeights)
 {
     int i = 0;
     Point *data = (Point*)malloc(sizeof(Point) * data_length);
     for(i = 0;i < data_length;i++)
         initialize_point(&(data[i]), points[i], dim_length);

     return EMStepM(data, data_length, gmm, membershipWeights);
 }
