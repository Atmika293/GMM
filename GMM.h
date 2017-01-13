#ifndef GMM_H
#define GMM_H

#include "DataPoints.h"
#include "DBSCAN.h"

typedef struct GaussianMixtureModel
{
    int clusterCount;
    int maxIterations;
    double deltaLogLikelihood;
    double previouslogLikelihood;
    Point *mean;
    double *variance;
    double *mixtureWeights;

}GaussianMixtureModel;

void initGMM(GaussianMixtureModel *gmm, int maxIterations,
             double deltaLogLikelihood);

int clusterData(Point *data, int data_length, enum metric m,
                 double distance_threshold, int min_points,
                 GaussianMixtureModel *gmm);

double calculateGaussianProbability(Point x, Point mean,
                                double variance);

int isConverged(Point *data, int data_length,
                GaussianMixtureModel *gmm);

int expectation(Point *data, int data_length,
                 GaussianMixtureModel *gmm,
                 double **membershipWeights);

int maximization(Point *data, int data_length,
                  double **membershipWeights,
                  GaussianMixtureModel *gmm);

int EMWithDBSCAN(Point *data, int data_length,
                  enum metric m, double distance_threshold,
                  int min_points, GaussianMixtureModel *gmm,
                 double **membershipWeights);

int EMStepE(Point *data, int data_length, int clusterCount,
            Point *initial_mean, double *initial_variance,
            double *initial_mixtureWeights, GaussianMixtureModel *gmm,
            double **membershipWeights);

int EMStepM(Point *data, int data_length,
            GaussianMixtureModel *gmm,
            double **membershipWeights);

int ClusterEM(double **points, int data_length, int dim_length,
               enum metric m, double distance_threshold,
               int min_points, GaussianMixtureModel *gmm,
               double **membershipWeights);

int ClusterEMStepE(double **points, int data_length, int dim_length, int clusterCount,
                    double **initial_mean, double *initial_variance,
                    double *initial_mixtureWeights, GaussianMixtureModel *gmm,
                   double **membershipWeights);

int ClusterEMStepM(double **points, int data_length, int dim_length,
                    GaussianMixtureModel *gmm,
                   double **membershipWeights);




#endif // GMM_H
