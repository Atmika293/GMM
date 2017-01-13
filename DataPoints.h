#ifndef DATAPOINTS_H
#define DATAPOINTS_H

#include <stdlib.h>

typedef struct PointFeatures
{
    double *vector;
    int length;

}PointFeatures;

typedef struct Point
{
    PointFeatures features;
    int visited;
    int cluster;
}Point;

enum metric
{
    EUCLIDEAN = 0,
    MANHATTAN,
    CHEBYSHEV
};

void initialize_point(Point *p, double *features, int length);

double distance(PointFeatures f1, PointFeatures f2, enum metric m);

int calculateMeanAndVariance(Point *data, int data_length,
                              int clusterCount, double **membership_weights, 
                                     Point *mean, double *variance);

#endif // DATAPOINTS_H
