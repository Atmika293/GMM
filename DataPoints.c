#include "DataPoints.h"

#include <string.h>
#include <math.h>

void initialize_point(Point *p, double *features, int length)
{
    memcpy(p->features.vector, features, length);
    p->cluster = -1; //UNMARKED
    p->visited = 0; //NOT VISITED
}

double distance(PointFeatures f1, PointFeatures f2, enum metric m)
{
    if(f1.length != f2.length)
        return -1;

    int i = 0;
    double sum = 0, max = -1, diff = 0;
    int length = f1.length;

    switch(m)
    {
        case EUCLIDEAN:
        {
            sum = 0;
            for(i = 0;i < length;i++)
                sum += pow((f1.vector[i] - f2.vector[i]), 2);
            return sqrt(sum);
        }

        case MANHATTAN:
        {
            sum = 0;
            for(i = 0;i < length;i++)
                sum += (f1.vector[i] < f2.vector[i])? (f2.vector[i] - f1.vector[i]) : (f1.vector[i] - f2.vector[i]);
            return sum;
        }

        case CHEBYSHEV:
        {
            max = -1;
            for(i = 0;i < length;i++)
            {
                diff = (f1.vector[i] < f2.vector[i])? (f2.vector[i] - f1.vector[i]) : (f1.vector[i] - f2.vector[i]);
                if(diff > max)
                    max = diff;
            }
            return max;
        }
    }

    return -1;
}

int calculateMeanAndVariance(Point *data, int data_length,
                              int clusterCount, double **membership_weights,
                              Point *mean, double *variance)
{
    int retVal = 0, row = 0, f = 0, k = 0;
    double var = 0, avg = 0, mw = 0;
    Point p;
    int dim = data[0].features.length;
    int *clusterStrength = (int*)calloc(sizeof(int), clusterCount);

    for(row = 0;row < data_length;row++)
    {
        p = data[row];
        mw = membership_weights[row][p.cluster - 1];
        clusterStrength[p.cluster - 1] += 1;

        if(dim != p.features.length)
        {
            retVal = -1;
            break;
        }

        for(f = 0;f < dim;f++)
            mean[p.cluster - 1].features.vector[f] += mw * p.features.vector[f];
    }

    if(retVal != -1)
    {
        for(k = 0;k < clusterCount;k++)
        {
            for(f = 0;f < dim;f++)
                mean[k].features.vector[f] = mean[k].features.vector[f]/clusterStrength[k];
        }

        for(row = 0;row < data_length;row++)
        {
            p = data[row];
            mw = membership_weights[row][p.cluster - 1];

            for(f = 0;f < dim;f++)
            {
                if(dim != mean[p.cluster - 1].features.length)
                    return -1;

                avg = mean[p.cluster - 1].features.vector[f];

                var = (p.features.vector[f] - avg) * (p.features.vector[f] - avg);

                variance[p.cluster - 1] += mw * var;
            }
        }

        for(k = 0;k < clusterCount;k++)
            variance[k] = variance[k]/clusterStrength[k];
    }

    free(clusterStrength);
    return retVal;
}
