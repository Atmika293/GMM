#ifndef DBSCAN_H
#define DBSCAN_H

#include "DataPoints.h"

int getNeighbourhoodMatrix(Point *data, int data_length, enum metric m,
                            double distance_threshold, int min_points, int **neighbours,
                            int *valid_clusters);

void expandCluster(int pointIdx, Point *data, int data_length,
                   int **neighbours, int *valid_clusters,
                   int clusterID);

int DBSCAN(Point *data, int data_length, enum metric m,
            double distance_threshold, int min_points);

int DBSCANCluster(double **points, int data_length, int dim_length,
           double distance_threshold, int min_points);

#endif // DBSCAN_H
