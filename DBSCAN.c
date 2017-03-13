#include "DBSCAN.h"

int getNeighbourhoodMatrix(Point *data, int data_length, enum metric m,
                            double distance_threshold, int min_points,
                            int **neighbours, int *valid_clusters)
{
    int i  = 0, j = 0, neighbours_count = 0;
    double dist;

    for(i = 0;i < data_length;i++)
    {
        neighbours_count = 0;
        for(j = 0;j < data_length;j++)
        {
            if(i != j)
            {
                dist = distance(data[i].features, data[j].features, m);
                if(dist < 0)
                    return -1;

                if(dist <= distance_threshold)
                {
                    neighbours[i][j] = 1;
                    neighbours_count += 1;
                }
                else
                    neighbours[i][j] = 0;
            }
            else
                neighbours[i][j] = 1;
        }

        if(neighbours_count >= min_points)
            valid_clusters[i] = 1;
        else
            valid_clusters[i] = 0;
    }

    return 0;
}

void expandCluster(int pointIdx, Point *data, int data_length,
                   int **neighbours, int *valid_clusters,
                   int clusterID)
{
    int neighbourIdx = 0;

    data->cluster = clusterID;

    for(neighbourIdx = 0;neighbourIdx < data_length;neighbourIdx++)
    {
        if(neighbours[pointIdx][neighbourIdx]) //traversing through neighbouring points
        {
            if(data[neighbourIdx].visited == 0)
            {
                data[neighbourIdx].visited = 1;
                expandCluster(neighbourIdx, data, data_length,
                              neighbours, cluster_point, clusterID);
            }

            if(data[neighbourIdx].cluster <= 0) //if the neighbour has been visited and not assigned
                data[neighbourIdx].cluster = clusterID; //or marked a noise point
        }
    }
}

int DBSCAN(Point *data, int data_length, enum metric m,
            double distance_threshold, int min_points)
{
    int i = 0;
    int **neighbours;
    int *valid_clusters;
    int clusterID = 0;

    neighbours = (int**)malloc(sizeof(int*) * data_length);
    for(i = 0;i < data_length;i++)
        neighbours[i] = (int*)malloc(sizeof(int) * data_length);
    valid_clusters = (int*)malloc(sizeof(int) * data_length);

    if(getNeighbourhoodMatrix(data, data_length, m,
                           distance_threshold, min_points,
                           neighbours, valid_clusters) != -1)
    {
        for(i = 0;i < data_length;i++)
        {
            if(data[i].visited)
                continue;

            data[i].visited = 1;
            if(valid_clusters[i])
            {
                clusterID += 1;
                expandCluster(i, data, data_length, neighbours, valid_clusters, clusterID);
            }
            else
                data[i].cluster = 0;//NOISE
        }

    }

    free(valid_clusters);
    for(i = 0;i < data_length;i++)
        free(neighbours[i]);
    free(neighbours);

    return clusterID;
}

int DBSCANCluster(double **points, int data_length, int dim_length,
                    double distance_threshold, int min_points)
{
    int i = 0;
    Point *data = (Point*)malloc(sizeof(Point) * data_length);
    for(i = 0;i < data_length;i++)
        initialize_point(&(data[i]), points[i], dim_length);

    return DBSCAN(data, data_length, EUCLIDEAN, distance_threshold, min_points);
}







