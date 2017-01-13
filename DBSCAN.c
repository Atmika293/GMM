#include "DBSCAN.h"

int getNeighbourhoodMatrix(Point *data, int data_length, enum metric m,
                            double distance_threshold, int min_points,
                            int **neighbours, int *valid_clusters)
{
    int i  = 0, j = 0, neighbours_count = 0, retVal = 0;
    int *row;
    double dist;

    neighbours = (int**)malloc(sizeof(int*) * data_length);
    valid_clusters = (int*)malloc(sizeof(int) * data_length);

    for(i = 0;i < data_length;i++)
    {
        row = (int*)malloc(sizeof(int) * data_length);
        neighbours_count = 0;
        for(j = 0;j < data_length;j++)
        {
            if(i != j)
            {
                dist = distance(data[i].features, data[j].features, m);
                if(dist < 0)
                {
                    retVal = -1;
                    break;
                }

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
        if(retVal != -1)
        {
            if(neighbours_count >= min_points)
            {
                neighbours[i] = row;
                valid_clusters[i] = 1;
            }
            else
            {
                neighbours[i] = NULL;
                valid_clusters[i] = 0;
                free(row);
            }
        }
        else
        {
            free(row);
            break;
        }

    }

    return retVal;
}

void expandCluster(int pointIdx, Point *data, int data_length,
                   int **neighbours, int *valid_clusters,
                   int clusterID)
{
    int neighbourIdx = 0;
    if(valid_clusters[pointIdx]) // Not a noise point
    {
        for(neighbourIdx = 0;neighbourIdx < data_length;neighbourIdx++)
        {

            if(neighbours[pointIdx][neighbourIdx]) //traversing through neighbouring points
            {
                if(data[neighbourIdx].visited == 0)
                {
                    data[neighbourIdx].visited = 1;
                    expandCluster(neighbourIdx, data, data_length,
                                  neighbours, valid_clusters, clusterID);
                }

                if(data[neighbourIdx].cluster <= 0) //if the neighbour has been visited and not assigned
                    data[neighbourIdx].cluster = clusterID; //or marked a noise point
            }
        }
    }
    else
    {
        data[pointIdx].visited = 1;
        data[pointIdx].cluster = clusterID;
    }
}

int DBSCAN(Point *data, int data_length, enum metric m,
            double distance_threshold, int min_points)
{
    int i = 0;
    int **neighbours;
    int *valid_clusters;
    int clusterID = 0;

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







