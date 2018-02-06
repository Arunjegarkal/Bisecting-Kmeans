#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <sys/time.h>
/*Function to find second cluster's centroid*/
void findCentroid_2(int *data1, int **data2, int dim,int size,int current_cluster_size);
/*Function to find i'th cluster's centroid*/
void findcentroid(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int ** temp_centroid);
/*Function to create cluster using bisecting Kmeans algorithm*/
void BisectingKmeans(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *searchData,int *clusterstart);
/*Function to assign data points to cluster*/
void ClusterAssign(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *temp_array);
/*Function to search a datapoint in a cluster*/
void Search(int dim,int ndata,int *data,int *searchData,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *clusterstart);
/*Function to clear Cluster size, cluster assign and cluster radious before calculating new centroids*/
void initializeZero(int ndata,int k,int *clustersize,int **cluster_assign,double *cluster_raidus);
/*Function to find Calculate SSE*/
int findSSE(int dim,int ndata,int *data,int k,int ** cluster_centroid, int **cluster_assign,int *clustersize,int cluster);
/*Function to calculate the centroid until stable midpoints are formed*/
void calculateCentroid(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *temp_array,int ** temp_centroid);
/*Sort the data array according to cluster assign*/
void sort(int dim,int ndata,int *data,int k,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *clusterstart);
void BisectingKmeans(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *searchData,int *clusterstart)
{
    int i=0,y,j,max_sseCluster=0,temp=0,loop;
    int x=((rand()%ndata));
    y=x-(x%dim);
    int **temp_centroid;
    int *temp_array=(int*)malloc(ndata * sizeof(int));
    temp_centroid=(int **)malloc(k * sizeof(int *));
    for (i=0; i<k; i++)
         *(temp_centroid+i) = (int *)malloc(dim * sizeof(int));
    for(i=0;i<k;i++)
    {
        for(j=0;j<dim;j++)
        {
            temp_centroid[i][j]=0;
        }
    }
    //Find cluster1 centroids
    for(i=0;i<dim;i++)
    {
        cluster_centroid[0][i]=data[y+i];
    }
    //Find centroid 2
    findCentroid_2(data,cluster_centroid,dim,ndata,1);
    //ClusterAssign function call to assign all the data points to assign to recently calculated clusters
    ClusterAssign(dim,ndata,data,2,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_array);
    int flag=0;
    //Loop untill K clusters are calculated
    for(loop=2;loop<k;loop++)
    {
        //printf("Calculating Clusters...\n");
        calculateCentroid(dim,ndata,data,loop,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_array,temp_centroid);
        //Calculate the sse for each clusters
        for(i=0;i<loop;i++)
        {
            int sse=findSSE(dim,ndata,data,loop,cluster_centroid,cluster_assign,clustersize,i);
            if(sse>temp)
                max_sseCluster=i;
            temp=sse;
        }
        int Clust1_sum[dim],Clust2_sum[dim],m,n;
        for(m=0;m<dim;m++)
        {
            Clust1_sum[m]=0;
            Clust2_sum[m]=0;
        }
        temp=0;
        for(j=0;j<1;j++)
        {
            for(i=0;i<clustersize[max_sseCluster];i++)
                if(data[cluster_assign[max_sseCluster][i]+j]>cluster_centroid[max_sseCluster][j])
                {
                    for(m=0;m<dim;m++)
                    {
                        Clust1_sum[m]+=data[cluster_assign[max_sseCluster][i]+m];
                    }
                }
                else
                {
                    for(m=0;m<dim;m++)
                    {
                        Clust2_sum[m]+=data[cluster_assign[max_sseCluster][i]+m];
                    }
                }
        }
        for(m=0;m<dim;m++)
        {
            Clust1_sum[m]=Clust1_sum[m]/clustersize[max_sseCluster];
            Clust2_sum[m]=Clust2_sum[m]/clustersize[max_sseCluster];
        }
        //add new centroids
        for(i=0;i<dim;i++)
        {
            cluster_centroid[max_sseCluster][i]=Clust1_sum[i];
            cluster_centroid[loop][i]=Clust2_sum[i];
        }
        if(loop+1<k)
            initializeZero(ndata,loop,clustersize,cluster_assign,cluster_raidus);
    }
    calculateCentroid(dim,ndata,data,loop,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_array,temp_centroid);
    printf("\nCluster Centroids\n");
    for(i=0;i<k;i++)
    {
        printf("Cluster %d  ",i);
        for(j=0;j<dim;j++)
        {
            printf("[%d]",cluster_centroid[i][j]);
        }
        printf("\n");
    }
    for(i=0;i<k;i++)
    {
        printf("Size of cluster %d is %d\n",i,clustersize[i]);
    }
    sort(dim,ndata,data,k,cluster_centroid,cluster_assign,clustersize,clusterstart);
    for(i=0;i<k;i++)
    {
        printf("Cluster start of cluster %d is %d\n",i,clusterstart[i]);
    }
    Search(dim,ndata,data,searchData,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,clusterstart);
}
void calculateCentroid(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *temp_array,int ** temp_centroid)
{
    int flag=0,i,j;
    //Repeat loop until clusters centroids are formed
    do
    {
        flag=0;
        initializeZero(ndata,k,clustersize,cluster_assign,cluster_raidus);
        ClusterAssign(dim,ndata,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_array);
        findcentroid(dim,ndata,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,temp_centroid);
        for(i=0;i<k;i++)
        {
            for(j=0;j<dim;j++)
            {
                if(cluster_centroid[i][j]!=temp_centroid[i][j])
                    flag=1;
            }
        }
        for(i=0;i<k;i++)
        {
            for(j=0;j<dim;j++)
            {
                cluster_centroid[i][j]=temp_centroid[i][j];
            }
        }
    }while(flag==1);

}

void sort(int dim,int ndata,int *data,int k,int ** cluster_centroid, int **cluster_assign,int *clustersize,int * clusterstart)
{
    int i,j,l,count=0;
    int *temp_array = (int*)malloc(ndata * sizeof(int));
    for(i=0;i<k;i++)
    {
        clusterstart[i]=count/dim;
        for(j=0;j<clustersize[i];j++)
        {
            for(l=0;l<dim;l++)
            {
                temp_array[count]=data[cluster_assign[i][j]+l];
                count++;
            }
        }
    }data=temp_array;
}
/*Function return SSE of the cluster*/
int findSSE(int dim,int ndata,int *data,int k,int ** cluster_centroid, int **cluster_assign,int *clustersize,int cluster)
{
    int i,j,sum_dim=0,temp=0,ssetemp=0;
    int mean_array[dim],sse[dim];
    //calculate mean
    for(j=0;j<dim;j++)
    {
        temp=0;
        for(i=0;i<clustersize[cluster];i++)
        {
            temp+=data[cluster_assign[cluster][i]+j];
        }
        mean_array[j]=temp/clustersize[cluster];
    }
    for(j=0;j<dim;j++)
    {
        temp=0;
        for(i=0;i<clustersize[cluster];i++)
        {
            temp+=(data[cluster_assign[cluster][i]+j]-mean_array[j])*(data[cluster_assign[cluster][i]+j]-mean_array[j]);
        }
        ssetemp+=temp;
    }return ssetemp;
}
void Search(int dim,int ndata,int *data,int *searchData,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *clusterstart)
{
    int i,ii,j,val1,val2,search_cluster[k],nearest_node,visited_nodes=0,search_cluster_size=0,startPos;
    double temp,dist,mindist=9999999999,search_cluster_dist[k];
    struct timeval stop, start;
    gettimeofday(&start, NULL);
	//Find the cluster in which search node fall in based on radius
    for(i=0;i<k;i++)
    {
        temp=0.0;
        for(j=0;j<dim;j++)
        {
            val1=cluster_centroid[i][j];
            val2=searchData[j];
            temp+=((val1-val2)*(val1-val2));
        }
        dist=sqrt(temp);
        //Store the cluster number if the search data point falls with in the radious of cluster
		if ((cluster_raidus[i]-dist)>0)
        {
            mindist=dist;
            search_cluster_dist[search_cluster_size]=mindist;
            search_cluster[search_cluster_size]=i;
            search_cluster_size++;
        }
        //visited_nodes++;
    }
    mindist=99999999999;
    //Search the nearest data point in all the clusers in which the search data point falls with in the radious of cluster
    for(ii=0;ii<search_cluster_size;ii++)
    {
        for(i=0;i<clustersize[search_cluster[ii]];i++)
        {
            temp=0.0;
            for(j=0;j<dim;j++)
            {
                val1=data[clusterstart[search_cluster[ii]]+i+j];
                val2=searchData[j];
                temp+=(val1-val2)*(val1-val2);
            }
            dist=sqrt(temp);
            if (mindist>dist)
            {
                mindist=dist;
                nearest_node=cluster_assign[search_cluster[ii]][i];
            }
            visited_nodes++;
        }
        if(((ii+1)<search_cluster_size) || (mindist<search_cluster_dist[ii+1]))
        {
            ii=ii+1;
        }
    }
    printf("\nSearch node is ");
    for(i=0;i<dim;i++)
    {
        printf("[%d] ",searchData[i]);
    }
    printf("\nNearest node is ");
    for(i=0;i<dim;i++)
    {
        printf("[%d] ",data[nearest_node+i]);
    }
    printf("Number of Nodes visited %d \nDistance is %f",visited_nodes,mindist);
}
void findcentroid(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int ** temp_centroid)
{
    int i,j,l,val=0;

    for(i=0;i<k;i++)
    {
        for(j=0;j<dim;j++)
        {
            temp_centroid[i][j]=0;
        }
    }
    /*Repeat for K clusters*/
    for(i=0;i<k;i++)
    {
        for(l=0;l<dim;l++)
        {
            for(j=0;j<clustersize[i];j++)
            {
                /*Calculate Sum of all data points in a cluster*/
                val=val+data[cluster_assign[i][j]+l];
            }
            if(val!=0)
            {
                /*Calculate mid point*/
                temp_centroid[i][l]=(val/clustersize[i]);
            }
            else
                temp_centroid[i][l]=0;
            val=0;
        }
    }
}
void ClusterAssign(int dim,int ndata,int *data,int k,double *cluster_raidus,int ** cluster_centroid, int **cluster_assign,int *clustersize,int *temp_array)
{

    int i=0,j,m,val1,val2;
    double temp,dist;
    for(j=0;j<ndata;j=j+dim)
    {
        double min_dist=99999999999;
        for(m=0;m<k;m++)
        {
            /*Loop to calculate the distance*/
            temp=0;
            for(i=0;i<dim;i++)
            {
                val1=data[j+i];
                val2=cluster_centroid[m][i];
                temp+=(val1-val2)*(val1-val2);
            }
            dist=sqrt(temp);
            /*store the distance if the calculated distance is less than previously stored distance*/
            if(min_dist>dist)
            {
                min_dist=dist;
                temp_array[j]=m;
                if(cluster_raidus[m]<min_dist)
                {
                    cluster_raidus[m]=min_dist;
                }
            }
        }
    }
    for(j=0;j<ndata;j=j+dim)
    {
        cluster_assign[temp_array[j]][clustersize[temp_array[j]]]=j;
        clustersize[temp_array[j]]=clustersize[temp_array[j]]+1;
    }
}
void initializeZero(int ndata,int k,int *clustersize,int **cluster_assign,double *cluster_raidus)
{
    int i,j;
    for(i=0;i<k;i++)
    {
        clustersize[i]=0;
        cluster_raidus[i]=0;
    }
    for(i=0;i<k;i++)
    {
        for(j=0;j<ndata;j++)
            cluster_assign[i][j]=0;
    }

}
/*Select the data point as the centroid for 2nd cluster which is far from cluster 1 centroid*/
void findCentroid_2(int *data1, int **data2, int dim,int size,int current_cluster_size)
{
    int i,j,*max_dist_pt;
    double dist,mindist=0,tmp=0.0;
    for(j=0;j<size;j=j+dim)
    {
        tmp=0.0;
        for(i=0;i<dim;i++)
        {
            int val=data1[i+j];
            int val2=data2[0][i];
            tmp+=((val-val2)*(val-val2));
        }
        dist=sqrt(tmp);
        if(mindist<dist)
        {
            mindist=dist;
            int m;
            for(m=0;m<dim;m++)
            {
                data2[current_cluster_size][m]=data1[j+m];
            }
        }

    }
}

void main()
{
    int ndata=1000,i,dim=2,*data,k=5,j,size;
    int **cluster_assign,**cluster_centroid,*searchData,*clustersize,*clusterstart;
    double *cluster_raidus;
    int *dist;
    int *ptr_data,*ptr_search;
    printf("Enter the Number of Data points\n");
    scanf("%d",&ndata);
    printf("Enter the Number of Dimention\n");
    scanf("%d",&dim);
    printf("Enter the Number of Clusters\n");
    scanf("%d",&k);
    size=ndata*dim;
    //Dynamic memory allocation
    data = (int*)malloc(size * sizeof(int));
    clusterstart = (int*)malloc(k * sizeof(int));
    cluster_assign = (int**)malloc(k * sizeof(int));
    for (i=0; i<k; i++)
    {
        *(cluster_assign+i) = (int *)malloc(size * sizeof(int));
    }
    cluster_raidus =(double*)malloc(k * sizeof(double));
    searchData = (int*)malloc(dim * sizeof(int));
    cluster_centroid=(int **)malloc(k * sizeof(int *));
    for (i=0; i<k; i++)
         *(cluster_centroid+i) = (int *)malloc(dim * sizeof(int));
    clustersize = (int*)malloc(k * sizeof(int));
    for(i=0;i<k;i++)
    {
        clustersize[i]=0;
    }
    //Generating random dataset
    for(i=0;i<size;i++)
    {
        data[i]=((rand()%1000)+(rand() / (double)RAND_MAX));
    }
    //Initializing all the radius to 0
    for(i=0;i<k;i++)
    {
        cluster_raidus[i]=0;
    }
    for(i=0;i<dim;i++)
        searchData[i]=((rand()%1000)+(rand() / (double)RAND_MAX));
    BisectingKmeans(dim,size,data,k,cluster_raidus,cluster_centroid,cluster_assign,clustersize,searchData,clusterstart);
}

