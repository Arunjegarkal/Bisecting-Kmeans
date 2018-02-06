# Bisecting-Kmeans

Implementations Detail:

•	Generated Random n data points and select a random point p as centroid of cluster C1.

•	Select the farthest point q from p as centroid of Cluster C2.

•	Assign data points to clusters based on min distance.

•	Calculate midpoint of clusters as new centroids.

•	Repeat steps 3 and 4 until cluster centroid does not change.

•	Calculate SSE for Clusters.

•	Bi-portion the cluster with highest SSE.

•	Repeat step 3 to 7 until K clusters are formed.

•	Sort the data based on cluster assignment.

•	Search cluster i’s data point in which search data points falls within cluster i’s radius.
