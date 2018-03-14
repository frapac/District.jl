using Clustering

#Input : a laplacian and a number of clusters
#Output : returns an assignment vector that associates each point to a cluster 
function spectralCluster(laplacian::Array{Float64,2}, nbClusters::Int64)

	eigvectors = eigvecs(laplacian)
	U = eigvectors[:,1:nbClusters]  #nbClusters first eigen vectors
	
	#Clustering of the rows of U among nbClusters clusters 
	return Clustering.kmeans(U', nbClusters)
end
