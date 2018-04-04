################################################################################
# District.jl
################################################################################
# Implement spectral clustering method.
################################################################################

export spectralclustering

using Clustering

#Input : a laplacian and a number of clusters
#Output : returns an assignment vector that associates each point to a cluster 
function spectralclustering(laplacian::Array{Float64,2}, nbclusters::Int64)

    eigvectors = eigvecs(laplacian)
    U = eigvectors[:,1:nbclusters]  #nbclusters first eigen vectors
    
    #Clustering of the rows of U among nbclusters clusters 
    return Clustering.kmeans(U', nbclusters)
end
