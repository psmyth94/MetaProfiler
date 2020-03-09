sourceCpp("qtclust.cpp")

qtclust <- function(data, cluster_by, radius, distance_method, group_by = NULL, centers = NULL, 
                    start = NULL, end = NULL, element_wise = F, progress = F) {
  if(is.null(centers))
  {
    centers = cluster_by
  }
  cols_needed <- unique(c(cluster_by, centers, group_by))
  gr <- copy(data[,..cols_needed])
  gr$idx <- 1:nrow(data)
  if(!is.null(group_by)) {
    gr <- gr[, id := .GRP, by = group_by]
    gr <- gr[order(id)]
  } else {
    gr$id <- 1
  }
  
  if(element_wise)
  {
    if(length(distance_method) != length(cluster_by) & length(distance_method) != 1) {
      stop(paste("length of distance_method must be either 1 or", length(cluster_by), "for clustering using element-wise distances"))
    }
    if(length(distance_method) == 1)
    {
      distance_method = rep(distance_method, length(cluster_by))
    }
    if (is.null(start) & is.null(end))
    {
      start = 1
      end = length(cluster_by)
    }
  } else {
    if (length(start) != length(distance_method) & length(end) != length(distance_method)) {
      stop("length of start and end must be the same as distance_method for clustering using row-wise distances")
    }
    if (is.null(start) & is.null(end))
    {
      start = 1
      end = length(cluster_by)
    }
  }
  clust <- qtclust_c (as.matrix(gr[,..cluster_by]), n_groups = max(gr$id), id = gr$id,
                      groups = as.matrix(gr[,..centers]),
                      radius = radius, method = distance_method, start = start - 1,
                      end = end - 1, element_wise = element_wise, verbose = progress)
  clust$centers <- setnames(as.data.table(clust$centers), c(centers, "N"))
  gr$cluster <- clust$cluster 
  clust$centers[,group_by] <- unique(gr[,c(group_by, "cluster"), with = F])[,-"cluster"]
  clust$SD <- setnames(as.data.table(clust$SD), centers)
  clust
}
