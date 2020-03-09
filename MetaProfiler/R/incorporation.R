cluster_rates <- function(Object,
                          var,
                          n,
                          timepoints = Object@timepoints,
                          cluster_method = "ward.D",
                          distance_method = "euclidean",
                          cluster_names = LETTERS[1:n],
                          sort = T)
{
  if(length(cluster_names) < n) {
    stop(paste("length of `cluster_names` must equal", n))
  }
  var_cols <- get_cols(var, timepoints)
  mat <- impute(Object, vars = var, timepoints = timepoints)@master_tbl
  mat_var1 <- mat[,..var_cols]
  if(distance_method == "correlation") {
    mat_dist <- as.dist(1 - cor(t(as.matrix(mat_var1))))
  } else {
    mat_dist <- dist(mat_var1, method = distance_method)
  }
  hc <- hclust(mat_dist, method = cluster_method)
  k <- as.factor(cluster_names[cutree(hc, n)])
  mat_var2 <- mat_var1[, lapply(.SD, median), by = k]
  if(sort)
  {
    ord <- mat_var2[order(rowSums(mat_var2[, mapply(function(x, y) x - median(y), .SD[,-"k"], mat_var1)])), k]
    levels(k) <- setNames(lapply(ord, as.character), cluster_names[1:n])
  }
  k
}


heatmap <- function(Object,
                    var,
                    timepoints = Object@timepoints,
                    cluster_names = NULL,
                    split = NULL,
                    cluster_method = "ward.D",
                    distance_method = "euclidean",
                    filename = NULL,
                    width = NA,
                    height = NA,
                    sort = T,
                    plot = T,
                    ...)
{
  if(is.null(split) & !is.null(cluster_names)) {
    split <- length(cluster_names)
  }
  mat <- impute(Object, var = var, timepoints)@master_tbl[,get_cols(var, timepoints),with=F]
  mat <- as.matrix(mat)
  colnames(mat) <- as.character(timepoints)
  q <- seq(0.25,0.75,length.out = 5)
  val <- sapply(q, quantile, x = mat)
  at <- c(seq(0, val[1], length.out = 4)[-4], val, seq(val[5], 100, length.out = 4)[-1])
  col_fun = colorRamp2(at, rev(brewer.pal(11, "Spectral")))
  require(ComplexHeatmap)
  dhm <- densityHeatmap(
    column_title = "",
    col = brewer.pal(9, "Blues"),
    mat, ylab = paste0(var, " (%)"),
    ylim = c(0,100),
    height = unit(0.93, "npc"),
    heatmap_legend_param = list(
      title_position = "topcenter",
      grid_width = unit(0.3, "in"),
      legend_direction = "horizontal",
      title = "Density"
    )
  )
  nam <- if(var == "Heavy RIA") {
    "Relative Isotopic Abundance"
  } else if (var == "RIF") {
    "Relative Isotopic Fraction"
  } else {
    "Labeling Ratio"
  }
  hm <- Heatmap(
    mat,
    name = nam,
    col = col_fun,
    row_split = split,
    cluster_row_slices = T,
    cluster_rows = T,
    cluster_columns = F,
    show_row_names = F,
    row_title_rot = 0,
    clustering_distance_rows = distance_method,
    clustering_method_rows = cluster_method,
    na_col = "black",
    heatmap_legend_param = list(
      title_position = "leftcenter-rot",
      grid_width = unit(2, "in"),
      col_fun = col_fun,
      legend_direction = "horizontal",
      side = "left"
    ),
    height = unit(5, "cm"),
    column_names_rot = 0
  )
  ord <- row_order(hm)
  panel_fun = function(index, nm) {
    xrange = c(Object@time_zero, max(timepoints))
    xaxis = annotation_axis_grob(at = round(seq(xrange[1], xrange[2], length.out = 5)[-1]),
                                 side = "bottom", facing = "outside",
                                 gp = gpar(cex = 0.6666667))
    yaxis = annotation_axis_grob(at = seq(0, 100, 25),
                                 side = "left", facing = "outside",
                                 gp = gpar(cex = 0.6666667))
    pushViewport(viewport(xscale = xrange, yscale = c(0, 100)))
    grid.rect(gp = gpar(col = "black", fill = NA))
    grid.draw(xaxis)
    grid.draw(yaxis)
    cols <- get_cols(var, timepoints)
    df <- Object@master_tbl[index, ..cols]
    x = seq(xrange[1], xrange[2], length.out = 100)
    cols <- paste0(rep(var, 100), " (", Object@time_unit, " ", x, ")")
    df2 <- as.data.table(Object@model(subset = index, var = var, x = x))
    colnames(df2) <- cols
    df[,paste0(var, " (", Object@time_unit, " ", 0, ")")] <- 0
    suppressWarnings(df <- melt(df))
    df <- df[, lapply(.SD, median, na.rm = T), by = "variable"]
    df$variable <- as.numeric(stringi::stri_extract_all_regex(df$variable, "\\d+"))
    df2[,paste0(var, " (", Object@time_unit, " ", 0, ")")] <- 0
    suppressWarnings(df2 <- melt(df2))
    df2 <- df2[, lapply(.SD, median, na.rm = T), by = "variable"]
    df2$variable <- as.numeric(stringi::stri_extract_all_regex(df2$variable, "[\\d\\.]+"))
    g1 <- ggplot(df, aes(x = variable, y = value)) +
      coord_cartesian(xlim = xrange, ylim = c(0,100)) +
      theme_void() + geom_point() + geom_line() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme(legend.position = "none", axis.title = element_blank())
    grid.draw(ggplotGrob(g1))
    popViewport()
  }
  anno = anno_zoom(align_to = ord, which = "row", panel_fun = panel_fun,
                   link_width = unit(6, "mm"),
                   size = unit(2, "cm"), width = unit(4, "cm"), gap = unit(5, "mm"),
                   link_gp = gpar(col = "grey", fill = NA), extend = unit(0, "cm"))
  
  hm <- Heatmap(
    mat,
    right_annotation = rowAnnotation(foo = anno),
    name = nam,
    col = col_fun,
    row_split = split,
    cluster_row_slices = T,
    cluster_rows = T,
    cluster_columns = F,
    show_row_names = F,
    row_title_rot = 0,
    clustering_distance_rows = distance_method,
    clustering_method_rows = cluster_method,
    na_col = "black",
    heatmap_legend_param = list(
      legend_direction = "horizontal",
      title_position = "topcenter",
      grid_width = unit(0.1, "in"),
      col_fun = col_fun
    ),
    height = unit(1, "npc"),
    column_names_rot = 0,
    column_names_centered = T
  )
  if(!is.null(cluster_names)){
    rnk <- 1:split
    if(sort) {
      rnk <- rank(rowSums(sapply(as.data.frame(mat), function(x) sapply(ord, function(y) median(x[y]) - median(x)))))
    }
    hm@row_title = cluster_names[rnk]
  }
  lst <- dhm %v% hm
  if(plot) {
    draw(lst, ...)
  }
  if(!is.null(filename)) {
    save_plot(filename = filename, plot = lst, width = width, height = height, column_gap = unit(0, "cm"), ...)
  }
  lst
}

get_enrichment_from_cluster <- function(Object,
                                        group,
                                        var,
                                        n,
                                        timepoints = Object@timepoints,
                                        cluster_names = LETTERS[1:n],
                                        cluster_method = "ward.D",
                                        distance_method = "euclidean",
                                        conf_pro_only = F,
                                        threshold = 0.05, 
                                        strategy = c("BH", "rFDR", "cFDR", "p-value"),
                                        sort = T) {
  
  strategy <- match.arg(strategy)
  data <- impute(Object, vars = var, timepoints)@master_tbl
  data$Cluster <- cluster_rates(Object,
                                n = n,
                                timepoints = timepoints,
                                cluster_method = cluster_method,
                                distance_method = distance_method,
                                var = var,
                                cluster_names = cluster_names, sort = sort)
  
  if(conf_pro_only) {
    if(grepl("GO", group)) {
      cols <- colnames(data)[!(colnames(data) %in% group)]
      f <- function(x) {
        stringi::stri_trim_both(unlist(stringi::stri_split_fixed(x, ";")))
      }
      data <- data[, lapply(.SD, f), by = cols]
    }
    keep <- data[, sum(unique), by = group]
    data <- data[get(group) %in% keep[[group]][keep$V1 > 1],]
    if(group == "NOG category")
    {
      data <- data[!grepl("Function unknown", get(group)),]
      data <- data[!grepl("predict", get(group)),]
      # data$`COG category` <- name2letter[data$`COG category`, letter, on = "name"]
    }
    if(group == "COG category")
    {
      data <- data[!grepl("Function unknown", get(group)),]
      data <- data[!grepl("predict", get(group)),]
      # data$`COG category` <- name2letter[data$`COG category`, letter, on = "name"]
    }
    mat <- dcast(data, get(group) ~ Cluster, fun = function(x) length(unique(x)), value.var = "protein_razor")
  } else {
    data2 <- data[, .(group = unlist(.SD), Cluster = Cluster), by = 1:nrow(data), .SDcols = group]
    mat <- dcast(data2, group ~ Cluster, fun = length, value.var = "group")
  }
  mat <- mat[!is.na(unlist(mat[,1])),]
  cols <- colnames(mat)[-1]
  mat$r_total <- rowSums(mat[,-"group"])
  mat$all_total <- nrow(data)
  data3 <- data[,.(N = .N), by = Cluster]
  mat <- cbind(mat, as.data.table(setNames(as.list(data3$N), paste0(data3$Cluster, "_tot"))))
  p_mat <- sapply(seq_along(cols), function(i) {
    all_total <- mat$all_total
    r_total <- mat$r_total
    c_total <- mat[[paste0(cols[i], "_tot")]]
    count <- mat[[cols[i]]]
    phyper(count - 1, r_total, all_total - r_total, c_total, lower.tail = F)
  })
  p <- sort(p_mat)
  if(strategy == "rFDR") {
    sig <- sapply(rank(p_mat, ties.method = "first"), function(i) {
      idx <- round(i * 1.6)
      if(idx > length(p)) {
        return(1)
      }
      fdr <- (p[idx]*length(p))/sum(p <= p[idx])
      if(fdr > 1) {
        return(1)
      }
      fdr
    })
  } else if(strategy == "cFDR") {
    sig <- sapply(rank(p_mat, ties.method = "first"), function(i) {
      c <- 0
      for(n in 1:i) {
        c <- c + 1/(i-n+1)
      }
      fdr <- c*((p[i]*length(p))/sum(p <= p[i]))
      if(fdr > 1) {
        return(1)
      }
      fdr
    })
  } else if (strategy == "BH") {
    sig <- p.adjust(p_mat, method = "BH")
    
  } else {
    sig <- p_mat
  }
  sig_mat <- matrix(sig, ncol = ncol(p_mat))
  keep <- rowSums(sig_mat < threshold, na.rm = T) > 0
  if(all(!keep)) {
    warning(paste0("No enrichment in \"", group, "\"."))
    return(NULL)
  }
  mat <- mat[keep,]
  sig_mat <- sig_mat[keep,]
  rownames(sig_mat) <- mat[,group]
  colnames(sig_mat) <- cols
  sig_mat[sig_mat >= threshold] <- NA_real_
  sig_mat
}

cool_plot <- function(Object, vars, by, on = NULL, name = NULL, add_points = T,
                      top_n_by = NULL, min_n_by = 3, top_n_on = NULL, min_n_on = NULL,
                      nrow = NULL, ncol = NULL,
                      filename = NULL, plot = T, width = NA, height = NA) {
  df <- Object@master_tbl
  grp <- unique(c(on, by))
  # df[[by]] <- gsub("\\W", "_", df[[by]])
  df2 <- df[,..grp]
  time_range = range(c(Object@time_zero, Object@timepoints))
  x = seq(time_range[1],time_range[2],length.out = 100)
  if(add_points) {
    x = sort(unique(c(x, Object@timepoints)))
  } else {
    mid_point = diff(time_range)/2
    x = sort(unique(c(x, mid_point)))
  }
  value <- lapply(vars, Object@model, x = x)
  value <- do.call(cbind, value)
  cols2 <- get_cols(vars, x)
  colnames(value) <- cols2
  df2 <- cbind(df2, value)
  cols1 <- get_cols(vars)
  df <- na.omit(df, grp)
  df2 <- na.omit(df2, grp)
  if(is.null(name)){
    grp2 = grp
    name <- unique(df[,..grp2])
  } else if (!inherits(name, "data.table")) {
    name2 = 
    grp2 = colnames(df2)[df2 %like% name2]
  }
  df <- df[name, , on = grp2]
  df2 <- df2[name, , on = grp2]
  if(!is.null(on)) {
    keep_on <- df[, .N, by = on]
    keep_on <- keep_on[order(N)]
    top_n_on <- ifelse(is.null(top_n_on), ifelse(is.null(min_n_on), nrow(keep_on), sum(keep_on$N >= min_n_on)), top_n_on)
    keep_on <- keep_on[,tail(.SD, top_n_on)]
    df <- df[keep_on, on = on]
    df2 <- df2[keep_on, on = on]
  }
  if(by %like% "category") {
    df <- df[get(by) != "Function unknown"]
    df2 <- df2[get(by) != "Function unknown"]
  }
  if(!is.null(on) && on %like% "category") {
    df <- df[get(on) != "Function unknown"]
    df2 <- df2[get(on) != "Function unknown"]
  }
  keep_by <- df[, .N, by = by]
  keep_by <- keep_by[order(N)]
  top_n_by <- ifelse(is.null(top_n_by), ifelse(is.null(min_n_by), nrow(keep_by), sum(keep_by$N >= min_n_by)), top_n_by)
  keep_by <- keep_by[,tail(.SD, top_n_by)]
  df <- df[keep_by, on = by]
  df2 <- df2[keep_by, on = by]
  
  df <- df[, lapply(.SD, median, na.rm = T), by = grp, .SDcols = cols1]
  df2 <- df2[, lapply(.SD, mean, na.rm = T), by = grp, .SDcols = cols2]
  df <- melt(df, id.vars = grp)
  df$x <- as.numeric(stringi::stri_extract_all_regex(df$variable, "\\d+"))
  df$variable <- stringi::stri_extract_first_regex(df$variable, paste0(vars, collapse = "|"))
  df2 <- melt(df2, id.vars = grp, value.name = "value1")
  df2$x <- as.numeric(stringi::stri_extract_all_regex(df2$variable, "[\\d\\.]+"))
  df2$variable <- stringi::stri_extract_first_regex(df2$variable, paste0(vars, collapse = "|"))
  cols3 <- setdiff(colnames(df2), "value1")
  if(!add_points) {
    df <- setnames(df2[,c(cols3, "value1"),with = F], c(cols3, "value"))
    df <- df[x == mid_point] 
  }
  df2 <- df[df2,, on = cols3]
  if(by %like% "category") {
    df2 = df2[df2[[by]] != "Function unknown"]
    Letter_NOG_Category <- fread("Letter_NOG_Category.csv")
    df2[[by]] <- Letter_NOG_Category[df2[[by]], Letter, on = "Name"]
  }
  df2$variable2 = df2$variable
  if(!is.null(on)) {
    df2$variable2 = df2[[on]]
  }
  df2$variable = paste0(df2$variable, " (%)")
  
  col = NULL
  n <- length(unique(df2$variable2))
  if(!is.null(on)) {
    col = c('#222222', '#F3C300', '#875692', '#F38400', '#A1CAF1', '#BE0032', '#C2B280', '#848482', '#008856', '#E68FAC', '#0067A5', '#F99379', '#604E97', '#F6A600', '#B3446C', '#DCD300', '#882D17', '#8DB600', '#654522', '#E25822', '#2B3D26')
    col = col[1:n]
    df2[["variable2"]] <- factor(df2[["variable2"]], levels = names(sort(table(df2[["variable2"]]), decreasing = F)))
  }
  g <- ggplot(df2, aes(x = x, y = value1, group = variable2)) +
    theme_classic() +
    # labs(title = stri_replace_first_fixed(name, " ", "\n")) +
    geom_line(eval(parse(text = ifelse(
      is.null(col),
      "aes(linetype = variable)",
      "aes(color = variable2)"
    )))) + 
    xlab("Time (day)") +
    geom_point(eval(parse(text = ifelse(
      is.null(col),
      "aes(x = x, y = value, shape = variable)",
      "aes(x = x, y = value, shape = variable2, color = variable2)"
    )))) +
    scale_shape_manual(values=1:n) +
    scale_color_manual(values = col) +
    coord_cartesian(ylim = c(0,100)) +
    guides(color = guide_legend(title = on),
           linetype = guide_legend(title = "y-axis"),
           shape = guide_legend(title = ifelse(is.null(on), "y-axis", on))) +
    facet_wrap(vars(get(by)), nrow = nrow, ncol = ncol, scales = "free",
               strip.position = "top") +
    theme(strip.background = element_blank(),
                     strip.text = element_text(size = 12),
                     text = element_text(family = "sans"),
                     strip.placement = "outside",
                     axis.title.y = element_blank(),
                     plot.title = element_text(hjust = 0.3, size = 12),
                     legend.position = ifelse(T, "right", "none"))
  if(!is.null(filename)) {
    val = ifelse(is.null(col), max(strwidth(df2$variable, units = "inches", cex = fontsize(10))), max(strwidth(df2$variable2, units = "inches", cex = fontsize(10))))
    if(is.na(height) & is.na(width)) {
      n_panels <- length(unique(ggplot_build(g)$data[[1]]$PANEL))
      dims = wrap_dims(n_panels)
      width = 3 * dims[2]
      height = 2.25 * dims[1]
    }
    save_plot(filename = filename, plot = g, width = width + val, height = height)
  }
  if(plot) {
    grid.draw(g)
  }
  g
}


