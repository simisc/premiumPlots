#' Pipe
#'
#' Imported from \code{\link[magrittr]{pipe}}
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
NULL

#' Summarise \code{rho}
#'
#' Summarise the posterior distribution of the variable selection
#'     parameter \code{rho}. Modified version of
#'     \code{\link[PReMiuM]{summariseVarSelectRho}}, allowing plotting
#'     functions to access \code{rho} summaries.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param rho_file Path to the rho.txt file written by \code{\link[PReMiuM]{profRegr}}.
#'     If NULL (the default), path is taken from the \code{ClusObjRunInfoObj} element
#'     of \code{riskProfObj}.
#' @return A \code{\link[tibble]{tibble}} containing the mean, median,
#'     lower (0.05) and upper (0.95) quantiles for each covariate,
#'     and a ranking from the covariate with the highest value of
#'     \code{rho} to the lowest.
tabulateVarSelectRho <- function (riskProfObj, rho_file = NULL) {

    for (i in 1:length(riskProfObj))
        assign(names(riskProfObj)[i],
               riskProfObj[[i]])
    for (i in 1:length(riskProfClusObj))
        assign(names(riskProfClusObj)[i],
               riskProfClusObj[[i]])
    for (i in 1:length(clusObjRunInfoObj))
        assign(names(clusObjRunInfoObj)[i],
               clusObjRunInfoObj[[i]])

    # Temporary work-around for models fitted without variable selection
    if (!varSelect) {
        warning("Model fitted without variable selection, assuming rho=1 for all covariates.")
        return(
            tibble::tibble(
                var = covNames,
                rhoMean = 1,
                rhoMedian = 1,
                rhoLowerCI = 1,
                rhoUpperCI = 1,
                rhoRank = 1)
        )
    }

    if (is.null(rho_file)) {
        # Allows rho_file to be specified in call (e.g. if WD different to the one in which profRegr was called)
        # Propagate this option to functions that call tabulateVarSelectRho()
        rho_file <- file.path(directoryPath, paste(fileStem, "_rho.txt", sep = ""))
    }

    rhoMat <- matrix(scan(rho_file, what = double(), quiet = T),
                     ncol = nCovariates,
                     byrow = T)
    firstLine <- ifelse(reportBurnIn, nBurn / nFilter + 2, 1)
    lastLine <-
        (nSweeps + ifelse(reportBurnIn, nBurn + 1, 0)) / nFilter
    rhoMat <- rhoMat[firstLine:lastLine, ]
    rhoMean <- apply(rhoMat, 2, mean)
    rhoMedian <- apply(rhoMat, 2, median)
    rhoLowerCI <- apply(rhoMat, 2, quantile, 0.05)
    rhoUpperCI <- apply(rhoMat, 2, quantile, 0.95)

    tibble::tibble(
        var = covNames,
        rhoMean = rhoMean,
        rhoMedian = rhoMedian,
        rhoLowerCI = rhoLowerCI,
        rhoUpperCI = rhoUpperCI
    ) %>%
        dplyr::mutate(rhoRank = rank(-rhoMean))
}

#' Plot covariate profiles by cluster
#'
#' Plots the (discrete) covariate profiles, facetted vertically by cluster and
#'     (optionally) horizontally by some grouping (\code{split}) of the covariates.
#'     Covariate categories that are more common within a cluster than they are
#'     overall are shaded dark; those that are not are shaded light.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param whichCovariates A vector of indices or a vector of strings
#'     corresponding to the covariates that are to be displayed. If
#'     a single integer, this is interpreted as the 'top N' covariates,
#'     as ranked by their value of \code{rho} (requires the model to
#'     have been fitted with variable selection).
#' @param rhoMinimum As an alternative to \code{whichCovariates},
#'     plot all covariates whose value of \code{rho} is greater than
#'     some threshold (requires the model to have been fitted with
#'     variable selection).
#' @param useProfileStar To be set equal to TRUE only if a variable
#'     selection procedure has been run. The definition of the star
#'     profile is given in Liverani, S., Hastie, D. I. and Richardson,
#'     S. (2013) PReMiuM: An R package for Bayesian profile regression.
#' @param covariate_info Optional list of details about the covariates,
#'     with (some of) the following named elements:
#'     \describe{
#'         \item{\code{title}:}{String describing the covariate, default
#'         to "Covariate category".}
#'         \item{\code{levels}:}{Vector of integer values taken by the
#'         (discrete) covariates, specifying the order of the facets.}
#'         \item{\code{labels}:}{Vector of strings giving labels for each
#'         level of the covariate, in the same order as \code{levels}.}
#'         \item{\code{split}:}{A character vector specifying the covariate
#'         grouping for horizontal facetting. The first element appears
#'         in a subset of covariates to be plotted in a separate horizontal
#'         facet. The second and third elements, if present, are facet labels
#'         for covariates that do and do not meet this criterion, respectively,
#'         e.g. \code{c("_citr", "Citrullinated", "Not citrullinated")}.}
#'     }
plotProfilesByCluster <- function (riskProfObj,
                                   whichCovariates = NULL,
                                   rhoMinimum = NULL,
                                   useProfileStar = TRUE,
                                   covariate_info = list(title = "Covariate\ncategory",
                                                         levels = NULL,
                                                         labels = NULL,
                                                         split = NULL)
                                   ) {

    profileDF <- tabulateCovariateProfiles(
        riskProfObj = riskProfObj,
        whichCovariates = whichCovariates,
        rhoMinimum = rhoMinimum,
        useProfileStar = useProfileStar
    )

    if (!is.null(covariate_info$levels) & !is.null(covariate_info$labels)) {
        profileDF <- profileDF %>%
            dplyr::mutate(category = factor(category,
                                            levels = covariate_info$levels,
                                            labels = covariate_info$labels))
    }

    empirical_proportions <- profileDF %>%
        dplyr::group_by(category, covname) %>%
        dplyr::summarise(x = mean(mean)) %>%
        dplyr::group_by(covname) %>%
        dplyr::arrange(dplyr::desc(category)) %>%
        dplyr::mutate(emp_propn = cumsum(x)) %>%
        dplyr::filter(emp_propn < max(emp_propn)) %>%
        dplyr::ungroup() %>%
        dplyr::select(-x)

    covtab <- profileDF %>%
        dplyr::left_join(empirical_proportions, by = c("covname", "category")) %>%
        dplyr::group_by(cluster, category, covname,
                        fillColor, rhoMean, rhoRank, emp_propn) %>%
        dplyr::summarise(prop = mean(est)) %>%
        dplyr::ungroup()

    if (!is.null(covariate_info$split)) {
        if(length(covariate_info$split) < 3) {
            covariate_info$split <- rep(covariate_info$split, 3 - length(covariate_info$split))
            covariate_info$split <- c(covariate_info$split, sprintf("not %s", covariate_info$split[2]))
        }
        covtab <- covtab %>%
            dplyr::mutate(type = stringr::str_detect(covname, covariate_info$split[1]),
                          type = dplyr::recode(
                              as.numeric(type), `1` = covariate_info$split[2], `0` = covariate_info$split[3]
                          ))
        facetting_layer <- list(
            ggplot2::facet_grid(cluster ~ type, scales = "free_x", space = "free_x")
        )
    } else {
        facetting_layer <- list(ggplot2::facet_grid(cluster ~ .))
    }

    ggplot2::ggplot(covtab,
                    ggplot2::aes(
                        x = reorder(covname, rhoRank),
                        y = prop,
                        fill = factor(category),
                        alpha = fillColor == "high"
                    )) +
        ggplot2::geom_bar(position = "fill", stat = "identity") +
        ggplot2::geom_point(ggplot2::aes(y = emp_propn, group = category),
                            col = "black", fill = "white", alpha = 1, shape = 18) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
        )) +
        ggplot2::labs(
            x = "Covariate",
            y = "Proportion (by cluster)",
            title = "Covariate profiles",
            subtitle = "Black pips: empirical proportions in each covariate category.\nDark fill: category more prevalent within cluster than overall.") +
        ggplot2::scale_fill_discrete(name = covariate_info$title) +
        ggplot2::scale_alpha_discrete(guide = FALSE,
                                      range = c(0.25, 1)) +
        facetting_layer
}

#' Plot covariate profiles by covariate
#'
#' Plots the (discrete) covariate profiles, facetted by covariate,
#'     similar to the layout in \code{\link[PReMiuM]{plotRiskProfile}}.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param whichCovariates A vector of indices or a vector of strings
#'     corresponding to the covariates that are to be displayed. If
#'     a single integer, this is interpreted as the 'top N' covariates,
#'     as ranked by their value of \code{rho} (requires the model to
#'     have been fitted with variable selection).
#' @param rhoMinimum As an alternative to \code{whichCovariates},
#'     plot all covariates whose value of \code{rho} is greater than
#'     some threshold (requires the model to have been fitted with
#'     variable selection). Ignored if \code{whichCovariates} is specified.
#' @param useProfileStar To be set equal to TRUE only if a variable
#'     selection procedure has been run. The definition of the star
#'     profile is given in Liverani, S., Hastie, D. I. and Richardson,
#'     S. (2013) PReMiuM: An R package for Bayesian profile regression.
#' @param covariate_levels Vector of integer values taken by the
#'     (discrete) covariates, specifying the order of the facets.
#' @param covariate_labels Vector of strings giving labels for each
#'     level of the covariate, in the same order as \code{covariate_levels}.
plotCovariateProfiles <- function (riskProfObj,
                                   whichCovariates = NULL,
                                   rhoMinimum = NULL,
                                   useProfileStar = TRUE,
                                   covariate_levels = NULL,
                                   covariate_labels = NULL) {
    profileDF <- tabulateCovariateProfiles(
        riskProfObj = riskProfObj,
        whichCovariates = whichCovariates,
        rhoMinimum = rhoMinimum,
        useProfileStar = useProfileStar
    )

    length(unique(profileDF$covname)) <= 10 ||
        stop("Don't plot so many covariates at once.")

    if (!is.null(covariate_levels) & !is.null(covariate_labels)) {
        profileDF <- profileDF %>%
            dplyr::mutate(category = factor(category,
                                            levels = covariate_levels,
                                            labels = covariate_labels))
    }

    cols <- c(high = "#CC0033",
              low = "#0066CC",
              avg = "#33CC66")

    ggplot2::ggplot(profileDF, ggplot2::aes(x = factor(cluster), col = factor(fillColor))) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean), linetype = "dashed") +
        ggplot2::geom_boxplot(ggplot2::aes(y = est, fill = factor(fillColor)),
                              col = "black",
                              outlier.size = 0.5) +
        ggplot2::geom_point(ggplot2::aes(y = lower), size = 1.5, shape = 4) +
        ggplot2::geom_point(ggplot2::aes(y = upper), size = 1.5, shape = 4) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::theme(legend.position = "none",
                       strip.text = ggplot2::element_text(size = 6)) +
        ggplot2::labs(x = "Cluster", title = "Covariate profiles", y = "Probability") +
        ggplot2::facet_grid(factor(category) ~ reorder(covname, rhoRank))
}

#' @export
tabulateCovariateProfiles <- function (riskProfObj,
                                       whichCovariates = NULL,
                                       rhoMinimum = NULL,
                                       useProfileStar = TRUE,
                                       rho_file = NULL) {
    # This now does most of the work for plotCovariateProfiles
    # Use same output for new covariate profiles plot

    for (i in 1:length(riskProfObj))
        assign(names(riskProfObj)[i],
               riskProfObj[[i]])
    for (i in 1:length(riskProfClusObj))
        assign(names(riskProfClusObj)[i],
               riskProfClusObj[[i]])
    for (i in 1:length(clusObjRunInfoObj))
        assign(names(clusObjRunInfoObj)[i],
               clusObjRunInfoObj[[i]])

    xModel == "Discrete" ||
        stop("Only implemented for discrete covariates.\nUse plotRiskProfile() for Normal- or Mixed-covariate models.")

    rhotab <- tabulateVarSelectRho(riskProfObj, rho_file = rho_file)

    if (varSelect) {

        if (useProfileStar) {
            profile <- profileStar
        }

        if (is.null(whichCovariates) & !is.null(rhoMinimum)) {
            # If specify minimum rho value instead of which covariates
            rhotab <- dplyr::filter(rhotab, rhoMean >= rhoMinimum)
            whichCovariates <- rhotab$var
        }

        if (is.integer(whichCovariates) & length(whichCovariates) == 1) {
            # If whichCovariates is a single integer, treat it as the 'top N'
            rhotab <- dplyr::arrange(rhotab, rhoRank)
            whichCovariates <- rhotab$var[1:whichCovariates]
        }

    } else {
        if (useProfileStar) {
            warning("No variable selection: can't use ProfileStar--using standard profile.")
        }
        if (is.null(whichCovariates) & !is.null(rhoMinimum)) {
            stop("No variable selection: can't use rhoMinimum--specify whichCovariates.")
        }
        if (is.integer(whichCovariates) & length(whichCovariates) == 1) {
            warning("No variable selection: can't rank covariates--returning single covariate with index N.")
        }
    }


    if (!is.null(whichCovariates)) {
        if (!is.numeric(whichCovariates)) {
            whichCovariates <- match(whichCovariates, covNames)
        }
        covNames <- covNames[whichCovariates]
        nCovariates <- length(whichCovariates)
        profile <- profile[, , whichCovariates, ]
        nCategories <- nCategories[whichCovariates]
    }

    orderStat <- apply(risk, 2, median)
    sortIndex <- order(clusterSizes, orderStat, decreasing = T)
    clusterSizes <- clusterSizes[sortIndex]
    profile <- profile[, sortIndex, , ]

    profDFlist = list()

    for (j in 1:nCovariates) {
        miniDFlist = list() # still growing lists, not ideal...

        for (k in 1:nCategories[j]) {
            probMat <- profile[, , j, k]
            probMeans <- apply(probMat, 2, mean)

            probMean <- sum(probMeans * clusterSizes) / sum(clusterSizes)
            probLower <- apply(probMat, 2, quantile, 0.05)
            probUpper <- apply(probMat, 2, quantile, 0.95)

            clusterDF <- tibble::tibble(
                cluster = 1:nClusters,
                category = k - 1,
                mean = probMean,
                lower = probLower,
                upper = probUpper
            )

            profileDF <-
                tibble::tibble(cluster = rep(1:nClusters, each = nrow(probMat)),
                               est = c(probMat)) %>%
                dplyr::left_join(clusterDF, by = "cluster")

            miniDFlist[[k]] <- profileDF
        }

        profileDF <- dplyr::bind_rows(miniDFlist)
        profDFlist[[covNames[j]]] <- profileDF
    }

    dplyr::bind_rows(profDFlist, .id = "covname") %>%
        dplyr::left_join(rhotab, by = c("covname" = "var")) %>%
        dplyr::mutate(
            fillColor = ifelse(
                lower > mean,
                "high",
                ifelse(upper < mean, "low", "avg")
            ),
            fillColor = as.character(fillColor)
        )
}

#' Plot response profiles
#'
#' Plots the (categorical) response profiles, similar to the layout
#'     of the 'risk' plot in \code{\link[PReMiuM]{plotRiskProfile}}.
#' @export
#' @param riskProfObj Object of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
#' @param response_levels Vector of integer values taken by the
#'     (categorical) response, specifying the order of the facets.
#' @param response_labels Vector of strings giving labels for each
#'     level of the response, in the same order as \code{response_levels}.
plotResponse <- function (riskProfObj,
                          response_levels = NULL,
                          response_labels = NULL) {
    # Restructured like tabulateCovariateProfile()
    # 33x times faster than before
    # Corrected ordering of clusters along x-axis!

    for (i in 1:length(riskProfObj))
        assign(names(riskProfObj)[i],
               riskProfObj[[i]])
    for (i in 1:length(riskProfClusObj))
        assign(names(riskProfClusObj)[i],
               riskProfClusObj[[i]])
    for (i in 1:length(clusObjRunInfoObj))
        assign(names(clusObjRunInfoObj)[i],
               clusObjRunInfoObj[[i]])

    if (!includeResponse) {
        stop("if !includeResponse, why are you trying to plotResponse()?")
    }

    if (yModel != "Categorical") {
        stop("Only for categorical response, use plotRiskProfile() for others...")
    }

    orderStat <- apply(risk, 2, median)
    sortIndex <- order(clusterSizes, orderStat, decreasing = T)
    clusterSizes <- clusterSizes[sortIndex]
    risk <- risk[, sortIndex, ]

    miniDFlist = list() # still growing lists, not ideal...

    for (k in 1:nCategoriesY) {
        probMat <- risk[, , k] # 2D matrix for categ. k
        probMeans <- apply(probMat, 2, mean, trim = 0.005)

        probMean <- sum(probMeans * clusterSizes) / sum(clusterSizes)
        probLower <- apply(probMat, 2, quantile, 0.05)
        probUpper <- apply(probMat, 2, quantile, 0.95)

        clusterDF <- tibble::tibble(
            cluster = 1:nClusters,
            category = k - 1,
            mean = probMean,
            lower = probLower,
            upper = probUpper
        ) %>%
            dplyr::mutate(
                fillColor = ifelse(lower > mean, "high",
                                   ifelse(upper < mean, "low", "avg")),
                fillColor = as.character(fillColor)
            )

        profileDF <-
            tibble::tibble(cluster = rep(1:nClusters, each = nrow(probMat)),
                           est = c(probMat)) %>%
            dplyr::left_join(clusterDF, by = "cluster")

        miniDFlist[[k]] <- profileDF
    }

    profileDF <- dplyr::bind_rows(miniDFlist)

    if (!is.null(response_levels) & !is.null(response_labels)) {
        profileDF <- profileDF %>%
            dplyr::mutate(category = factor(category,
                                            levels = response_levels,
                                            labels = response_labels))
    }

    cols <- c(high = "#CC0033",
              low = "#0066CC",
              avg = "#33CC66")

    ggplot2::ggplot(profileDF, ggplot2::aes(x = factor(cluster), col = factor(fillColor))) +
        ggplot2::geom_hline(ggplot2::aes(yintercept = mean), linetype = "dashed") +
        ggplot2::geom_boxplot(ggplot2::aes(y = est, fill = factor(fillColor)),
                              col = "black",
                              outlier.size = 0.5) +
        ggplot2::geom_point(ggplot2::aes(y = lower), size = 1.5, shape = 4) +
        ggplot2::geom_point(ggplot2::aes(y = upper), size = 1.5, shape = 4) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "Cluster", title = "Response profile", y = "Probability") +
        ggplot2::facet_grid(. ~ factor(category))
}

#' Plot cluster sizes
#'
#' Plots the cluster sizes from one or more PReMiuM models
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotClusterSizes <- function (...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- lapply(riskprofs, function(r) {

        nClusters <- r$riskProfClusObj$nClusters
        clusterSizes <- r$riskProfClusObj$clusterSizes

        orderStat <- apply(r$risk, 2, median)
        sortIndex <- order(clusterSizes, orderStat, decreasing = T)
        clusterSizes <- clusterSizes[sortIndex]

        tibble::tibble(cluster = 1:nClusters,
                       clusterSize = clusterSizes)

    })
    data <- dplyr::bind_rows(data, .id = "model")

    ggplot2::ggplot(data, ggplot2::aes(
        x = factor(cluster),
        y = clusterSize,
        col = model
    )) +
        ggplot2::geom_point(size = 3,
                            position = ggplot2::position_dodge(width = .5)) +
        ggplot2::labs(title = "Cluster size",
                      y = "Number of subjects",
                      x = "Cluster") +
        ggplot2::expand_limits(y = 0) +
        scale_colour_discrete(guide = length(riskprofs) > 1)
}

#' Plot similarity matrix
#'
#' Plots the similarity matrix, in the original order and sorted by
#'     (hierarchical) cluster to illustrate cluster structure.
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotSimilarityMatrix <- function(...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- lapply(riskprofs, function(d) {

        disSimMat <- d$riskProfClusObj$clusObjDisSimMat

        n <- d$riskProfClusObj$clusObjRunInfoObj$nSubjects

        dsmat <- matrix(0, nrow = n, ncol = n)
        dsmat[lower.tri(dsmat)] <- disSimMat
        dsmat[upper.tri(dsmat)] <- t(dsmat)[upper.tri(dsmat)]

        dsmat <- 1 - vec2mat(disSimMat, nrow = n)

        ddist <- as.dist(1 - dsmat)
        hc <- hclust(ddist)
        dsmat2 <- dsmat[hc$order, hc$order]

        dsmat[lower.tri(dsmat)] <- NA
        dsmat2[lower.tri(dsmat2)] <- NA

        m <- tibble::as_tibble(reshape2::melt(dsmat, na.rm = TRUE))
        m2 <- tibble::as_tibble(reshape2::melt(dsmat2, na.rm = TRUE))
        m$type <- "Subjects in original order"
        m2$type <- "Subjects reordered by hclust"
        m3 <- dplyr::bind_rows(m, m2)
    })
    data <- dplyr::bind_rows(data, .id = "model")

    ggplot2::ggplot(data, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
        ggplot2::geom_tile() +
        ggplot2::facet_grid(model ~ type) +
        ggplot2::scale_y_reverse() +
        ggplot2::theme(
            line = ggplot2::element_blank(),
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
        ) +
        ggplot2::scale_fill_gradient(
            name = "",
            low = "white",
            high = "black",
            limits = c(0, 1)
        ) +
        ggplot2::coord_fixed() +
        ggplot2::labs(title = "Similarity matrix")
}

#' Plot cumulative distribution of \code{rho}
#'
#' Plots the cumulative distribution of the variable selection
#'     parameter \code{rho}. Can be interpreted as the proportion of
#'     covariates that would be dropped if only those above a given
#'     threshold of \code{rho} were included.
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotVarSelectRho <- function(...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- lapply(riskprofs, function(m) {
        rho <- tabulateVarSelectRho(m)$rhoMean # Add option for rho_file location
        ecd <- ecdf(rho)(rho) # ecdf(rho) returns a function!
        tibble::tibble(rho = rho, ecd = ecd)
    })
    data <- dplyr::bind_rows(data, .id = "model")
    ggplot2::ggplot(data, ggplot2::aes(x = rho, y = ecd, col = model)) +
        ggplot2::geom_line() +
        ggplot2::labs(title = "ECDF of rho",
                      x = "rho (mean for each covariate)",
                      y = "cumulative proportion of covariates")
}

#' Plot distribution of \code{rho} per covariate
#'
#' Plots the mean, median and CI of the variable selection
#'     parameter \code{rho} for each covariate.
#' @export
#' @param ... Object(s) of type \code{riskProfObj}, output of
#'     \code{\link[PReMiuM]{calcAvgRiskAndProfile}}.
plotRhoDistributions <- function(...) {

    riskprofs <- list(...)
    if (is.null(names(riskprofs))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(riskprofs) <- names
    }

    data <- lapply(riskprofs, function(m) {
        tabulateVarSelectRho(m) # Add option for rho_file location
    })
    data <- dplyr::bind_rows(data, .id = "model") %>%
        dplyr::rename(mean = rhoMean, median = rhoMedian) %>%
        tidyr::gather(centre, value, mean, median)
    ggplot2::ggplot(data, ggplot2::aes(x = factor(var), col = model, group = model)) +
        ggplot2::geom_linerange(ggplot2::aes(ymin = rhoLowerCI, ymax = rhoUpperCI), position = ggplot2::position_dodge(width = 0.25)) +
        ggplot2::geom_point(ggplot2::aes(y = value, shape = centre), position = ggplot2::position_dodge(width = 0.25)) +
        ggplot2::labs(x = "Covariate", y = "rho") +
        ggplot2::scale_shape_manual(name = NULL,
                                    values = c(mean = 16, median = 4))
}

#' Make a coda object PReMiuM samples
#'
#' Assemble MCMC samples from \code{\link[PReMiuM]{profRegr}} into
#'     a \code{\link[coda]{mcmc.list}} object with multiple chains.
#'     Useful for convergence diagnositics.
#' @export
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @return A \code{\link[coda]{mcmc.list}}.
codaFromPremium <- function(global.parameter, ...) {

    models <- list(...)
    if (is.null(names(models))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(models) <- names
    }

    data <- lapply(models, function(m) {
        for (i in 1:length(m)) {
            assign(names(m)[i], m[[i]])
        }
        parFileName <- file.path(directoryPath,
                                 paste(fileStem, "_", global.parameter, ".txt", sep = ""))
        parData <- coda::as.mcmc(read.table(parFileName)$V1)
    })
    coda::mcmc.list(data)
}

#' Get PReMiuM hyperparameters
#'
#' Get hyperparameters used to fit a PReMiuM model.
#' @export
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @return A matrix with a named column for each model and a named
#'     row for each hyperparameter.
getHyperparams <- function(...) {

    # Must be a cleaner way to do this... hyperparameters are only
    # written to file if they differ from the defaults...

    models <- list(...)
    if (is.null(names(models))) {
        # Name models if ... arguments are not supplied with names
        dots <- substitute(list(...))[-1]
        names <- sapply(dots, deparse)
        names(models) <- names
    }

    data <- sapply(models, function(m) {
        for (i in 1:length(m)) {
            assign(names(m)[i], m[[i]])
        }
        parFileName <- file.path(directoryPath,
                                 paste(fileStem, "_hyper.txt", sep = ""))
        hypstr <- read.table(parFileName)$V1
        hypstr <- stringr::str_split(hypstr, "=")
        hyptab <- t(sapply(hypstr, function (h) {
            h
        }))
        hyp <- as.numeric(hyptab[, 2])
        names(hyp) <- hyptab[, 1]
        hyp
    })
    data
}

#' Render PReMiuM report
#'
#' Render a standardised HTML report on PReMiuM model(s) using
#'     markdown, with covariate/response profiles and MCMC diagnostics.
#' @export
#' @param ... Object(s) of type \code{runInfoObj}, output of
#'     \code{\link[PReMiuM]{profRegr}}.
#' @param filename File name (.html) to save markdown report.
#' @param rmd.template An alternate Rmd template to populate using
#'     the models in \code{...}.
#' @return An HTML file written to \code{filename}.
renderPremiumReport <-
    function(...,
             filename = "premium_report.html",
             rmd.template = NULL) {

        premium.models <- list(...)
        if (is.null(names(premium.models))) {
            # Name models if ... arguments are not supplied with names
            dots <- substitute(list(...))[-1]
            names <- sapply(dots, deparse)
            names(premium.models) <- names
        }

        length(premium.models) >= 1 ||
            stop("Supply at least one premium model.")
        if (is.null(rmd.template)) {
            rmd.template <-
                system.file("rmd/premium_report.Rmd", package = "premiumPlots")
        }
        knitr.root.directory <- getwd()
        filepath <- file.path(knitr.root.directory, filename)
        rmarkdown::render(rmd.template, output_file = filepath)
    }
