# =============================================================================
# Functions for calculating and plotting removal metrics
# =============================================================================

#' Load results from a previous `remove_erroneous_obs()` run
#'
#' @param run_id Character. The `run_id` used in [remove_erroneous_obs()].
#'
#' @return The results list as saved by [remove_erroneous_obs()].
#'
#' @export
load_removal_results <- function(run_id) {
  if (!is.character(run_id)) stop("`run_id` must be a character string.")
  rds_file <- paste0(run_id, "_results.RDS")
  if (!file.exists(rds_file)) stop(paste("Results file not found:", rds_file))
  results <- readRDS(rds_file)
  expected <- c("par", "rem", "phi", "rmse")
  missing  <- setdiff(expected, names(results))
  if (length(missing) > 0) {
    stop(paste("Results missing components:", paste(missing, collapse = ", ")))
  }
  results
}


#' Calculate the pseudo-objective function value (pOFV) change across iterations
#'
#' The pOFV is defined as the sum of individual \eqn{-2\log L_i} (OBJ) values
#' restricted to subjects that had **no** observations removed across all
#' iterations.  The change in pOFV from one iteration to the next reflects the
#' genuine model-fit improvement attributable to each removal.
#'
#' @param results A results list from [remove_erroneous_obs()] or
#'   [load_removal_results()].
#'
#' @return A data frame with columns `ITERATION`, `pofv`, and `diff_pofv`.
#'
#' @examples
#' \dontrun{
#' results <- load_removal_results("err_rem1")
#' pofv_df <- calculate_pofv_change(results)
#' }
#'
#' @export
calculate_pofv_change <- function(results) {
  .validate_results(results)
  final_phi <- results$phi
  final_rem <- results$rem

  missing_ids <- setdiff(final_rem$ID, final_phi$ID)
  if (length(missing_ids) > 0) {
    stop(paste(
      "IDs in removed observations not found in phi data:",
      paste(missing_ids, collapse = ", ")
    ))
  }

  # Only subjects with no removals across all iterations
  clean_phi <- final_phi[!final_phi$ID %in% final_rem$ID, ]

  n <- max(clean_phi$ITERATION)
  pofv_list <- vector("list", n)
  for (i in seq_len(n)) {
    iter_phi   <- clean_phi[clean_phi$ITERATION == i, ]
    pofv_list[[i]] <- data.frame(pofv = sum(iter_phi$OBJ), ITERATION = i)
  }

  final_pofv <- dplyr::bind_rows(pofv_list) |>
    dplyr::mutate(
      diff_pofv = .data$pofv - dplyr::lag(
        .data$pofv,
        default = dplyr::first(.data$pofv),
        order_by = .data$ITERATION
      )
    )
  final_pofv
}


#' Plot metrics from the observation removal process
#'
#' Creates diagnostic plots tracking the effect of iterative observation
#' removal on model fit and parameter stability.
#'
#' @param results A results list from [remove_erroneous_obs()] or
#'   [load_removal_results()].
#' @param metric Character. Which metric to plot:
#'   \describe{
#'     \item{`"pOFV"`}{Change in pseudo-OFV at each iteration.}
#'     \item{`"thetas"`}{THETA estimates normalised to the final estimate.}
#'     \item{`"omegas"`}{OMEGA estimates normalised to the final estimate.}
#'     \item{`"sigmas"`}{SIGMA estimates normalised to the final estimate.
#'       Returns \code{NULL} with a message if all SIGMA parameters are fixed.}
#'     \item{`"rmse"`}{Root mean square error over iterations.}
#'   }
#' @param verbose Logical. Print progress messages?  Default `TRUE`.
#'
#' @return A `ggplot` object.
#'
#' @section Parameter labels:
#' When the results object contains a `param_labels` element (populated
#' automatically by [remove_erroneous_obs()]), the THETA and OMEGA plot panels
#' are labelled with the human-readable names parsed from the \code{.mod} file.
#'
#' To enable labelling, add inline comments to your \code{$THETA} and
#' \code{$OMEGA} blocks following the format \code{; <index>. <Label>}:
#'
#' \preformatted{
#' $THETA
#'   (0, 11.6) ; 1. CL
#'   (0,  9.2) ; 2. V1
#'
#' $OMEGA
#'   0.09      ; 1. IIV CL
#'   0.04      ; 2. IIV V1
#'
#' $OMEGA BLOCK(1)
#'   0.03      ; 3. IOV CL
#' $OMEGA BLOCK(1) SAME
#' }
#'
#' \code{BLOCK(1) SAME} entries are detected automatically and excluded from
#' the OMEGA plot regardless of whether a comment is present.  Range notation
#' (\code{; 3-7. IOV CL}) is also supported; only the first diagonal element
#' is plotted.
#'
#' @examples
#' \dontrun{
#' results <- load_removal_results("err_rem1")
#' plot_removal_metrics(results, metric = "pOFV")
#' plot_removal_metrics(results, metric = "thetas")
#' }
#'
#' @export
plot_removal_metrics <- function(
  results,
  metric  = c("pOFV", "thetas", "omegas", "sigmas", "rmse"),
  verbose = TRUE
) {
  .validate_results(results)
  metric <- match.arg(metric)

  if (metric == "pOFV") {
    final_pofv <- calculate_pofv_change(results)
    return(.plot_pofv(final_pofv))
  }

  n <- max(results$phi$ITERATION)

  labels <- results$param_labels

  if (metric == "thetas") {
    return(.plot_thetas(.process_parameters(results$par, "THETA", n), labels))
  }
  if (metric == "omegas") {
    return(.plot_omegas(.process_parameters(results$par, "OMEGA", n), labels))
  }
  if (metric == "sigmas") {
    params_obj <- .process_parameters(results$par, "SIGMA", n)
    if (is.null(params_obj)) {
      message("All SIGMA parameters are fixed; no SIGMA plot generated.")
      return(invisible(NULL))
    }
    return(.plot_sigmas(params_obj, labels))
  }
  if (metric == "rmse") {
    return(.plot_rmse(results$rmse))
  }
}


# =============================================================================
# Internal helpers
# =============================================================================

#' @keywords internal
.calculate_nrmse <- function(obs, pred) {
  rmse <- sqrt(mean((pred - obs)^2, na.rm = TRUE))
  rmse / mean(obs, na.rm = TRUE)
}

#' @keywords internal
.validate_results <- function(results) {
  expected <- c("par", "rem", "phi", "rmse")
  missing  <- setdiff(expected, names(results))
  if (length(missing) > 0) {
    stop(paste(
      "Results object missing components:",
      paste(missing, collapse = ", ")
    ))
  }
}

#' @keywords internal
.process_parameters <- function(final_par, param_type, n) {
  params <- final_par[, grepl(param_type, names(final_par)), drop = FALSE]

  # Exclude fixed parameters (value identical across all iterations)
  is_fixed <- vapply(params, function(x) {
    v <- x[!is.na(x)]
    length(v) == 0L || all(v == v[1L])
  }, logical(1L))
  params <- params[, !is_fixed, drop = FALSE]
  if (ncol(params) == 0L) return(NULL)

  # Store raw first and last values before normalising (for plot labels)
  raw_first <- setNames(as.numeric(params[1L, ]),  names(params))
  raw_last  <- setNames(as.numeric(params[n,  ]),  names(params))

  # Normalise each column to the value at the last iteration
  last_val <- raw_last
  last_val[last_val == 0] <- NA_real_
  norm <- as.data.frame(mapply(function(col, ref) col / ref, params, last_val))
  norm$ITERATION <- seq_len(n)

  list(data = norm, raw_first = raw_first, raw_last = raw_last)
}

#' @keywords internal
.irx_colors <- list(
  primary  = "#003B4F",
  green    = "#20794D",
  blue     = "#00769E",
  olive    = "#657422",
  red      = "#AD0000",
  gray     = "#5E5E5E",
  charcoal = "#111111",
  grid     = "#F0F0F0"
)

#' @keywords internal
.clean_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = "#F0F0F0", linewidth = 0.3
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(
        colour = "#5E5E5E", linewidth = 0.4
      ),
      strip.background = ggplot2::element_rect(fill = "#003B4F", colour = NA),
      strip.text       = ggplot2::element_text(
        colour = "white", size = 8, face = "bold"
      ),
      plot.title       = ggplot2::element_text(
        colour = "#003B4F", size = 11, face = "bold"
      ),
      plot.subtitle    = ggplot2::element_text(colour = "#5E5E5E", size = 9),
      axis.text        = ggplot2::element_text(colour = "#5E5E5E", size = 8),
      axis.title       = ggplot2::element_text(colour = "#111111", size = 9),
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 8)
    )
}

#' @keywords internal
.plot_pofv <- function(final_pofv) {
  pofv_nz <- dplyr::filter(final_pofv, .data$diff_pofv != 0)
  n   <- nrow(pofv_nz)
  n50 <- max(floor(n * 0.5), 2L)
  n25 <- max(floor(n * 0.25), 2L)

  # Three subsets: first 25%, first 50%, and 100% of removals.
  # Each gets its own facet panel with a free y-scale so the wide CI
  # on the full-data fit does not obscure the early-removal panels.
  df_fit <- rbind(
    cbind(pofv_nz[pofv_nz$ITERATION <= n25, , drop = FALSE], group = "25%"),
    cbind(pofv_nz[pofv_nz$ITERATION <= n50, , drop = FALSE], group = "50%"),
    cbind(pofv_nz,                                            group = "100%")
  )
  df_fit$group <- factor(df_fit$group, levels = c("25%", "50%", "100%"))

  ggplot2::ggplot(
    df_fit,
    ggplot2::aes(.data$ITERATION, .data$diff_pofv)
  ) +
    ggplot2::geom_point(colour = "#5E5E5E") +
    ggplot2::geom_line(alpha = 0.4, linetype = 2, colour = "#5E5E5E") +
    ggplot2::stat_smooth(
      method  = "lm",
      formula = y ~ log(x),
      colour  = "#003B4F",
      fill    = "#003B4F",
      alpha   = 0.15
    ) +
    ggplot2::facet_wrap(~group, scales = "free_y", nrow = 1) +
    .clean_theme() +
    ggplot2::ggtitle(
      "Change in pOFV per iteration",
      subtitle = paste0(
        "Restricted to subjects with no observations removed  |  ",
        "log-linear fit (y ~ log(x)) with 95\u0025 CI  |  ",
        "Facets show fit on first 25\u0025, 50\u0025, and all removals"
      )
    ) +
    ggplot2::ylab("\u0394 pOFV") +
    ggplot2::xlab("Iteration")
}

#' @keywords internal
.plot_thetas <- function(params_obj, labels = NULL) {
  if (is.null(params_obj)) return(NULL)
  theta_df  <- params_obj$data
  raw_first <- params_obj$raw_first
  raw_last  <- params_obj$raw_last
  n         <- max(theta_df$ITERATION)

  # Apply human-readable labels and remap raw_ lookup keys to match
  if (!is.null(labels) && length(labels$thetas) > 0) {
    lbl <- labels$thetas
    names(raw_first) <- ifelse(names(raw_first) %in% names(lbl),
                               lbl[names(raw_first)], names(raw_first))
    names(raw_last)  <- ifelse(names(raw_last)  %in% names(lbl),
                               lbl[names(raw_last)],  names(raw_last))
  }

  long <- theta_df |>
    tidyr::pivot_longer(-"ITERATION", names_to = "THETA") |>
    dplyr::arrange(.data$THETA)
  if (!is.null(labels) && length(labels$thetas) > 0) {
    lbl <- labels$thetas
    long$THETA <- ifelse(
      long$THETA %in% names(lbl), lbl[long$THETA], long$THETA
    )
  }

  # First and last iteration rows for endpoint labels
  first_rows       <- long[long$ITERATION == 1L, ]
  first_rows$label <- format(signif(raw_first[first_rows$THETA], 3L),
                             scientific = FALSE)
  last_rows        <- long[long$ITERATION == n, ]
  last_rows$label  <- format(signif(raw_last[last_rows$THETA], 3L),
                             scientific = FALSE)
  label_data <- if (n == 1L) last_rows else rbind(first_rows, last_rows)

  phantom_df <- data.frame(
    THETA     = rep(unique(long$THETA), each = 2L),
    value     = rep(c(0.75, 1.5), times = length(unique(long$THETA))),
    ITERATION = 1L,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(long, ggplot2::aes(.data$ITERATION, .data$value)) +
    ggplot2::geom_blank(data = phantom_df) +
    ggplot2::geom_point(colour = "#003B4F", show.legend = FALSE) +
    ggplot2::geom_line(colour = "#003B4F") +
    ggplot2::geom_hline(
      yintercept = 1, linetype = 2, colour = "#00769E", linewidth = 0.6
    ) +
    ggplot2::geom_label(
      data        = label_data,
      ggplot2::aes(label = .data$label),
      size        = 3,
      vjust       = -0.4,
      hjust       = ifelse(label_data$ITERATION == 1L, 0, 1),
      colour      = "#111111",
      fill        = "white",
      label.size  = 0,
      label.padding = ggplot2::unit(0.12, "lines")
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.15))
    ) +
    .clean_theme() +
    ggplot2::facet_wrap(~THETA, scales = "free_y") +
    ggplot2::ggtitle(
      "THETA estimates across iterations",
      subtitle = paste0(
        "Normalised to final estimate (dashed line = 1.0)  |  ",
        "Fixed parameters excluded  |  ",
        "First and last estimates annotated"
      )
    ) +
    ggplot2::ylab("Normalised parameter value")
}

#' @keywords internal
.plot_omegas <- function(params_obj, labels = NULL) {
  if (is.null(params_obj)) return(NULL)
  omega_df  <- params_obj$data
  raw_first <- params_obj$raw_first
  raw_last  <- params_obj$raw_last
  n         <- max(omega_df$ITERATION)

  # Show diagonal variance terms only
  all_cols  <- names(omega_df)
  diag_cols <- all_cols[
    grepl("^OMEGA\\.\\d+\\.\\d+$", all_cols) &
      vapply(strsplit(all_cols, "\\."), function(x) {
        length(x) == 3L && x[2L] == x[3L]
      }, logical(1L))
  ]
  if (length(diag_cols) == 0L) diag_cols <- all_cols[all_cols != "ITERATION"]

  # When labels are available, restrict to labelled columns only.
  # This drops IOV SAME entries (identical values, redundant panels).
  if (!is.null(labels) && length(labels$omegas) > 0) {
    labelled <- intersect(diag_cols, names(labels$omegas))
    if (length(labelled) > 0) diag_cols <- labelled
  }

  # Filter raw_ vectors to the surviving columns
  raw_first <- raw_first[names(raw_first) %in% diag_cols]
  raw_last  <- raw_last[names(raw_last)   %in% diag_cols]

  # Apply human-readable labels and remap raw_ lookup keys
  if (!is.null(labels) && length(labels$omegas) > 0) {
    lbl <- labels$omegas
    names(raw_first) <- ifelse(names(raw_first) %in% names(lbl),
                               lbl[names(raw_first)], names(raw_first))
    names(raw_last)  <- ifelse(names(raw_last)  %in% names(lbl),
                               lbl[names(raw_last)],  names(raw_last))
  }

  long <- omega_df[, c("ITERATION", diag_cols), drop = FALSE] |>
    tidyr::pivot_longer(-"ITERATION", names_to = "OMEGA") |>
    dplyr::arrange(.data$OMEGA)
  if (!is.null(labels) && length(labels$omegas) > 0) {
    lbl <- labels$omegas
    long$OMEGA <- ifelse(
      long$OMEGA %in% names(lbl), lbl[long$OMEGA], long$OMEGA
    )
  }

  # First and last iteration rows for endpoint labels
  first_rows       <- long[long$ITERATION == 1L, ]
  first_rows$label <- format(signif(raw_first[first_rows$OMEGA], 3L),
                             scientific = FALSE)
  last_rows        <- long[long$ITERATION == n, ]
  last_rows$label  <- format(signif(raw_last[last_rows$OMEGA], 3L),
                             scientific = FALSE)
  label_data <- if (n == 1L) last_rows else rbind(first_rows, last_rows)

  phantom_df <- data.frame(
    OMEGA     = rep(unique(long$OMEGA), each = 2L),
    value     = rep(c(0.75, 1.5), times = length(unique(long$OMEGA))),
    ITERATION = 1L,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(long, ggplot2::aes(.data$ITERATION, .data$value)) +
    ggplot2::geom_blank(data = phantom_df) +
    ggplot2::geom_point(colour = "#003B4F", show.legend = FALSE) +
    ggplot2::geom_line(colour = "#003B4F") +
    ggplot2::geom_hline(
      yintercept = 1, linetype = 2, colour = "#00769E", linewidth = 0.6
    ) +
    ggplot2::geom_label(
      data        = label_data,
      ggplot2::aes(label = .data$label),
      size        = 3,
      vjust       = -0.4,
      hjust       = ifelse(label_data$ITERATION == 1L, 0, 1),
      colour      = "#111111",
      fill        = "white",
      label.size  = 0,
      label.padding = ggplot2::unit(0.12, "lines")
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.15))
    ) +
    .clean_theme() +
    ggplot2::facet_wrap(~OMEGA, scales = "free_y") +
    ggplot2::ggtitle(
      "OMEGA variance estimates across iterations",
      subtitle = paste0(
        "Diagonal elements only, normalised to final estimate",
        " (dashed = 1.0)  |  ",
        "Fixed parameters excluded  |  First and last estimates annotated"
      )
    ) +
    ggplot2::ylab("Normalised parameter value")
}

#' @keywords internal
.plot_rmse <- function(final_rmse) {
  final_rmse$rmse <- final_rmse$rmse * 100
  ggplot2::ggplot(final_rmse, ggplot2::aes(.data$ITERATION, .data$rmse)) +
    ggplot2::geom_point(colour = "#003B4F") +
    ggplot2::geom_line(alpha = 0.5, linetype = 2, colour = "#003B4F") +
    .clean_theme() +
    ggplot2::ggtitle(
      "Normalised RMSE (nRMSE) across iterations",
      subtitle = "nRMSE = RMSE / mean(DV) x 100"
    ) +
    ggplot2::ylab("nRMSE (%)")
}

#' @keywords internal
.plot_sigmas <- function(params_obj, labels = NULL) {
  if (is.null(params_obj)) {
    message("All SIGMA parameters are fixed; skipping SIGMA plot.")
    return(NULL)
  }
  sigma_df  <- params_obj$data
  raw_first <- params_obj$raw_first
  raw_last  <- params_obj$raw_last
  n         <- max(sigma_df$ITERATION)

  # Show diagonal variance terms only
  all_cols  <- names(sigma_df)
  diag_cols <- all_cols[
    grepl("^SIGMA\\.\\d+\\.\\d+$", all_cols) &
      vapply(strsplit(all_cols, "\\."), function(x) {
        length(x) == 3L && x[2L] == x[3L]
      }, logical(1L))
  ]
  if (length(diag_cols) == 0L) diag_cols <- all_cols[all_cols != "ITERATION"]

  if (!is.null(labels) && length(labels$sigmas) > 0) {
    labelled <- intersect(diag_cols, names(labels$sigmas))
    if (length(labelled) > 0) diag_cols <- labelled
  }

  raw_first <- raw_first[names(raw_first) %in% diag_cols]
  raw_last  <- raw_last[names(raw_last)   %in% diag_cols]

  if (!is.null(labels) && length(labels$sigmas) > 0) {
    lbl <- labels$sigmas
    names(raw_first) <- ifelse(names(raw_first) %in% names(lbl),
                               lbl[names(raw_first)], names(raw_first))
    names(raw_last)  <- ifelse(names(raw_last)  %in% names(lbl),
                               lbl[names(raw_last)],  names(raw_last))
  }

  long <- sigma_df[, c("ITERATION", diag_cols), drop = FALSE] |>
    tidyr::pivot_longer(-"ITERATION", names_to = "SIGMA") |>
    dplyr::arrange(.data$SIGMA)

  if (!is.null(labels) && length(labels$sigmas) > 0) {
    lbl <- labels$sigmas
    long$SIGMA <- ifelse(
      long$SIGMA %in% names(lbl), lbl[long$SIGMA], long$SIGMA
    )
  }

  first_rows       <- long[long$ITERATION == 1L, ]
  first_rows$label <- format(signif(raw_first[first_rows$SIGMA], 3L),
                             scientific = FALSE)
  last_rows        <- long[long$ITERATION == n, ]
  last_rows$label  <- format(signif(raw_last[last_rows$SIGMA], 3L),
                             scientific = FALSE)
  label_data <- if (n == 1L) last_rows else rbind(first_rows, last_rows)

  phantom_df <- data.frame(
    SIGMA     = rep(unique(long$SIGMA), each = 2L),
    value     = rep(c(0.75, 1.5), times = length(unique(long$SIGMA))),
    ITERATION = 1L,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(long, ggplot2::aes(.data$ITERATION, .data$value)) +
    ggplot2::geom_blank(data = phantom_df) +
    ggplot2::geom_point(colour = "#003B4F", show.legend = FALSE) +
    ggplot2::geom_line(colour = "#003B4F") +
    ggplot2::geom_hline(
      yintercept = 1, linetype = 2, colour = "#00769E", linewidth = 0.6
    ) +
    ggplot2::geom_label(
      data        = label_data,
      ggplot2::aes(label = .data$label),
      size        = 3,
      vjust       = -0.4,
      hjust       = ifelse(label_data$ITERATION == 1L, 0, 1),
      colour      = "#111111",
      fill        = "white",
      label.size  = 0,
      label.padding = ggplot2::unit(0.12, "lines")
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.15))
    ) +
    .clean_theme() +
    ggplot2::facet_wrap(~SIGMA, scales = "free_y") +
    ggplot2::ggtitle(
      "SIGMA estimates across iterations",
      subtitle = paste0(
        "Diagonal elements only, normalised to final estimate",
        " (dashed = 1.0)  |  ",
        "Fixed parameters excluded  |  First and last estimates annotated"
      )
    ) +
    ggplot2::ylab("Normalised parameter value")
}

#' @keywords internal
.parse_param_name <- function(x) {
  if (is.null(x)) return(NULL)
  # Already a valid column name — return as-is
  m <- regmatches(x, regexec(
    "^(THETA|OMEGA|SIGMA)\\((\\d+)\\)$", x, perl = TRUE
  ))[[1L]]
  if (length(m) == 3L) {
    type <- m[2L]
    idx  <- m[3L]
    if (type == "THETA") return(paste0("THETA", idx))
    # SIGMA(2) -> SIGMA.2.2 , OMEGA(3) -> OMEGA.3.3
    return(paste0(type, ".", idx, ".", idx))
  }
  x
}

#' @keywords internal
.find_error_params <- function(param_labels, par_df, verbose = FALSE) {
  prop_col   <- NULL
  add_col    <- NULL
  prop_label <- NULL
  add_label  <- NULL
  source_note <- NULL

  # 1. Search THETA labels (most models encode error as THETAs)
  theta_labels <- if (!is.null(param_labels) && length(param_labels$thetas) > 0)
    param_labels$thetas else character(0L)

  if (length(theta_labels) > 0L) {
    prop_m <- names(theta_labels)[
      grepl("(?i)(prop|proportional|ruv_prop|ruv.prop)", theta_labels, perl = TRUE)
    ]
    add_m  <- names(theta_labels)[
      grepl("(?i)(\\badd\\b|additive|ruv_add|ruv.add)", theta_labels, perl = TRUE)
    ]
    if (length(prop_m) > 0L) { prop_col <- prop_m[1L]; prop_label <- theta_labels[prop_col] }
    if (length(add_m)  > 0L) { add_col  <- add_m[1L];  add_label  <- theta_labels[add_col]  }
    if (!is.null(prop_col) || !is.null(add_col)) source_note <- "THETA label"
  }

  # 2. Search SIGMA labels (classic sigma-based error models)
  sigma_labels <- if (!is.null(param_labels) && length(param_labels$sigmas) > 0)
    param_labels$sigmas else character(0L)

  if (length(sigma_labels) > 0L) {
    free_sigma <- .free_sigma_cols(par_df)
    if (is.null(prop_col)) {
      prop_s <- intersect(
        names(sigma_labels)[
          grepl("(?i)(prop|proportional|ruv_prop|ruv.prop)", sigma_labels, perl = TRUE)
        ],
        free_sigma
      )
      if (length(prop_s) > 0L) {
        prop_col <- prop_s[1L]; prop_label <- sigma_labels[prop_col]
        source_note <- if (is.null(source_note)) "SIGMA label" else source_note
      }
    }
    if (is.null(add_col)) {
      add_s <- intersect(
        names(sigma_labels)[
          grepl("(?i)(\\badd\\b|additive|ruv_add|ruv.add)", sigma_labels, perl = TRUE)
        ],
        free_sigma
      )
      if (length(add_s) > 0L) {
        add_col <- add_s[1L]; add_label <- sigma_labels[add_col]
        source_note <- if (is.null(source_note)) "SIGMA label" else source_note
      }
    }
  }

  # 3. Fall back to positional free-SIGMA columns (no labels available)
  if (is.null(prop_col) || is.null(add_col)) {
    free_sigma <- .free_sigma_cols(par_df)
    if (is.null(prop_col) && length(free_sigma) >= 1L) {
      prop_col <- free_sigma[1L]; source_note <- "positional (first free SIGMA)"
    }
    if (is.null(add_col) && length(free_sigma) >= 2L) {
      add_col  <- free_sigma[2L]
    }
  }

  if (verbose) {
    prop_str <- if (!is.null(prop_col)) {
      if (!is.null(prop_label))
        sprintf('"%s" (label: "%s")', prop_col, prop_label)
      else
        sprintf('"%s"', prop_col)
    } else {
      "not detected"
    }
    add_str <- if (!is.null(add_col)) {
      if (!is.null(add_label))
        sprintf('"%s" (label: "%s")', add_col, add_label)
      else
        sprintf('"%s"', add_col)
    } else {
      "not detected"
    }
    if (!is.null(prop_col) || !is.null(add_col)) {
      message(sprintf(
        paste0(
          "irxclean: Auto-detected error parameters (source: %s):\n",
          "  Proportional : %s\n",
          "  Additive     : %s\n",
          "  Override via prop_error_param / add_error_param",
          " (e.g. \"THETA(7)\" or \"SIGMA(1)\")."
        ),
        source_note, prop_str, add_str
      ))
    } else {
      message(paste0(
        "irxclean: Could not auto-detect error parameters.\n",
        "  No param_labels found or no prop/add keywords matched.\n",
        "  Supply prop_error_param and add_error_param manually",
        " (e.g. \"THETA(7)\" or \"SIGMA(1)\")."
      ))
    }
  }

  list(prop_col = prop_col, add_col = add_col)
}

#' @keywords internal
.free_sigma_cols <- function(par_df) {
  if (is.null(par_df)) return(character(0L))
  sigma_cols <- grep("^SIGMA", names(par_df), value = TRUE)
  if (length(sigma_cols) == 0L) return(character(0L))
  is_fixed <- vapply(par_df[, sigma_cols, drop = FALSE], function(x) {
    v <- x[!is.na(x)]
    length(v) == 0L || all(v == v[1L])
  }, logical(1L))
  sigma_cols[!is_fixed]
}

#' @keywords internal
.plot_param_trajectory <- function(results, param_col, title, subtitle = NULL) {
  par_df <- results$par
  if (!param_col %in% names(par_df)) return(NULL)

  n      <- max(par_df$ITERATION, na.rm = TRUE)
  ord    <- order(par_df$ITERATION)
  vals   <- par_df[[param_col]][ord]
  last_val <- vals[n]
  if (is.na(last_val) || last_val == 0) return(NULL)

  norm_vals <- vals / last_val
  raw_first <- vals[1L]
  raw_last  <- last_val

  # Resolve human-readable display name from labels
  display_name <- param_col
  labels <- results$param_labels
  if (!is.null(labels)) {
    all_lbl <- c(
      if (length(labels$thetas) > 0) labels$thetas else character(0L),
      if (length(labels$omegas) > 0) labels$omegas else character(0L),
      if (!is.null(labels$sigmas) && length(labels$sigmas) > 0)
        labels$sigmas else character(0L)
    )
    if (param_col %in% names(all_lbl)) display_name <- all_lbl[[param_col]]
  }

  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Normalised to final estimate (dashed = 1.0)  |  Parameter: ",
      display_name
    )
  }

  df <- data.frame(ITERATION = seq_len(n), value = norm_vals)

  first_row       <- df[df$ITERATION == 1L, , drop = FALSE]
  last_row        <- df[df$ITERATION == n,  , drop = FALSE]
  first_row$label <- format(signif(raw_first, 3L), scientific = FALSE)
  last_row$label  <- format(signif(raw_last,  3L), scientific = FALSE)
  label_data <- if (n == 1L) last_row else rbind(first_row, last_row)

  phantom_df <- data.frame(ITERATION = 1L, value = c(0.75, 1.5))

  ggplot2::ggplot(df, ggplot2::aes(.data$ITERATION, .data$value)) +
    ggplot2::geom_blank(data = phantom_df) +
    ggplot2::geom_point(colour = "#003B4F") +
    ggplot2::geom_line(colour = "#003B4F") +
    ggplot2::geom_hline(
      yintercept = 1, linetype = 2, colour = "#00769E", linewidth = 0.6
    ) +
    ggplot2::geom_label(
      data        = label_data,
      ggplot2::aes(label = .data$label),
      size        = 3,
      vjust       = -0.4,
      hjust       = ifelse(label_data$ITERATION == 1L, 0, 1),
      colour      = "#111111",
      fill        = "white",
      label.size  = 0,
      label.padding = ggplot2::unit(0.12, "lines")
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.15, 0.15))
    ) +
    .clean_theme() +
    ggplot2::ggtitle(title, subtitle = subtitle) +
    ggplot2::ylab("Normalised parameter value") +
    ggplot2::xlab("Iteration")
}
