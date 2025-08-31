# =========================================================
# OccuPast (OCP) Toolkit — Harmonization, Aoristics & Occupancy
# =========================================================

# ---- deps ----
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(digest)
library(withr)   
library(mgcv)    

# ===============================
# Utilities & Terminology
# ===============================


#' ocp_save_plot — Save a ggplot with good defaults
#'
#' Saves a ggplot object to disk, creating directories as needed.
#' If the file ends with `.png` and `ragg` is available, uses the
#' high-quality `ragg::agg_png` device for sharper text/lines.
#' Otherwise falls back to `ggplot2::ggsave`. Width/height are inches.
ocp_save_plot <- function(p, path, width = 8, height = 4.5, dpi = 300) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (requireNamespace("ragg", quietly = TRUE) && grepl("\\.png$", path, ignore.case = TRUE)) {
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi, device = ragg::agg_png)
  } else {
    ggplot2::ggsave(path, p, width = width, height = height, dpi = dpi)
  }
}

#' ocp_write_csv — Write a CSV with UTF-8 encoding
#'
#' Writes a data frame to CSV (via `readr::write_excel_csv`) and
#' creates parent directories if required. The "excel" variant
#' ensures locale issues (e.g. special characters) are handled sensibly.
ocp_write_csv <- function(df, path, .id = NULL) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  
  if (is.null(df)) {
    stop("ocp_write_csv: object is NULL for path: ", path, call. = FALSE)
  }
  
  # If it's a list, try to bind if all elements are data.frames
  if (inherits(df, "list") && !is.data.frame(df)) {
    all_df <- vapply(df, is.data.frame, logical(1))
    if (length(df) > 0 && all(all_df)) {
      df <- dplyr::bind_rows(df, .id = .id %||% "source")
    } else {
      bad_classes <- unique(vapply(df, function(x) paste(class(x), collapse = "/"), character(1)))
      stop(
        "ocp_write_csv: got a non-data.frame list for path: ", path, "\n",
        "  classes in list: ", paste(bad_classes, collapse = ", "), "\n",
        "  → Pick a specific component or bind into a data.frame before writing.",
        call. = FALSE
      )
    }
  }
  
  # Coerce common almost-DFs
  if (is.matrix(df)) df <- as.data.frame(df)
  if (is.atomic(df)) df <- data.frame(value = df, stringsAsFactors = FALSE)
  
  if (!is.data.frame(df)) {
    stop("ocp_write_csv: Expected a data.frame; got: ", paste(class(df), collapse = "/"),
         " for path: ", path, call. = FALSE)
  }
  
  readr::write_excel_csv(df, path)
}


#' oc_set/get/_seed & _derive_seed - ensure uniform seed handling
#'
#' When used consistently, makes sure that only one seed or local derivations
#' of that seed are used across a project. Should be set at top of file.
ocp_set_global_seed <- function(seed) {
  options(ocp.seed = as.integer(seed)); invisible(seed)
}
ocp_get_global_seed <- function(default = 2025L) {
  as.integer(getOption("ocp.seed", default))
}
ocp_derive_seed <- function(..., base = ocp_get_global_seed()) {
  # deterministic sub-seeds from a base seed + tags
  digest::digest2int(paste(c(base, ...), collapse = "::"))
}

#' check_required — Assert presence of required columns
#'
#' Throws a descriptive error if any of `required` columns are missing
#' from `data`. Helps fail-fast in early pipeline stages.
check_required <- function(data, required, where="data frame") {
  miss <- setdiff(required, names(data))
  if (length(miss)) stop(sprintf("Missing required columns in %s: %s", where, paste(miss, collapse=", ")), call.=FALSE)
}

# Traffic-light helper (low is good): green ≤ g, red ≥ r, amber in-between
ocp_status <- function(x, green_thr, red_thr) {
  dplyr::case_when(
    !is.finite(x)            ~ "grey",
    x <= green_thr           ~ "green",
    x >= red_thr             ~ "red",
    TRUE                     ~ "orange"
  )
}

#' fmt_manifest — Format a named list as key=value bullets
#'
#' Converts a named list into a compact string like
#' "x=1 · y=a,b · z=TRUE". Useful for plot subtitles/manifests.
fmt_manifest <- function(lst) {
  kv <- purrr::imap_chr(lst, ~ sprintf("%s=%s", .y, paste0(.x, collapse=",")))
  paste(kv, collapse=" · ")
}

#' bucket_coverage_fraction — Fraction of bucket covered by [gs, ge)
#'
#' Computes, for each bucket b (width = w, offset = o), the fraction
#' \eqn{c(b) = | [b, b+w) ∩ [g_s, g_e) | / w}, where \eqn{g_s} and \eqn{g_e}
#' are the global min/max horizons across the harmonization table.
#' Output: tibble(horizon_bucket, cover_frac).
bucket_coverage_fraction <- function(harm_tab, width, offset = 0, eps = 1e-9) {
  stopifnot(all(c("horizon_start","horizon_end") %in% names(harm_tab)))
  gs <- min(harm_tab$horizon_start, na.rm = TRUE)
  ge <- max(harm_tab$horizon_end,   na.rm = TRUE)
  b0 <- floor((gs - offset)/width) * width + offset
  b1 <- floor(((ge - eps) - offset)/width) * width + offset
  bks <- seq(b0, b1, by = width)
  cover <- pmax(0, pmin(ge, bks + width) - pmax(gs, bks)) / width
  tibble::tibble(horizon_bucket = as.numeric(bks), cover_frac = cover)
}

#' apply_edge_correction — Trim or rescale by coverage
#'
#' Applies edge treatment to a pooled per-bucket series:
#' - "trim": keep only buckets with coverage ≥ `min_frac`.
#' - "coverage": rescale by \eqn{k = 1/c(b)} where c(b) is coverage.
#'   Means and standard errors are multiplied by k (variances by k^2).
#' Leaves input unchanged for `method="none"`.
#' Should be tested and evaluated on test- and real-world-data - "trim" may 
#' often be the best option.
apply_edge_correction <- function(pooled,
                                  harm_tab,
                                  width,
                                  offset = 0,
                                  method = c("none","trim","coverage"),
                                  min_frac = 1.0) {
  method <- match.arg(method)
  if (method == "none") return(pooled)
  
  cf <- bucket_coverage_fraction(harm_tab, width, offset)
  out <- pooled %>% dplyr::left_join(cf, by = "horizon_bucket")
  
  if (method == "trim") {
    return(out %>% dplyr::filter(cover_frac >= min_frac))
  }
  
  if (method == "coverage") {
    # scale factor; guard tiny/NA
    out <- out %>%
      dplyr::mutate(k = dplyr::case_when(
        is.finite(cover_frac) & cover_frac > 0 ~ 1/cover_frac,
        TRUE ~ NA_real_
      ))
    
    # Keep originals for CI reconstruction
    has_Q  <- "Q_bar"     %in% names(pooled)
    has_lo <- "CI_lower"  %in% names(pooled)
    has_hi <- "CI_upper"  %in% names(pooled)
    has_se <- "se"        %in% names(pooled)
    
    if (has_Q) {
      orig_Q <- pooled$Q_bar
      out$Q_bar <- out$Q_bar * out$k
      if (has_lo && has_hi) {
        half_low  <- orig_Q - pooled$CI_lower
        half_high <- pooled$CI_upper - orig_Q
        out$CI_lower <- out$Q_bar - half_low  * out$k
        out$CI_upper <- out$Q_bar + half_high * out$k
      }
    }
    if (has_se) out$se <- out$se * out$k
    
    return(out %>% dplyr::select(-k))
  }
  
  out
}

# # Half-open bucket sequence aligned by width/offset
# bucket_seq_with_offset <- function(s, e, width, offset = 0, eps = 1e-9) {
#   if (any(is.na(c(s, e))) || e <= s) return(numeric(0))
#   b0 <- floor((s - offset)/width) * width + offset
#   b1 <- floor(((e - eps) - offset)/width) * width + offset
#   seq(b0, b1, by = width)
# }

#' row_overlap_long — Overlap length and weights per bucket
#'
#' For a single interval [eff_start, eff_end) and a bucket grid,
#' returns a long table of buckets with overlap length and normalized
#' weight \eqn{w_b = \text{overlap}_b / (e - s)}. Used for diagnostics.
row_overlap_long <- function(eff_start, eff_end, width, offset) {
  dur <- eff_end - eff_start
  if (!is.finite(dur) || dur <= 0) {
    return(tibble::tibble(horizon_bucket = numeric(0),
                          overlap_len = numeric(0),
                          weight_entry = numeric(0)))
  }
  bks <- bucket_seq_with_offset(eff_start, eff_end, width, offset)
  if (!length(bks)) {
    return(tibble::tibble(horizon_bucket = numeric(0),
                          overlap_len = numeric(0),
                          weight_entry = numeric(0)))
  }
  ov <- vapply(bks, function(b) max(0, min(eff_end, b + width) - max(eff_start, b)), numeric(1))
  tibble::tibble(horizon_bucket = as.numeric(bks),
                 overlap_len = as.numeric(ov),
                 weight_entry = as.numeric(ov / dur))
}


# ===============================
# Weight profiles (within-interval)
# ===============================

#' weight_profile — Within-interval temporal profile f(u)
#'
#' Returns a function \eqn{f(u)} on \eqn{u ∈ [0,1]} used to weight
#' within-interval mass. Options:
#' - uniform: \eqn{f(u)=1}
#' - triangular(peak p): ramps up to p then down
#' - beta(α,β): \eqn{f(u) ∝ \text{Beta}(u; α, β)} 
#' Profiles are later evaluated at overlap midpoints. All profiles operate on 
#' normalized position u in [0,1] from eff_start to eff_end. 
weight_profile <- function(type=c("uniform","triangular","beta"),
                           peak=0.5, alpha=2.5, beta=5) {
  type <- match.arg(type)
  switch(type,
         "uniform"   = function(u) rep(1, length(u)),
         "triangular"= function(u) {
           # peak in (0,1); linearly up to peak, then down
           p <- peak
           up   <- ifelse(u <= p, u/p, 0)
           down <- ifelse(u >  p, (1 - (u - p)/(1 - p)), 0)
           (up + down)
         },
         "beta"      = function(u) {
           # Beta(alpha, beta): early peak & long taper when alpha>1, beta>alpha
           dbeta(pmin(pmax(u, 0), 1), shape1=alpha, shape2=beta)
         }
  )
}

#' founding_rise_taper — Convenience beta profile
#'
#' Shorthand for a beta-shaped profile with early rise and long taper.
#' Used to model 'founding period' or 'object use period' style curves as 
#' theorised for archaeological material.
#' Equivalent to `weight_profile("beta", alpha, beta)`.
founding_rise_taper <- function(alpha=2.5, beta=6) weight_profile("beta", alpha=alpha, beta=beta)

# ===============================
# Prepare mortuary tables
# ===============================

#' prepare_mortuary_individual — Clean/standardize individual data
#'
#' Validates required columns, trims character fields, and tags the
#' source as "individual". Returns a tidy frame for downstream joins.
prepare_mortuary_individual <- function(indiv_df, debug=FALSE) {
  req <- c("ID","burial_id","sex_gender","age","phase_id","site_id")
  check_required(indiv_df, req, "individual mortuary")
  out <- indiv_df %>%
    mutate(sex_gender=str_trim(as.character(sex_gender)),
           age=str_trim(as.character(age)),
           source="individual")
  if (debug) message(sprintf("[indiv] rows=%d | unique burial_id=%d", nrow(out), n_distinct(out$burial_id)))
  out
}

#' prepare_mortuary_aggregated — Expand counts to pseudo-individuals
#'
#' Validates required columns, drops non-positive amounts, and expands
#' aggregated counts (`amount`) into pseudo-rows via `tidyr::uncount`.
#' Generates synthetic `burial_id`s and tags source="aggregated".
prepare_mortuary_aggregated <- function(agg_df, debug=FALSE) {
  req <- c("ID","site_id","sex_gender","age","amount","phase_id")
  check_required(agg_df, req, "aggregated mortuary")
  before <- nrow(agg_df)
  pruned <- agg_df %>% filter(!is.na(amount) & amount > 0L)
  if (debug && nrow(pruned) < before) message(sprintf("[aggreg] dropped %d rows with NA/<=0 amount", before - nrow(pruned)))
  expanded <- pruned %>%
    mutate(amount=as.integer(amount)) %>%
    tidyr::uncount(weights=amount, .remove=FALSE, .id="expand_seq") %>%
    mutate(sex_gender=str_trim(as.character(sex_gender)),
           age=str_trim(as.character(age)),
           burial_id=paste0("agg-", site_id, "-", ID, "-", expand_seq),
           source="aggregated") %>%
    select(ID, burial_id, sex_gender, age, phase_id, site_id, source)
  if (debug) message(sprintf("[aggreg] %d -> %d rows after expansion", before, nrow(expanded)))
  expanded
}

# ===============================
# Meta intervals (fallback)
# ===============================

#' site_meta_intervals — Derive interval bounds from meta table
#'
#' Builds site-level fallback intervals by mapping meta phase IDs (start/end)
#' through the harmonization table to numeric horizon bounds:
#' \eqn{[s,e] = [\min(\text{fade_in}, \text{start}), \max(\text{fade_out}, \text{end})]}.
#' Produces in-period and pre/post-period bounds for last-ditch meta fallback allocation.
site_meta_intervals <- function(meta_df, harm_tab) {
  req_meta <- c("site_id","cemetery_start_ia","cemetery_end_ia","cemetery_start_pre","cemetery_end_post")
  req_harm <- c("phase_id","horizon_start","horizon_end","fade_in_start","fade_out_end")
  check_required(meta_df, req_meta, "metadata for meta fallback")
  check_required(harm_tab, req_harm, "harmonization table")
  phase_to_interval <- function(pid) {
    if (is.na(pid)) return(c(NA_real_, NA_real_))
    hit <- harm_tab %>% filter(.data$phase_id==pid) %>%
      select(fade_in_start,horizon_start,fade_out_end,horizon_end)
    if (nrow(hit)!=1) return(c(NA_real_, NA_real_))
    s <- ifelse(!is.na(hit$fade_in_start), hit$fade_in_start, hit$horizon_start)
    e <- ifelse(!is.na(hit$fade_out_end),  hit$fade_out_end,  hit$horizon_end)
    c(min(s,e, na.rm=TRUE), max(s,e, na.rm=TRUE))
  }
  ia_bounds <- t(vapply(seq_len(nrow(meta_df)), function(i) {
    s <- meta_df$cemetery_start_ia[i]; e <- meta_df$cemetery_end_ia[i]
    if (is.na(s) || is.na(e)) return(c(NA_real_, NA_real_))
    rs <- phase_to_interval(s); re <- phase_to_interval(e)
    c(min(rs[1], re[1], na.rm=TRUE), max(rs[2], re[2], na.rm=TRUE))
  }, numeric(2)))
  pp_bounds <- t(vapply(seq_len(nrow(meta_df)), function(i) {
    s <- meta_df$cemetery_start_pre[i]; e <- meta_df$cemetery_end_post[i]
    if (is.na(s) || is.na(e)) return(c(NA_real_, NA_real_))
    rs <- phase_to_interval(s); re <- phase_to_interval(e)
    c(min(rs[1], re[1], na.rm=TRUE), max(rs[2], re[2], na.rm=TRUE))
  }, numeric(2)))
  meta_df %>% mutate(
    ia_start=ia_bounds[,1], ia_end=ia_bounds[,2],
    pp_start=pp_bounds[,1], pp_end=pp_bounds[,2]
  )
}

# ===============================
# Aoristic allocation core
# ===============================

#' bucket_seq_with_offset — Half-open bucket grid [s, e)
#'
#' Returns the sequence of bucket left-edges aligned to width `w` and
#' offset `o` that intersect the half-open interval [s, e). To ensure
#' consistent tiling, the right edge is nudged left by a tiny epsilon.
bucket_seq_with_offset <- function(s, e, width, offset = 0) {
  if (any(is.na(c(s, e)))) return(numeric(0))
  # Treat intervals as [s, e) by nudging e just left of the boundary
  eps <- 1e-9
  b0 <- floor((s - offset) / width) * width + offset
  b1 <- floor(((e - eps) - offset) / width) * width + offset
  seq(b0, b1, by = width)
}

#' weighted_overlaps — Profile-weighted overlaps per bucket at bucket centers
#'
#' For interval [s,e), bucket set {b}, and width w, computes overlap
#' length \eqn{\ell_b = | [b, b+w) ∩ [s,e) |}. The midpoint of the
#' overlap maps to \eqn{u_b = (m_b - s)/(e - s)}; the weight is
#' \eqn{w_b = \ell_b · f(u_b)} where f is the chosen profile.
#' Returns {w_b} with zeros for empty overlaps.
weighted_overlaps <- function(bks, s, e, width, f) {
  if (length(bks) == 0 || any(is.na(c(s, e))) || e <= s) return(numeric(0))
  
  # Overlap length per bucket
  lefts   <- pmax(s, bks)
  rights  <- pmin(e, bks + width)
  lens    <- pmax(0, rights - lefts)
  
  # Midpoint of the overlap segment (not the bucket center)
  mids    <- (lefts + rights) / 2
  u       <- (mids - s) / (e - s)           # normalize to [0,1] within the interval
  u       <- pmin(pmax(u, 0), 1)
  
  prof    <- f(u)
  # If no overlap or profile returns NA, treat as 0
  prof[!is.finite(prof)] <- 0
  
  # Profiled overlap; note: where lens==0 this is 0
  w <- lens * prof
  w[!is.finite(w)] <- 0
  w
}

#' choose_one_bucket — Deterministic or stochastic placement
#'
#' Given buckets {b} and nonnegative weights {w}, forms probabilities
#' \eqn{p_i = w_i / Σ_j w_j}. 
#' 
#' Modes:
#' - "deterministic_max": choose argmax p (earliest tie-break).
#' - "stochastic": sample one b ~ Multinomial(1, p) with a row-stable
#'   RNG scope (\code{withr::with_seed(digest(global_seed,row_id))}).
#'   
#' Returns list(bucket, p) for the selected bucket.
choose_one_bucket <- function(bks, wts, mode=c("deterministic_max","stochastic"),
                              global_seed=ocp_get_global_seed(), row_id=NULL) {
  mode <- match.arg(mode)
  if (length(bks)==0) return(list(bucket=NA_real_, p=NA_real_))
  wts[is.na(wts)] <- 0
  if (sum(wts) <= 0) wts <- rep(1, length(bks))
  p <- wts / sum(wts)
  if (mode=="deterministic_max") {
    mx <- max(p)
    sel <- which(p==mx)
    b  <- min(bks[sel])      # earliest tie-break
    return(list(bucket=b, p=mx))
  } else {
    # stochastic, localize RNG to avoid polluting global stream
    base <- paste0(global_seed, "::", if (is.null(row_id)) "" else as.character(row_id))
    withr::with_seed(digest::digest2int(base), {
      idx <- sample.int(length(bks), size=1, prob=p)
      list(bucket=bks[idx], p=p[idx])
    })
  }
}

# ===============================
# Harmonize & allocate to buckets
# ===============================

#' harmonize_chronologies — Main allocation engine
#'
#' Joins mortuary records to harmonization phases; builds effective
#' intervals \eqn{[s,e]} via fade-in/out vs start/end; constructs the
#' bucket grid (width+offset); allocates mass per record:
#' - "fractional": distribute \eqn{w_b / Σ w_b} across all overlapped buckets.
#' - "stochastic"/"deterministic_max": place into one bucket using
#'   `choose_one_bucket`.
#' Optional meta fallback uses site-level intervals with profile f_meta,
#' scaled by `meta_weight`, and capped by `meta_cap_fraction` per site.
#' Returns long table with per-record bucket weights plus site metadata.
harmonize_chronologies <- function(
    mort_indiv,
    mort_agg,
    harm_tab,
    meta_data,
    horizon_bracket_size = 10L,
    horizon_offset       = 0,
    sampling             = c("deterministic_max","stochastic","fractional"),
    seed                 = NULL,
    # within-interval profile for primary + meta
    primary_profile      = weight_profile("uniform"),
    meta_profile         = founding_rise_taper(alpha=2.5, beta=6),
    meta_fallback        = TRUE,
    meta_weight          = 0.25,
    meta_cap_fraction    = 0.30,   # cap meta share per site (0..1). NULL to disable
    debug                = FALSE
) {
  sampling <- match.arg(sampling)
  req_harm <- c("phase_id","system_name","horizon_start","horizon_end","fade_in_start","fade_out_end")
  req_meta_core <- c("site_id","site_name","coord_x","coord_y","cemetery_size","cemetery_region")
  check_required(harm_tab, req_harm, "harmonization table")
  check_required(meta_data, req_meta_core, "metadata core")
  
  mort_all <- bind_rows(
    prepare_mortuary_individual(mort_indiv %||% tibble(), debug=FALSE),
    prepare_mortuary_aggregated(mort_agg %||% tibble(), debug=FALSE)
  )
  if (nrow(mort_all)==0) stop("[harmonize] No mortuary rows provided.", call.=FALSE)
  check_required(mort_all, c("burial_id","phase_id","site_id","source"), "mortuary (combined)")
  
  # Join phases
  joined <- mort_all %>% left_join(harm_tab, by="phase_id")
  # Effective interval
  with_int <- joined %>% mutate(
    eff_start0 = coalesce(fade_in_start, horizon_start),
    eff_end0   = coalesce(fade_out_end,  horizon_end),
    eff_start  = pmin(eff_start0, eff_end0),
    eff_end    = pmax(eff_start0, eff_end0)
  )
  
  # PRIMARY allocation
  alloc_primary <- with_int %>%
    mutate(
      buckets       = pmap(list(eff_start, eff_end),
                           ~ bucket_seq_with_offset(..1, ..2, horizon_bracket_size, horizon_offset)),
      w_overlaps    = pmap(list(buckets, eff_start, eff_end),
                           ~ weighted_overlaps(..1, ..2, ..3, horizon_bracket_size, primary_profile)),
      num_buckets   = lengths(buckets)
    )
  
  if (sampling == "fractional") {
    normalize_w <- function(bks, w) {
      if (length(bks) == 0) return(numeric(0))
      if (length(w) == 0)   return(rep(1/length(bks), length(bks)))
      sw <- sum(w, na.rm = TRUE)
      if (!is.finite(sw) || sw <= 0) rep(1/length(bks), length(bks)) else (w / sw)
    }
    
    primary_long <- alloc_primary %>%
      mutate(
        p_vec = purrr::map2(buckets, w_overlaps, normalize_w),
        pairs = purrr::map2(buckets, p_vec,
                            ~ tibble::tibble(horizon_bucket = .x, weight_entry = .y))
      ) %>%
      select(burial_id, site_id, source, eff_start, eff_end, pairs) %>%
      tidyr::unnest(pairs) %>%
      dplyr::filter(!is.na(horizon_bucket) & weight_entry > 0) %>%
      dplyr::mutate(assign_source = "primary")
  } else {
    # single bucket per record
    picks <- pmap(list(alloc_primary$buckets, alloc_primary$w_overlaps, alloc_primary$burial_id),
                  ~ choose_one_bucket(..1, ..2, mode=sampling, global_seed=seed, row_id=..3))
    primary_long <- alloc_primary %>%
      mutate(horizon_bucket = map_dbl(picks, "bucket"),
             p_selected     = map_dbl(picks, "p"),
             weight_entry   = p_selected) %>%   # unbiased wrt fractional baseline
      filter(!is.na(horizon_bucket)) %>%
      transmute(burial_id, site_id, source, eff_start, eff_end, horizon_bucket, weight_entry,
                assign_source="primary")
  }
  
  # META fallback (only unassigned in non-fractional; in fractional, we add additional rows if primary missing)
  meta_rows <- tibble()
  if (meta_fallback) {
    meta_bounds <- site_meta_intervals(meta_data, harm_tab) %>% select(site_id, ia_start, ia_end, pp_start, pp_end)
    
    # Which rows need meta? For fractional: rows with zero buckets. For others: those dropped above.
    need_meta <- alloc_primary %>%
      mutate(has_primary = lengths(buckets) > 0) %>%
      filter(!has_primary) %>%
      select(burial_id, site_id, source, eff_start, eff_end)
    
    if (nrow(need_meta) > 0) {
      meta_aug <- need_meta %>%
        left_join(meta_bounds, by="site_id") %>%
        mutate(m_start = coalesce(ia_start, pp_start),
               m_end   = coalesce(ia_end,   pp_end)) %>%
        mutate(
          m_buckets    = pmap(list(m_start, m_end),
                              ~ bucket_seq_with_offset(..1, ..2, horizon_bracket_size, horizon_offset)),
          m_w_overlaps = pmap(list(m_buckets, m_start, m_end),
                              ~ weighted_overlaps(..1, ..2, ..3, horizon_bracket_size, meta_profile))
        )
      
      if (sampling == "fractional") {
        meta_rows <- meta_aug %>%
          mutate(
            p_vec = purrr::map2(m_buckets, m_w_overlaps, normalize_w),
            pairs = purrr::map2(m_buckets, p_vec,
                                ~ tibble::tibble(horizon_bucket = .x, w = .y))
          ) %>%
          select(burial_id, site_id, source, pairs) %>%
          tidyr::unnest(pairs) %>%
          dplyr::filter(!is.na(horizon_bucket) & w > 0) %>%
          dplyr::mutate(weight_entry = meta_weight * w,
                        assign_source = "meta") %>%
          dplyr::select(burial_id, site_id, source, horizon_bucket, weight_entry, assign_source)
      } else {
        picks_m <- pmap(list(meta_aug$m_buckets, meta_aug$m_w_overlaps, meta_aug$burial_id),
                        ~ choose_one_bucket(..1, ..2, mode=sampling, global_seed=seed, row_id=paste0(..3,"::meta")))
        meta_rows <- meta_aug %>%
          mutate(horizon_bucket = map_dbl(picks_m, "bucket"),
                 p_selected     = map_dbl(picks_m, "p"),
                 weight_entry   = meta_weight * p_selected) %>%
          filter(!is.na(horizon_bucket)) %>%
          transmute(burial_id, site_id, source, horizon_bucket, weight_entry, assign_source="meta")
      }
    }
  }
  
  # Combine
  combined <- bind_rows(primary_long, meta_rows)
  
  # Apply meta cap per site (overall across buckets)
  if (!is.null(meta_cap_fraction)) {
    site_totals <- combined %>%
      dplyr::group_by(site_id) %>%
      dplyr::summarise(
        primary = sum(weight_entry[assign_source == "primary"], na.rm = TRUE),
        meta    = sum(weight_entry[assign_source == "meta"],     na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        total = primary + meta,
        cap   = meta_cap_fraction * total,
        scale = dplyr::if_else(meta > 0 & meta > cap, cap / meta, 1)
      ) %>%
      dplyr::select(site_id, scale)
    
    combined <- combined %>%
      dplyr::left_join(site_totals, by = "site_id") %>%
      tidyr::replace_na(list(scale = 1)) %>%
      dplyr::mutate(
        weight_entry = dplyr::if_else(assign_source == "meta",
                                      weight_entry * scale, weight_entry)
      ) %>%
      dplyr::select(-scale)
  }
  
  
  # Attach site metadata core
  meta_core <- meta_data %>%
    select(site_id, site_name, coord_x, coord_y, cemetery_size, cemetery_region)
  final <- combined %>%
    left_join(meta_core, by="site_id") %>%
    mutate(horizon_bracket_size = horizon_bracket_size,
           horizon_offset       = horizon_offset)
  
  if (debug) {
    cat(sprintf("[harmonize] rows=%d | sites=%d | Pseudo Horizons=%d\n",
                nrow(final), n_distinct(final$site_id), n_distinct(final$horizon_bucket)))
    tab <- table(final$assign_source, useNA="ifany")
    cat(sprintf("[harmonize] sources: primary=%d | meta=%d\n",
                unname(tab["primary"] %||% 0L), unname(tab["meta"] %||% 0L)))
  }
  
  list(
    data       = final,
    drop_table = tibble(), # (drops are handled earlier; here we return placed rows)
    params     = list(
      bracket_size = horizon_bracket_size, offset = horizon_offset,
      sampling = sampling, seed = seed,
      meta_fallback = meta_fallback, meta_weight = meta_weight,
      meta_cap_fraction = meta_cap_fraction,
      profile_primary = deparse(substitute(primary_profile)),
      profile_meta    = deparse(substitute(meta_profile))
    )
  )
}


# ===============================
# Normalization helpers
# ===============================

#' %||% — Null coalescing helper
#'
#' Returns `a` unless it is NULL, otherwise returns `b`. Useful for
#' optional arguments and defaults.
`%||%` <- function(a,b) if (is.null(a)) b else a

#' compute_active_duration — Active PHU (Pseudo-Horizon-Units) per site
#'
#' For site×bucket data, computes \eqn{\text{duration\_phu} =
#' (\#\text{distinct buckets with mass}) × \text{bracket\_size}}.
compute_active_duration <- function(harmonized_data) {
  harmonized_data %>%
    group_by(site_id) %>%
    summarise(
      buckets = n_distinct(horizon_bucket),
      bracket_size = first(horizon_bracket_size),
      duration_phu = buckets * bracket_size,
      .groups="drop"
    )
}

#' add_normalized_entries — Normalize entries by size/duration
#'
#' Adds `normalized_entries = entries_per_period / denom`, where `denom`
#' is chosen by `normalization` - available modes:
#' - "size": cemetery_size
#' - "duration": duration_phu (computed from positive-mass buckets)
#' - "size_duration": product of size&duration
#' - "none": denom=1
#' Guards against zero/NA denominators: Computes duration_phu from the 
#' site×bucket summary using entries_per_period > 0
add_normalized_entries <- function(df_site_bucket,
                                   normalization = c("size","duration","size_duration","none")) {
  normalization <- match.arg(normalization)
  df <- df_site_bucket
  
  if (normalization == "none") {
    return(df %>%
             dplyr::mutate(
               denom = 1,
               normalized_entries = entries_per_period
             ))
  }
  ## SIZE
  if (normalization == "size") {
    return(df %>%
             dplyr::mutate(
               denom = cemetery_size,
               normalized_entries = dplyr::if_else(denom > 0, entries_per_period / denom, NA_real_)
             ))
  }
  
  # PHU duration derivation from positive-mass buckets
  req <- c("site_id", "horizon_bucket", "bracket_size", "entries_per_period")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("[add_normalized_entries] Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  dur_tbl <- df %>%
    dplyr::filter(is.finite(entries_per_period), entries_per_period > 0) %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(
      buckets      = dplyr::n_distinct(horizon_bucket),
      bracket_size = dplyr::first(bracket_size),
      duration_phu = buckets * bracket_size,
      .groups = "drop"
    )
  
  df <- df %>% dplyr::left_join(dur_tbl, by = "site_id")
  
  ## DURATION
  if (normalization == "duration") {
    return(df %>%
             dplyr::mutate(
               denom = duration_phu,
               normalized_entries = dplyr::if_else(denom > 0, entries_per_period / denom, NA_real_)
             ))
  }
  
  ## SIZE_DURATION
  if (normalization == "size_duration") {
    return(df %>%
             dplyr::mutate(
               denom = cemetery_size * duration_phu,
               normalized_entries = dplyr::if_else(denom > 0, entries_per_period / denom, NA_real_)
             ))
  }
  
  stop("[add_normalized_entries] Unknown normalization: ", normalization, call. = FALSE)
}

# ===============================
# Analysis 
# ===============================

#' analyze_occupancy — Aggregate site×bucket and compute normalization
#'
#' Aggregates assigned weights to site×bucket totals
#' (`entries_per_period`), carries site covariates, computes duration,
#' and adds `normalized_entries`. Also returns a diagnostic total mass 
#' per bucket (primary+meta).
analyze_occupancy <- function(harmonized_data,
                              normalization = c("size","duration","size_duration","none"),
                              horizon_bracket_size,  # kept for interface symmetry
                              horizon_offset = 0,
                              debug = FALSE) {
  normalization <- match.arg(normalization)
  data <- harmonized_data
  
  # Aggregate existing bucket allocations; no re-overlap here.
  req_cols <- c("site_id","cemetery_size","coord_x","coord_y",
                "horizon_bucket","weight_entry","assign_source")
  miss <- setdiff(req_cols, names(data))
  if (length(miss)) stop("[analyze_occupancy] Missing columns: ", paste(miss, collapse = ", "))
  
  has_region <- "cemetery_region" %in% names(data)
  has_bsize  <- "horizon_bracket_size" %in% names(data)
  
  # Ensure numeric buckets
  data <- data %>% dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket))
  
  # Site × bucket: total mass and split by assign_source ("primary"/"meta")
  site_bucket <- data %>%
    dplyr::group_by(site_id, horizon_bucket) %>%
    dplyr::summarise(
      entries_per_period = sum(weight_entry, na.rm = TRUE),
      primary_share      = sum(dplyr::if_else(assign_source == "primary", weight_entry, 0), na.rm = TRUE),
      meta_share         = sum(dplyr::if_else(assign_source == "meta",    weight_entry, 0), na.rm = TRUE),
      cemetery_size      = dplyr::first(cemetery_size),
      coord_x            = dplyr::first(coord_x),
      coord_y            = dplyr::first(coord_y),
      !!(if (has_region) rlang::sym("cemetery_region") else rlang::sym("tmp")) :=
        if (has_region) dplyr::first(.data$cemetery_region) else NA_character_,
      bracket_size       = if (has_bsize) dplyr::first(horizon_bracket_size) else horizon_bracket_size,
      .groups = "drop"
    ) %>%
    { if (has_region) dplyr::rename(., cemetery_region = `cemetery_region`) else
      dplyr::select(., -tmp) }
  
  # Per-site active duration (PHU) = (#distinct positive-mass buckets) × width
  dur_tbl <- site_bucket %>%
    dplyr::filter(entries_per_period > 0) %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(duration_phu = dplyr::n_distinct(horizon_bucket) * dplyr::first(bracket_size),
                     .groups = "drop")
  
  site_bucket <- site_bucket %>%
    dplyr::left_join(dur_tbl, by = "site_id") %>%
    dplyr::mutate(
      denom = dplyr::case_when(
        normalization == "size"          ~ cemetery_size,
        normalization == "duration"      ~ duration_phu,
        normalization == "size_duration" ~ cemetery_size * duration_phu,
        TRUE                             ~ 1
      ),
      normalized_entries = dplyr::if_else(denom > 0, entries_per_period / denom, NA_real_)
    )
  
  # Horizon-level total (for diagnostics), includes meta via entries_per_period
  horizon_source_table <- site_bucket %>%
    dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(total = sum(entries_per_period, na.rm = TRUE), .groups = "drop")
  
  diag_plot_sources <- ggplot2::ggplot(horizon_source_table,
                                       ggplot2::aes(x = as.numeric(horizon_bucket), y = total)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Pseudo Horizon (bucket)", y = "Total weight",
                  title = "Total (primary+meta) mass per bucket") +
    ggplot2::theme_minimal()
  
  list(
    site_bucket = site_bucket,
    horizon_source_table = horizon_source_table,
    diag_plot_sources = diag_plot_sources
  )
}


#' plot_norm_vs_size — Size effect diagnostic
#'
#' Scatter of `normalized_entries` vs cemetery size with LOESS smoother.
#' Uses log10 x-scale to reveal systematic over/under-normalization.
plot_norm_vs_size <- function(site_bucket) {
  ggplot(site_bucket, aes(x=cemetery_size, y=normalized_entries)) +
    geom_point(alpha=0.4) +
    geom_smooth(method="loess", se=TRUE, span=0.8) +
    scale_x_continuous(trans="log10") +
    labs(x="Cemetery size (log10)", y="Normalized entries",
         title="Normalized entries vs Cemetery size") +
    theme_minimal()
}

#' fit_gam_with_size_covariate — Model entries ~ f(H) + f(log size)
#'
#' Fits \eqn{E[y] = s(H) + s(\log(1+size))} using `mgcv::gam` with REML.
#' Tests whether cemetery size explains residual structure beyond time
#' (H = numeric bucket). Useful to check adequacy of simple division.
fit_gam_with_size_covariate <- function(site_bucket) {
  # mean entries per site×bucket; use log(size) covariate and smooth over horizon
  dat <- site_bucket %>% mutate(log_size = log1p(cemetery_size),
                                H = as.numeric(horizon_bucket))
  mgcv::gam(entries_per_period ~ s(H, k=10) + s(log_size, k=5),
            data=dat, family=gaussian(), method="REML")
}

#' within_variance_block_boot — Region block bootstrap (within-rep Var)
#'
#' Estimates, for each bucket H, the variance of the regional average
#' \eqn{\hat μ(H)} by resampling regions with replacement (counts = mult),
#' then recomputing the region-weighted mean. The empirical variance
#' across B bootstrap draws approximates within-replicate uncertainty
#' from regional composition. Returns tibble(horizon_bucket, var_hat).
within_variance_block_boot <- function(site_bucket,
                                       region_col = "cemetery_region",
                                       B = 200L,
                                       seed = ocp_get_global_seed()) {
  reg_sym <- rlang::sym(region_col)
  req <- c("horizon_bucket", "normalized_entries", region_col)
  miss <- setdiff(req, names(site_bucket))
  if (length(miss)) stop("[within_variance_block_boot] Missing: ", paste(miss, collapse=", "))
  
  sb <- site_bucket %>%
    dplyr::mutate(H = as.numeric(horizon_bucket)) %>%
    dplyr::select(!!reg_sym, H, normalized_entries)
  
  regions <- sb %>% dplyr::distinct(!!reg_sym) %>% dplyr::pull(1)
  Hvals   <- sort(unique(sb$H))
  regH <- sb %>%
    dplyr::group_by(!!reg_sym, H) %>%
    dplyr::summarise(mean_norm = mean(normalized_entries, na.rm = TRUE), .groups="drop")
  
  boot_once <- function() {
    samp   <- sample(regions, size = length(regions), replace = TRUE)
    counts <- tibble::tibble(!!reg_sym := samp) %>% dplyr::count(!!reg_sym, name="mult")
    merged <- regH %>% dplyr::left_join(counts, by = region_col)
    merged %>% dplyr::group_by(H) %>%
      dplyr::summarise(
        mu = ifelse(sum(mult, na.rm=TRUE) > 0,
                    sum(mean_norm * mult, na.rm=TRUE) / sum(mult, na.rm=TRUE),
                    NA_real_),
        .groups="drop")
  }
  
  boot_mat <- matrix(NA_real_, nrow=length(Hvals), ncol=B, dimnames=list(H=Hvals, b=seq_len(B)))
  withr::with_seed(seed %||% ocp_derive_seed("boot"), {
    for (b in seq_len(B)) {
      boot <- boot_once()
      boot_mat[match(boot$H, Hvals), b] <- boot$mu
    }
  })
  var_hat <- apply(boot_mat, 1, stats::var, na.rm = TRUE)
  tibble::tibble(horizon_bucket = as.numeric(names(var_hat)), var_hat = as.numeric(var_hat))
}

#' run_one_replicate — One full run (harmonize→analyze→variance)
#'
#' Sets RNG to a given `seed`, harmonizes chronologies, aggregates to
#' site×bucket (`analyze_occupancy`), computes per-bucket mean of
#' `normalized_entries` and within variance via block bootstrap, and
#' applies requested edge correction (trim/coverage/none).
#' Also returns diagnostics, coverage table, and manifest of settings.
run_one_replicate <- function(seed = ocp_get_global_seed(),
                              mort_indiv, mort_agg, harm_tab, meta_data,
                              horizon_bracket_size = 10L,
                              horizon_offset       = 0,
                              sampling             = c("deterministic_max","stochastic","fractional"),
                              primary_profile      = weight_profile("uniform"),
                              meta_profile         = founding_rise_taper(),
                              meta_weight          = 0.25,
                              meta_cap_fraction    = 0.30,
                              normalization        = c("size","duration","size_duration","none"),
                              region_col           = "cemetery_region",
                              B_blocks             = 200L,
                              edge_correction      = c("none","trim","coverage"),
                              min_coverage_frac    = 1.0,
                              debug                = FALSE) {
  
  sampling        <- match.arg(sampling)
  normalization   <- match.arg(normalization)
  edge_correction <- match.arg(edge_correction)
  
  # ---- Harmonize & analyze ----
  withr::with_seed(seed, {
    harm_out <- harmonize_chronologies(
      mort_indiv          = mort_indiv,
      mort_agg            = mort_agg,
      harm_tab            = harm_tab,
      meta_data           = meta_data,
      horizon_bracket_size= horizon_bracket_size,
      horizon_offset      = horizon_offset,
      sampling            = sampling,
      seed                = seed,
      primary_profile     = primary_profile,
      meta_profile        = meta_profile,
      meta_fallback       = TRUE,
      meta_weight         = meta_weight,
      meta_cap_fraction   = meta_cap_fraction,
      debug               = debug
    )
    
    ana <- analyze_occupancy(
      harmonized_data      = harm_out$data,
      normalization        = normalization,
      horizon_bracket_size = horizon_bracket_size,
      horizon_offset       = horizon_offset,
      debug                = FALSE
    )
  # Per-bucket mean and within-replicate variance via region block bootstrap
    
    sb  <- ana$site_bucket
    per_bucket <- sb %>%
      dplyr::group_by(horizon_bucket) %>%
      dplyr::summarise(mean_entries = mean(normalized_entries, na.rm=TRUE),
                       n_sites = dplyr::n(), .groups="drop")

    vhat <- within_variance_block_boot(sb, region_col = region_col, B = B_blocks,
                                       seed = ocp_derive_seed("boot", seed))
    per_bucket <- per_bucket %>% dplyr::left_join(vhat, by = "horizon_bucket")
    per_bucket$var_hat[is.na(per_bucket$var_hat)] <- 0  
    
  # ---- Edge handling ---------------------------------------------
  # Coverage (fraction of each bucket inside the global PHU domain)
  coverage_tbl <- bucket_coverage_fraction(
    harm_tab = harm_tab,
    width    = horizon_bracket_size,
    offset   = horizon_offset
  )
  
  per_bucket_raw <- per_bucket  # keep unmodified for transparency
  
  per_bucket <- per_bucket %>%
    dplyr::left_join(coverage_tbl, by = "horizon_bucket")
  
  if (edge_correction == "trim") {
    kept_before <- nrow(per_bucket)
    per_bucket  <- per_bucket %>%
      dplyr::filter(cover_frac >= min_coverage_frac)
    if (debug) message(
      sprintf("[run_one_replicate] edge=trim; kept %d/%d buckets (min_frac=%.3f)",
              nrow(per_bucket), kept_before, min_coverage_frac)
    )
    per_bucket <- per_bucket %>% dplyr::select(-cover_frac)
  } else if (edge_correction == "coverage") {
    # scale factor k = 1/cover_frac (guard tiny/NA)
    per_bucket <- per_bucket %>%
      dplyr::mutate(
        k = dplyr::case_when(is.finite(cover_frac) & cover_frac > 0 ~ 1/cover_frac,
                             TRUE ~ NA_real_)
      )
    # Scale the mean and variance accordingly
    per_bucket <- per_bucket %>%
      dplyr::mutate(
        mean_entries = mean_entries * k,
        var_hat      = var_hat * (k^2)
      ) %>%
      dplyr::select(-k, -cover_frac)
    if (debug) message("[run_one_replicate] edge=coverage; rescaled means/variances by 1/coverage")
  } else {
    # edge_correction == "none" → drop helper column
    per_bucket <- per_bucket %>% dplyr::select(-cover_frac)
  }
  
  # ---- Return --------------------
  list(
    per_bucket       = per_bucket,        # corrected (or unchanged if 'none')
    per_bucket_raw   = per_bucket_raw,    # raw, pre-correction
    diagnostics      = list(
      horizon_source_table = ana$horizon_source_table,
      diag_plot_sources    = ana$diag_plot_sources,
      site_bucket          = sb
    ),
    coverage_tbl     = coverage_tbl,      # bucket → coverage fraction
    edge_correction  = edge_correction,
    min_coverage_frac = if (edge_correction == "trim") min_coverage_frac else NA_real_,
    horizon_bracket_size = horizon_bracket_size,
    horizon_offset       = horizon_offset,
    sampling             = sampling,
    normalization        = normalization,
    seed                 = seed
  )
  })
}



#' pool_replicates — Rubin-style pooling across M runs
#'
#' Combines per-bucket estimates from M replicates:
#' \deqn{Q̄ = mean(Q_m),\quad W = mean(W_m),\quad B = var(Q_m),}
#' \deqn{T = W + (1 + 1/M)B,\quad SE = √T.}
#' Degrees of freedom (Barnard–Rubin):
#' \deqn{df = (M-1)\left(1 + \frac{W}{(1 + 1/M)B}\right)^2.}
#' 95% CI uses t-quantiles with `df` (or z=1.96 if `use_t_df=FALSE`).
pool_replicates <- function(reps_list, use_t_df=TRUE) {
  all <- dplyr::bind_rows(
    purrr::imap(reps_list, ~ dplyr::mutate(.x$per_bucket, .rep = .y))
  )
  pooled <- all %>%
    group_by(horizon_bucket) %>%
    summarise(
      M      = n(),
      Q_bar  = mean(mean_entries, na.rm=TRUE),
      W      = mean(var_hat,      na.rm=TRUE),
      B      = stats::var(mean_entries, na.rm=TRUE),
      T      = W + (1 + 1/M) * B,
      df     = ifelse(B > 0,
                      (M - 1) * (1 + W/((1 + 1/M) * B))^2,
                      Inf),
      SE     = sqrt(T),
      crit   = ifelse(use_t_df, qt(0.975, df=pmax(df,1)), 1.96),
      CI_lower = Q_bar - crit * SE,
      CI_upper = Q_bar + crit * SE,
      .groups="drop"
    )
  pooled
}

#' run_stochastic_ensemble — Grid of (width, offset) with M replicates
#'
#' For each (bucket width, offset) pair, runs `M` replicates with seeds
#' drawn from `base_seed`, pools them (Rubin), applies edge correction,
#' and also pools regional curves (facets). Returns per-setting results
#' and ready-to-plot ggplots for overall and regional trends.
#' NOTE: global RNG is seeded once for seed generation; per-replicate
#' runs use their own (derived) fixed seeds for reproducibility.
run_stochastic_ensemble <- function(M,
                                    mort_indiv, mort_agg, harm_tab, meta_data,
                                    widths            = 10L,
                                    offsets           = 0,
                                    sampling          = c("deterministic_max","stochastic","fractional"),
                                    base_seed         = ocp_get_global_seed(),
                                    primary_profile   = weight_profile("uniform"),
                                    meta_profile      = founding_rise_taper(),
                                    meta_weight       = 0.25,
                                    meta_cap_fraction = 0.30,
                                    normalization     = c("size","duration","size_duration","none"),
                                    region_col        = "cemetery_region",
                                    B_blocks          = 200L,
                                    use_t_df          = TRUE,
                                    edge_correction   = c("none","trim","coverage"),
                                    min_coverage_frac = 1.0,
                                    # OSI thresholds
                                    thr_green         = 0.05,
                                    thr_orange        = 0.10,
                                    # --- Regional facets ---
                                    do_regional       = TRUE,
                                    regional_loess_span = 0.75,
                                    # --- Spatial slices ---
                                    do_spatial          = TRUE,
                                    spatial_selected_buckets = c(30, 60, 90, 120),
                                    spatial_use_normalized   = FALSE,
                                    spatial_bins        = 7,
                                    spatial_bbox        = c(7, 15, 51, 58),
                                    spatial_point_color = "red3",
                                    spatial_point_alpha = 0.6,
                                    spatial_anim        = TRUE,
                                    spatial_anim_bucket_bin = 10) {
  
  sampling        <- match.arg(sampling)
  normalization   <- match.arg(normalization)
  edge_correction <- match.arg(edge_correction)
  
  # ensure integer grid (prevents type/unique() hiccups downstream)
  widths  <- as.integer(widths)
  offsets <- as.integer(offsets)
  grid    <- tidyr::expand_grid(width = widths, offset = offsets)
  
  # spatial background once (if requested)
  if (do_spatial) {
    if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
        !requireNamespace("sf", quietly = TRUE)) {
      stop("Packages 'rnaturalearth' and 'sf' are required for spatial maps.")
    }
  }
  world_sf <- if (do_spatial) rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") else NULL
  
  # ------------- flattened triplets: (width, offset, rep_idx, seed) -------------
  triplets <- tidyr::crossing(grid, rep_idx = seq_len(M)) %>%
    dplyr::mutate(
      seed = digest::digest2int(paste(base_seed, width, offset, rep_idx, sep = "::"))
    )
  
  # worker for one triplet (explicitly force edge_correction="none" inside replicate)
  .one_triplet <- function(width, offset, seed) {
    run_one_replicate(
      seed                 = seed,
      mort_indiv           = mort_indiv,
      mort_agg             = mort_agg,
      harm_tab             = harm_tab,
      meta_data            = meta_data,
      horizon_bracket_size = as.integer(width),
      horizon_offset       = as.integer(offset),
      sampling             = sampling,
      primary_profile      = primary_profile,
      meta_profile         = meta_profile,
      meta_weight          = meta_weight,
      meta_cap_fraction    = meta_cap_fraction,
      normalization        = normalization,
      region_col           = region_col,
      B_blocks             = B_blocks,
      edge_correction      = "none",
      min_coverage_frac    = 1.0,
      debug                = FALSE
    )
  }
  
  # map over triplets (parallel if furrr available)
  use_furrr <- requireNamespace("furrr", quietly = TRUE)
  if (use_furrr) {
    res_list <- furrr::future_pmap(
      list(triplets$width, triplets$offset, triplets$seed),
      ~ .one_triplet(..1, ..2, ..3),
      .progress = FALSE
    )
  } else {
    res_list <- purrr::pmap(
      list(triplets$width, triplets$offset, triplets$seed),
      ~ .one_triplet(..1, ..2, ..3)
    )
  }
  
  # attach replicate outputs to triplets
  trip_out <- dplyr::bind_cols(triplets, tibble::tibble(.rep_out = res_list))
  
  # ------------- regroup by (width, offset) and build per-cell results -------------
  results <- trip_out %>%
    dplyr::group_by(width, offset) %>%
    dplyr::group_map(function(df, keys) {
      reps <- df$.rep_out
      
      # pool on raw replicate series, then apply chosen edge-correction once
      reps_for_pool <- lapply(reps, function(r) list(per_bucket = if ("per_bucket_raw" %in% names(r)) r$per_bucket_raw else r$per_bucket))
      pooled_raw <- pool_replicates(reps_for_pool, use_t_df = use_t_df)
      pooled_corr <- apply_edge_correction(
        pooled   = pooled_raw,
        harm_tab = harm_tab,
        width    = keys$width,
        offset   = keys$offset,
        method   = edge_correction,
        min_frac = min_coverage_frac
      )
      
      coverage_tbl <- bucket_coverage_fraction(harm_tab, keys$width, keys$offset)
      
      # overall plot
      edge_note <- switch(edge_correction,
                          none     = "edge: none",
                          trim     = paste0("edge: trim (min=", signif(min_coverage_frac, 3), ")"),
                          coverage = "edge: coverage")
      overall_plot <- ggplot2::ggplot(pooled_corr,
                                      ggplot2::aes(x = as.numeric(horizon_bucket), y = Q_bar)) +
        ggplot2::geom_point() +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
        ggplot2::geom_smooth(method = "loess", se = TRUE, span = 0.75) +
        ggplot2::labs(
          x = "Pseudo Horizon (bucket index)",
          y = "Normalized occupancy (pooled)",
          title = "Pooled trend across Pseudo Horizons",
          subtitle = paste0("M=", nrow(df), " · width=", keys$width, " · offset=", keys$offset,
                            " · ", edge_note, " · n=", nrow(pooled_corr))
        ) +
        ggplot2::theme_minimal()
      
      # regional pooling & plot (interior only)
      region_pooled <- NULL
      region_plot   <- NULL
      if (isTRUE(do_regional)) {
        ib <- interior_mask(coverage_tbl, min_frac = min_coverage_frac)
        reg_long <- purrr::imap_dfr(reps, function(rep, irep) {
          sb <- rep$diagnostics$site_bucket %>%
            dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
            dplyr::inner_join(ib, by = "horizon_bucket") %>%
            dplyr::filter(is_interior)
          if (!nrow(sb) || !region_col %in% names(sb)) return(tibble::tibble())
          sb %>%
            dplyr::group_by(.data[[region_col]], horizon_bucket) %>%
            dplyr::summarise(
              mean_entries = mean(normalized_entries, na.rm = TRUE),
              n_sites      = dplyr::n_distinct(site_id),
              var_hat_reg  = {
                v <- stats::var(normalized_entries, na.rm = TRUE)
                if (!is.finite(v)) 0 else v / pmax(n_sites, 1)
              },
              .groups = "drop"
            ) %>%
            dplyr::transmute(region = as.character(.data[[region_col]]),
                             horizon_bucket = as.numeric(horizon_bucket),
                             mean_entries, var_hat_reg, .rep = irep)
        })
        
        if (nrow(reg_long)) {
          region_pooled <- reg_long %>%
            dplyr::group_by(region, horizon_bucket) %>%
            dplyr::summarise(
              M     = dplyr::n(),
              Q_bar = mean(mean_entries, na.rm = TRUE),
              W     = mean(var_hat_reg,  na.rm = TRUE),
              B     = stats::var(mean_entries, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            dplyr::mutate(
              T  = W + (1 + 1/M) * B,
              df = dplyr::if_else(B > 0,
                                  (M - 1) * (1 + W/((1 + 1/M) * B))^2,
                                  Inf),
              SE = sqrt(T),
              crit = if (use_t_df) stats::qt(0.975, df = pmax(df, 1)) else 1.96,
              CI_lower = Q_bar - crit * SE,
              CI_upper = Q_bar + crit * SE
            ) %>%
            dplyr::left_join(coverage_tbl, by = "horizon_bucket") %>%
            { if (edge_correction == "coverage") {
              dplyr::mutate(.,
                            k = dplyr::if_else(is.finite(cover_frac) & cover_frac > 0, 1/cover_frac, NA_real_),
                            Q_bar = Q_bar * k, SE = SE * k,
                            CI_lower = Q_bar - crit * SE, CI_upper = Q_bar + crit * SE
              ) %>% dplyr::select(-k, -cover_frac)
            } else dplyr::select(., -cover_frac)
            }
          
          region_plot <- ggplot2::ggplot(region_pooled,
                                         ggplot2::aes(x = as.numeric(horizon_bucket), y = Q_bar)) +
            ggplot2::geom_point(size = 1.6) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, alpha = 0.7) +
            ggplot2::geom_smooth(method = "loess", se = TRUE, span = regional_loess_span) +
            ggplot2::facet_wrap(~ region, scales = "free_y") +
            ggplot2::labs(
              title = "Regional temporal trends (pooled across replicates)",
              subtitle = paste0("width=", keys$width, " · offset=", keys$offset, " · ", edge_note),
              x = "Pseudo Horizon (bucket)", y = "Normalized occupancy"
            ) +
            ggplot2::theme_minimal()
        }
      }
      
      # spatial slices (per (width, offset)) — same behaviour as your current version
      spatial <- NULL
      if (isTRUE(do_spatial)) {
        sb <- reps[[1]]$diagnostics$site_bucket
        need <- c("horizon_bucket","entries_per_period","normalized_entries","coord_x","coord_y","site_id","cemetery_size")
        miss <- setdiff(need, names(sb))
        if (!length(miss)) {
          val_col <- if (isTRUE(spatial_use_normalized)) "normalized_entries" else "entries_per_period"
          sb_export <- sb %>%
            dplyr::select(site_id, entries_per_period, normalized_entries, coord_x, coord_y, cemetery_size) %>%
            dplyr::distinct() %>%
            dplyr::left_join(meta_data %>% dplyr::select(site_id, site_name), by = "site_id") %>%
            dplyr::relocate(site_name, .after = site_id)
          
          slice_base <- sb %>%
            dplyr::filter(.data$horizon_bucket %in% spatial_selected_buckets) %>%
            dplyr::filter(is.finite(.data$coord_x), is.finite(.data$coord_y)) %>%
            dplyr::mutate(value = .data[[val_col]])
          
          aes_density <- ggplot2::aes(x = .data$coord_x, y = .data$coord_y,
                                      fill = after_stat(level), alpha = after_stat(level))
          slices <- list()
          for (hb in spatial_selected_buckets) {
            d <- dplyr::filter(slice_base, .data$horizon_bucket == hb)
            if (!nrow(d)) next
            p <- ggplot2::ggplot() +
              ggplot2::geom_sf(data = world_sf, fill = "gray95", color = "gray80") +
              ggplot2::stat_density_2d(
                data = d, mapping = aes_density, geom = "polygon",
                contour = TRUE, contour_var = "ndensity", bins = spatial_bins,
                show.legend = c(fill = TRUE, alpha = FALSE, size = FALSE)
              ) +
              ggplot2::scale_fill_viridis_c(option = "C") +
              ggplot2::scale_alpha(range = c(0.2, 0.7), guide = "none") +
              ggplot2::geom_point(
                data = d,
                ggplot2::aes(x = .data$coord_x, y = .data$coord_y, size = .data$value),
                color = spatial_point_color, alpha = spatial_point_alpha,
                show.legend = c(size = TRUE)
              ) +
              ggplot2::scale_size_continuous(
                name = if (spatial_use_normalized) "Normalized occupancy" else "Expected burials",
                guide = ggplot2::guide_legend(
                  override.aes = list(shape = 19, colour = spatial_point_color, fill = "white", alpha = 1, stroke = 0)
                )
              ) +
              ggplot2::labs(title = sprintf("Horizon Bucket: %s", hb),
                            x = "Longitude", y = "Latitude", fill = "Density level") +
              ggplot2::coord_sf(xlim = c(spatial_bbox[1], spatial_bbox[2]),
                                ylim = c(spatial_bbox[3], spatial_bbox[4]), expand = FALSE) +
              ggplot2::guides(alpha = "none", fill = ggplot2::guide_colorbar(order = 1)) +
              ggplot2::theme_minimal() +
              ggplot2::theme(
                legend.position       = "right",
                legend.background     = ggplot2::element_rect(fill = "white", colour = NA),
                legend.box.background = ggplot2::element_rect(fill = "white", colour = NA),
                legend.key            = ggplot2::element_rect(fill = "white", colour = NA)
              )
            slices[[as.character(hb)]] <- p
          }
          
          combined <- NULL
          if (length(slices) > 1 &&
              requireNamespace("cowplot", quietly = TRUE) &&
              requireNamespace("patchwork", quietly = TRUE)) {
            with_leg <- lapply(slices, function(p) {
              leg   <- cowplot::get_legend(p)
              panel <- p + ggplot2::theme(legend.position = "none")
              cowplot::plot_grid(panel, leg, ncol = 2, rel_widths = c(1, 0.35), align = "h", axis = "tb")
            })
            combined <- patchwork::wrap_plots(with_leg, ncol = if (length(with_leg) == 4) 2 else min(3, length(with_leg))) +
              patchwork::plot_annotation(title = "Spatial density by horizon bucket")
          }
          
          anim_plot <- NULL
          if (isTRUE(spatial_anim) && requireNamespace("gganimate", quietly = TRUE)) {
            anim_base <- sb %>%
              dplyr::filter(is.finite(.data$coord_x), is.finite(.data$coord_y)) %>%
              dplyr::mutate(value = .data[[val_col]],
                            combined_bucket = floor(.data$horizon_bucket / spatial_anim_bucket_bin) * spatial_anim_bucket_bin) %>%
              dplyr::group_by(.data$combined_bucket) %>%
              dplyr::filter(dplyr::n() > 1) %>%
              dplyr::ungroup()
            if (nrow(anim_base)) {
              anim_plot <- ggplot2::ggplot(anim_base, ggplot2::aes(x = .data$coord_x, y = .data$coord_y)) +
                ggplot2::geom_sf(data = world_sf, fill = "gray95", color = "gray80", inherit.aes = FALSE) +
                ggplot2::coord_sf(xlim = c(spatial_bbox[1], spatial_bbox[2]),
                                  ylim = c(spatial_bbox[3], spatial_bbox[4]), expand = FALSE) +
                ggplot2::stat_density_2d(
                  ggplot2::aes(fill = after_stat(level), alpha = after_stat(level)),
                  geom = "polygon", contour = TRUE, contour_var = "ndensity", bins = spatial_bins,
                  show.legend = c(fill = TRUE, alpha = FALSE, size = FALSE)
                ) +
                ggplot2::scale_fill_viridis_c(option = "C") +
                ggplot2::scale_alpha(range = c(0.1, 0.9), guide = "none") +
                ggplot2::geom_point(
                  ggplot2::aes(size = .data$value),
                  color = spatial_point_color, alpha = spatial_point_alpha,
                  show.legend = c(size = TRUE)
                ) +
                ggplot2::scale_size_continuous(
                  name = if (spatial_use_normalized) "Normalized occupancy" else "Expected burials",
                  guide = ggplot2::guide_legend(
                    override.aes = list(shape = 19, colour = spatial_point_color, fill = "white", alpha = 1, stroke = 0)
                  )
                ) +
                ggplot2::labs(title = "Horizon Bucket: {closest_state}",
                              x = "Longitude", y = "Latitude", fill = "Density level") +
                ggplot2::guides(alpha = "none", fill = ggplot2::guide_colorbar(order = 1)) +
                ggplot2::theme_minimal() +
                ggplot2::theme(
                  legend.key           = ggplot2::element_rect(fill = "white", colour = NA),
                  legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
                  legend.box.background= ggplot2::element_rect(fill = "white", colour = NA)
                ) +
                gganimate::transition_states(as.factor(.data$combined_bucket),
                                             transition_length = 15, state_length = 30) +
                gganimate::ease_aes("circular-in-out")
            }
          }
          
          spatial <- list(
            slices             = slices,
            combined           = combined,
            anim_plot          = anim_plot,
            site_bucket_export = sb_export
          )
        }
      }
      
      list(
        width        = as.integer(keys$width),
        offset       = as.integer(keys$offset),
        seeds        = df$seed,
        reps         = reps,
        pooled_raw   = pooled_raw,
        pooled       = pooled_corr,
        coverage_tbl = coverage_tbl,
        overall_plot = overall_plot,
        region_pooled= region_pooled,
        region_plot  = region_plot,
        spatial      = spatial
      )
    })
  
  # convenience lists
  plots        <- purrr::map(results, "overall_plot")
  region_plots <- purrr::map(results, "region_plot")
  
  # ------------- OSI (rank-aligned) + LOSO medians, per width -------------
  by_width      <- split(results, vapply(results, function(r) as.character(r$width), character(1)))
  diag_osi      <- list()
  diag_loso_med <- list()
  
  for (w in names(by_width)) {
    rr   <- by_width[[w]]
    wnum <- as.numeric(w)
    
    # rank-aligned long table across offsets
    long <- purrr::map_dfr(rr, function(r) {
      ib <- interior_mask(r$coverage_tbl, min_frac = min_coverage_frac)
      r$pooled %>%
        dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
        dplyr::inner_join(ib, by = "horizon_bucket") %>%
        dplyr::filter(is_interior) %>%
        dplyr::arrange(horizon_bucket) %>%
        dplyr::mutate(rank = dplyr::row_number(),
                      offset = r$offset) %>%
        dplyr::select(offset, rank, Q_bar)
    })
    
    if (!nrow(long)) {
      diag_osi[[w]] <- list(
        overall      = tibble::tibble(width = wnum, OSI_abs = NA_real_, OSI_rel = NA_real_,
                                      thr_green = thr_green, thr_orange = thr_orange),
        per_offset   = tibble::tibble(),
        env_long     = tibble::tibble(),
        median_curve = tibble::tibble()
      )
      diag_loso_med[[w]] <- tibble::tibble()
      next
    }
    
    env_long <- long %>%
      dplyr::group_by(rank) %>%
      dplyr::summarise(
        Q_min    = min(Q_bar, na.rm = TRUE),
        Q_max    = max(Q_bar, na.rm = TRUE),
        Q_median = stats::median(Q_bar, na.rm = TRUE),
        Q_mean   = mean(Q_bar, na.rm = TRUE),
        .groups = "drop"
      )
    
    OSI_abs <- max(env_long$Q_max - env_long$Q_min, na.rm = TRUE)
    denom   <- mean(env_long$Q_mean, na.rm = TRUE)
    OSI_rel <- if (is.finite(denom) && denom > 0) OSI_abs/denom else NA_real_
    
    overall_row <- tibble::tibble(
      width = wnum, OSI_abs = OSI_abs, OSI_rel = OSI_rel,
      thr_green = thr_green, thr_orange = thr_orange
    )
    
    per_off <- long %>%
      dplyr::left_join(env_long %>% dplyr::select(rank, Q_median, Q_mean), by = "rank") %>%
      dplyr::mutate(dev = abs(Q_bar - Q_median)) %>%
      dplyr::group_by(offset) %>%
      dplyr::summarise(
        OSI_Linf = max(dev, na.rm = TRUE),
        OSI_q95  = stats::quantile(dev, 0.95, na.rm = TRUE, names = FALSE),
        OSI_RMSE = sqrt(mean(dev^2, na.rm = TRUE)),
        denom    = mean(Q_mean, na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(
        OSIrel_Linf = dplyr::if_else(is.finite(denom) & denom > 0, OSI_Linf/denom, NA_real_),
        OSIrel_q95  = dplyr::if_else(is.finite(denom) & denom > 0, OSI_q95 /denom, NA_real_),
        OSIrel_RMSE = dplyr::if_else(is.finite(denom) & denom > 0, OSI_RMSE/denom, NA_real_)
      )
    
    # LOSO medians across offsets (rank aligned)
    loso_w <- purrr::map_dfr(rr, function(r) {
      bw <- ocp_loso_bucketwise(r, min_frac = min_coverage_frac)
      if (!nrow(bw)) return(tibble::tibble())
      bw %>%
        dplyr::mutate(offset = r$offset) %>%
        dplyr::arrange(horizon_bucket) %>%
        dplyr::mutate(rank = dplyr::row_number()) %>%
        dplyr::select(offset, rank, bucket_influence_mean = infl_mean, bucket_influence_max = infl_max)
    })
    loso_med <- if (nrow(loso_w)) {
      loso_w %>%
        dplyr::group_by(rank) %>%
        dplyr::summarise(
          infl_mean_med = stats::median(bucket_influence_mean, na.rm = TRUE),
          infl_max_med  = stats::median(bucket_influence_max,  na.rm = TRUE),
          .groups = "drop"
        )
    } else tibble::tibble()
    
    diag_osi[[w]]      <- list(overall = overall_row, per_offset = per_off,
                               env_long = env_long,
                               median_curve = env_long %>% dplyr::select(rank, Q_median) %>%
                                 dplyr::rename(Q = Q_median))
    diag_loso_med[[w]] <- loso_med
  }
  
  list(
    grid                = grid,
    results             = results,
    plots               = plots,
    region_plots        = region_plots,
    sampling            = sampling,
    normalization       = normalization,
    edge_correction     = edge_correction,
    min_coverage_frac   = min_coverage_frac,
    diag = list(
      osi          = diag_osi,
      loso_median  = diag_loso_med,
      spatial_ctrl = list(
        selected_buckets = spatial_selected_buckets,
        use_normalized   = spatial_use_normalized,
        bins             = spatial_bins,
        bbox             = spatial_bbox,
        point_color      = spatial_point_color,
        point_alpha      = spatial_point_alpha,
        anim             = spatial_anim,
        anim_bucket_bin  = spatial_anim_bucket_bin
      )
    )
  )
}



#'ocp_median_curve_across_offsets
#'
#' Rank-align interior buckets across offsets for a fixed width and build
#' a median curve (and spread) of Q_bar.
ocp_median_curve_across_offsets <- function(ens, width, min_frac = 1.0) {
  Rs <- purrr::keep(ens$results, ~ .x$width == width)
  if (!length(Rs)) return(tibble::tibble())
  
  series <- purrr::map(Rs, function(r) {
    ib <- interior_mask(r$coverage_tbl, min_frac)
    r$pooled %>%
      dplyr::mutate(H = as.numeric(horizon_bucket)) %>%
      dplyr::inner_join(ib, by = "horizon_bucket") %>%
      dplyr::filter(is_interior) %>%
      dplyr::arrange(H) %>%
      dplyr::transmute(idx = dplyr::row_number(), H = H, Q = Q_bar, offset = r$offset)
  })
  series <- purrr::keep(series, ~ nrow(.x) > 0)
  if (!length(series)) return(tibble::tibble())
  K_min <- min(vapply(series, nrow, integer(1)))
  if (K_min <= 0) return(tibble::tibble())
  
  long <- purrr::map_dfr(series, ~ dplyr::filter(.x, idx <= K_min))
  long %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(
      H_center = stats::median(H, na.rm = TRUE),
      Q_median = stats::median(Q, na.rm = TRUE),
      Q_q25    = stats::quantile(Q, 0.25, na.rm = TRUE, names = FALSE),
      Q_q75    = stats::quantile(Q, 0.75, na.rm = TRUE, names = FALSE),
      Q_min    = min(Q, na.rm = TRUE),
      Q_max    = max(Q, na.rm = TRUE),
      n_offsets= dplyr::n(),
      .groups  = "drop"
    ) %>%
    dplyr::mutate(width = width)
}

#' ocp_plot_median_across_offsets 
#'
#'
#'
ocp_plot_median_across_offsets <- function(ens, width, min_frac = 1.0,
                                           out_path = NULL, print_plot = TRUE) {
  med <- ocp_median_curve_across_offsets(ens, width, min_frac)
  if (nrow(med) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::labs(title = "Median across offsets",
                    subtitle = sprintf("width=%s · no interior/common ranks", width),
                    x = "Pseudo Horizon (median center)", y = "Normalized occupancy") +
      ggplot2::theme_minimal()
    if (!is.null(out_path)) ocp_save_plot(p, out_path, 8, 4.5)
    if (print_plot) print(p)
    return(invisible(p))
  }
  # Reference: first offset points
  r0 <- purrr::detect(ens$results, ~ .x$width == width)
  ref <- NULL
  if (!is.null(r0)) {
    ib0 <- interior_mask(r0$coverage_tbl, min_frac)
    tmp <- r0$pooled %>%
      dplyr::mutate(H = as.numeric(horizon_bucket)) %>%
      dplyr::inner_join(ib0, by = "horizon_bucket") %>%
      dplyr::filter(is_interior) %>%
      dplyr::arrange(H) %>%
      dplyr::transmute(idx = dplyr::row_number(), H, Q = Q_bar)
    ref <- dplyr::filter(tmp, idx <= max(med$idx))
  }
  
  p <- ggplot2::ggplot(med, ggplot2::aes(H_center, Q_median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Q_q25, ymax = Q_q75), alpha = 0.20) +
    ggplot2::geom_line(linewidth = 0.9) +
    { if (!is.null(ref) && nrow(ref)) ggplot2::geom_point(data = ref, ggplot2::aes(H, Q),
                                                          inherit.aes = FALSE, alpha = 0.35, size = 1.4) } +
    ggplot2::labs(title = "Median across offsets (fixed width)",
                  subtitle = sprintf("width=%s · IQR ribbon across offsets", width),
                  x = "Pseudo Horizon (median center across offsets)", y = "Normalized occupancy") +
    ggplot2::theme_minimal()
  if (!is.null(out_path)) ocp_save_plot(p, out_path, 8.5, 4.8)
  if (print_plot) print(p)
  invisible(p)
}

#' fit_and_plot_gam — GAM over pooled means with 1/SE² weights
#'
#' Fits \eqn{Q̄(H) = s(H)} using weights \eqn{w = 1/SE^2} (downweighting
#' noisy buckets), then plots the fitted curve with a pointwise ribbon
#' (±1.96·SE_fit). Returns list(model, plot).
fit_and_plot_gam <- function(pooled, k=10, manifest=list()) {
  dat <- pooled %>%
    mutate(H = as.numeric(horizon_bucket),
           w = ifelse(SE>0, 1/SE^2, 0))
  m <- mgcv::gam(Q_bar ~ s(H, k=k), weights = w, data=dat, method="REML")
  pred <- tibble(H = sort(unique(dat$H)))
  pr   <- predict(m, newdata=pred, se.fit=TRUE)
  pred$fit <- pr$fit; pred$se <- pr$se.fit
  crit <- 1.96
  pred <- pred %>% mutate(lo=fit - crit*se, hi=fit + crit*se)
  
  p <- ggplot() +
    geom_point(data=dat, aes(x=H, y=Q_bar)) +
    geom_errorbar(data=dat, aes(x=H, ymin=CI_lower, ymax=CI_upper), width=0.2, alpha=0.5) +
    geom_line(data=pred, aes(x=H, y=fit)) +
    geom_ribbon(data=pred, aes(x=H, ymin=lo, ymax=hi), alpha=0.2) +
    labs(x="Pseudo Horizon (bucket index)", y="Normalized occupancy",
         title="GAM (weighted by 1/SE²) over pooled means",
         subtitle=paste0(fmt_manifest(manifest), " · n=", nrow(dat))) +
    theme_minimal()
  list(model=m, plot=p)
}

#' plot_loess_eda — LOESS smoother for quick EDA
#'
#' Convenience visualization of pooled per-bucket means with error bars
#' and a LOESS trend. Useful for detecting broad shape without modeling
#' commitments.
plot_loess_eda <- function(pooled, manifest=list()) {
  ggplot(pooled, aes(x=as.numeric(horizon_bucket), y=Q_bar)) +
    geom_point() +
    geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0.2) +
    geom_smooth(method="loess", se=TRUE, span=0.75) +
    labs(x="Pseudo Horizon (bucket index)", y="Normalized occupancy",
         title="LOESS EDA over pooled means",
         subtitle=paste0(fmt_manifest(manifest), " · n=", nrow(pooled))) +
    theme_minimal()
}

#' selection_diagnostics — Source composition over horizons
#'
#' Summarizes counts and mass by (horizon, assign_source, mortuary source),
#' and plots a stacked area of assignment source ("primary"/"meta") over
#' horizons, revealing phase ranges dominated by meta fallback.
selection_diagnostics <- function(harmonized_data) {
  # Table: counts by horizon × assignment source × mortuary source
  tab <- harmonized_data %>%
    group_by(horizon_bucket, assign_source, source) %>%
    summarise(n = n(), w = sum(weight_entry, na.rm=TRUE), .groups="drop")
  
  # Graph: stacked area by source type over horizons
  p <- tab %>%
    group_by(horizon_bucket, assign_source) %>%
    summarise(w=sum(w), .groups="drop") %>%
    ggplot(aes(x=as.numeric(horizon_bucket), y=w, fill=assign_source)) +
    geom_area(alpha=0.8) +
    labs(x="Pseudo Horizon (bucket index)", y="Total weight",
         title="Assignment source mix over Pseudo Horizons") +
    theme_minimal()
  
  list(table = tab, plot = p)
}

# ===============================
## Model Testing Functions
# Provide metrics for all data modelling operations, including a fixed uniform
# dataset - as well as metrics for the real-world dataset. Saves as .csv
# and outputs relevant plots for inspection.
# ===============================

#' interior_mask — Mark interior buckets by coverage threshold
#'
#' Returns a mask of buckets whose coverage fraction \eqn{c(b) ≥ min_frac}.
#' Commonly used to remove edge-affected buckets from summaries/metrics.
interior_mask <- function(coverage_tbl, min_frac = 1.0) {
  stopifnot(all(c("horizon_bucket","cover_frac") %in% names(coverage_tbl)))
  coverage_tbl %>% transmute(horizon_bucket, is_interior = cover_frac >= min_frac)
}

#' ocp_make_uniform_fixture — Synthetic calibration dataset
#'
#' Generates an idealized dataset with evenly tiled phases and uniform
#' per-Phase counts so the true occupancy curve is flat under perfect
#' conditions. Useful for validation of offset/edge behavior and QC.
#' Produces (harm_tab, meta, indiv, agg).
ocp_make_uniform_fixture <- function(
    n_sites           = 3L,
    n_phases          = 10L,        # phases tile the axis exactly
    per_site_total    = 100L,       # total burials per site (will be snapped to a multiple of n_phases)
    phase_length_phu  = 10L,        # PHU per phase; use 10 to align with width=10 tests
    start_phu         = 0L,         # first phase starts here (allows testing negative shifts if desired)
    system_name       = "sys_A",
    cemetery_size_mult= 1.20,       # cemetery_size = per_site_total * multiplier (rounded)
    sex_value         = "unknown",  # keep constant to avoid noise
    age_value         = "unknown",
    seed              = ocp_get_global_seed()  # used only for reproducible coords/labels (not counts)
) {
  stopifnot(n_sites >= 1, n_phases >= 1, per_site_total >= 0, phase_length_phu >= 1)
  
  # Snap totals so each phase gets the same integer count
  per_site_total <- as.integer(round(per_site_total / n_phases) * n_phases)
  per_phase_amount <- as.integer(per_site_total / n_phases)
  
  ## Harmonization (single system)
  
  harm_tab <- tibble::tibble(
    phase_id       = seq_len(n_phases),
    system_name    = system_name,
    horizon_start  = start_phu + (phase_id - 1L) * phase_length_phu,
    horizon_end    = start_phu +  phase_id      * phase_length_phu,
    fade_in_start  = NA_integer_,
    fade_out_end   = NA_integer_
  )
  
  ## Meta (n_sites, identical span)
  
  withr::with_seed(seed, {
    regions   <- c("North","Central","South","East","West")
    countries <- c("X","Y","Z")
    meta <- tibble::tibble(
      coord_y = 50 + stats::rnorm(n_sites, 0, 0.2),
      coord_x = 10 + stats::rnorm(n_sites, 0, 0.2),
      site_id            = seq_len(n_sites),
      site_name          = paste0("Site_", sprintf("%02d", site_id)),
      cemetery_lit1      = NA_character_,
      cemetery_lit2      = NA_character_,
      cemetery_lit3      = NA_character_,
      cemetery_lit4      = NA_character_,
      cemetery_country   = rep_len(countries, n_sites),
      cemetery_region    = rep_len(regions,   n_sites),
      cemetery_admin     = rep_len(c("AdminA","AdminB"), n_sites),
      cemetery_start_ia  = 1L,
      cemetery_start_pre = 1L,
      cemetery_end_ia    = n_phases,
      cemetery_end_post  = n_phases,
      cemetery_size      = as.integer(round(per_site_total * cemetery_size_mult)),
      dig_date           = rep_len(c("1990-1995","1996-2000","2001-2005","2006-2010"), n_sites)
    )
  
  ## Aggregated mortuary (strictly uniform)

  agg <- tidyr::crossing(
    site_id  = meta$site_id,
    phase_id = harm_tab$phase_id
  ) |>
    dplyr::arrange(site_id, phase_id) |>
    dplyr::mutate(
      ID         = dplyr::row_number(),
      sex_gender = sex_value,
      age        = age_value,
      amount     = per_phase_amount
    ) |>
    dplyr::select(ID, site_id, sex_gender, age, amount, phase_id)
  
  ## Individual mortuary (expanded uniform)

  indiv <- agg |>
    tidyr::uncount(weights = amount) |>
    dplyr::arrange(site_id, phase_id) |>
    dplyr::group_by(site_id) |>
    dplyr::mutate(seq = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      ID         = dplyr::row_number(),
      burial_id  = paste0("S", site_id, "I", sprintf("%04d", seq)),
      sex_gender = sex_gender,
      age        = age,
      phase_id   = phase_id,
      site_id    = site_id
    )
  })
  list(harm_tab = harm_tab, meta = meta, indiv = indiv, agg = agg)
}

# ------------------------
# Metrics 
# ------------------------

#' metric_flatness — Is the interior curve flat (calibration)?
#'
#' On interior buckets, computes mean, sd, CV, and the slope/R² of a
#' linear fit \eqn{Q̄ ~ H}. For a uniform fixture, sd≈0 and slope≈0.
#' Flags edge or allocation biases if not near-flat.
metric_flatness <- function(res, min_frac = 1.0) {
  ib <- interior_mask(res$coverage_tbl, min_frac)
  pooled <- res$pooled %>%
    dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
    dplyr::inner_join(ib, by="horizon_bucket") %>%
    dplyr::filter(is_interior)
  
  m <- mean(pooled$Q_bar, na.rm = TRUE)
  s <- stats::sd(pooled$Q_bar,   na.rm = TRUE)
  
  if (!is.finite(s) || s == 0) {
    slope <- 0
    r2    <- NA_real_
    note  <- "Constant series on interior (sd=0); slope≈0; R2 is NA by design."
  } else {
    fit   <- stats::lm(Q_bar ~ as.numeric(horizon_bucket), data = pooled)
    slope <- unname(stats::coef(fit)[2])
    r2    <- summary(fit)$r.squared
    note  <- NA_character_
  }
  
  tibble::tibble(
    width = res$width, offset = res$offset,
    n_buckets = nrow(pooled),
    mean_Q = m, sd_Q = s, cv_Q = ifelse(m > 0, s/m, NA_real_),
    slope_lm = slope, r2_lm = r2, note = note
  )
}

#' metric_offset_sensitivity — Envelope across offsets
#'
#' For a fixed width, collects interior buckets per offset and computes
#' the envelope: per-bucket min/max of Q̄ across offsets and then
#' \eqn{\mathrm{OSI}_{abs} = \max_H (Q_{max}-Q_{min})} and
#' \eqn{\mathrm{OSI}_{rel} = \mathrm{OSI}_{abs}/\text{mean}(Q̄)}.
#' Requires a common interior set; otherwise returns NA with message.
metric_offset_sensitivity <- function(ens, width, min_frac = 1.0, grid_step = NULL) {
  # collect (width==width) results
  res_list <- Filter(function(z) isTRUE(z$width == width), ens$results)
  if (length(res_list) < 2) {
    msg <- sprintf("[metric_offset_sensitivity] Need >=2 offsets for width=%s.", width)
    message(msg)
    return(tibble::tibble(width = width, OSI_abs = NA_real_, OSI_rel = NA_real_,
                          issue = TRUE, msg = msg))
  }
  w <- width
  # interior continuous ranges per offset
  ranges <- lapply(res_list, function(r) {
    ib <- r$coverage_tbl %>% dplyr::filter(cover_frac >= min_frac)
    if (nrow(ib) == 0) return(NULL)
    c(min_start = min(as.numeric(ib$horizon_bucket), na.rm = TRUE),
      max_end   = max(as.numeric(ib$horizon_bucket), na.rm = TRUE) + w)
  })
  ranges <- Filter(Negate(is.null), ranges)
  if (!length(ranges)) {
    msg <- sprintf("[metric_offset_sensitivity] No interior buckets for width=%s after filtering.", width)
    message(msg)
    return(tibble::tibble(width = width, OSI_abs = NA_real_, OSI_rel = NA_real_,
                          issue = TRUE, msg = msg))
  }
  ranges_mat <- do.call(rbind, ranges)
  x_min <- max(ranges_mat[, "min_start"])
  x_max <- min(ranges_mat[, "max_end"])
  if (!(is.finite(x_min) && is.finite(x_max)) || x_max <= x_min) {
    msg <- sprintf("[metric_offset_sensitivity] No continuous overlap across offsets for width=%s.", width)
    message(msg)
    return(tibble::tibble(width = width, OSI_abs = NA_real_, OSI_rel = NA_real_,
                          issue = TRUE, msg = msg))
  }
  # grid step: default to about 5 samples per bucket (>=1)
  if (is.null(grid_step)) grid_step <- max(1, floor(w / 5))
  x_grid <- seq(x_min, x_max - 1e-9, by = grid_step)
  
  # evaluate each offset's pooled step function on the common grid
  eval_on_grid <- function(r, xg) {
    o <- r$offset; w <- r$width
    # bucket start that contains x: floor((x - o)/w)*w + o
    starts <- floor((xg - o) / w) * w + o
    df <- r$pooled %>% dplyr::mutate(hb = as.numeric(horizon_bucket))
    idx <- match(starts, df$hb)
    df$Q_bar[idx]
  }
  Q_mat <- do.call(cbind, lapply(res_list, eval_on_grid, xg = x_grid))
  # keep rows where all offsets have values
  keep <- rowSums(!is.na(Q_mat)) == ncol(Q_mat)
  if (!any(keep)) {
    msg <- sprintf("[metric_offset_sensitivity] Overlap exists but pooled series are missing on the common grid (width=%s).", width)
    message(msg)
    return(tibble::tibble(width = width, OSI_abs = NA_real_, OSI_rel = NA_real_,
                          issue = TRUE, msg = msg))
  }
  Q_sub <- Q_mat[keep, , drop = FALSE]
  row_min  <- apply(Q_sub, 1, min, na.rm = TRUE)
  row_max  <- apply(Q_sub, 1, max, na.rm = TRUE)
  row_mean <- apply(Q_sub, 1, mean, na.rm = TRUE)
  
  OSI_abs <- max(row_max - row_min, na.rm = TRUE)
  OSI_rel <- OSI_abs / mean(row_mean, na.rm = TRUE)
  
  tibble::tibble(width = width, OSI_abs = OSI_abs, OSI_rel = OSI_rel,
                 issue = FALSE, msg = NA_character_)
}

#' metric_osi_reference — Pairwise Offset Sensitivity to a reference offset 
#' 
#' Compares each offset's pooled mean to a chosen reference on their
#' SHARED interior buckets only (no "common across ALL" requirement).
#' Returns OSI_abs = max|Q_ref - Q_off| and OSI_rel = OSI_abs / mean(Q_ref).
metric_osi_reference <- function(ens, width, ref_offset = NULL, min_frac = 1.0) {
  # pick reference (first result with matching width if not given)
  cand <- purrr::keep(ens$results, ~ .x$width == width)
  if (!length(cand)) {
    return(tibble::tibble(width = width, offset = numeric(0), OSI_abs = numeric(0), OSI_rel = numeric(0)))
  }
  if (is.null(ref_offset)) ref_offset <- cand[[1]]$offset
  ref <- purrr::detect(cand, ~ .x$offset == ref_offset)
  if (is.null(ref)) ref <- cand[[1]]
  
  ref_ib <- interior_mask(ref$coverage_tbl, min_frac)
  ref_curve <- ref$pooled %>%
    dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
    dplyr::inner_join(ref_ib, by = "horizon_bucket") %>%
    dplyr::filter(is_interior) %>%
    dplyr::select(horizon_bucket, Q_ref = Q_bar)
  
  out <- purrr::map_dfr(cand, function(r) {
    off_ib <- interior_mask(r$coverage_tbl, min_frac)
    off_curve <- r$pooled %>%
      dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
      dplyr::inner_join(off_ib, by = "horizon_bucket") %>%
      dplyr::filter(is_interior) %>%
      dplyr::select(horizon_bucket, Q_off = Q_bar)
    
    shared <- dplyr::inner_join(ref_curve, off_curve, by = "horizon_bucket")
    if (nrow(shared) == 0) {
      tibble::tibble(width = r$width, offset = r$offset, OSI_abs = NA_real_, OSI_rel = NA_real_)
    } else {
      d <- abs(shared$Q_ref - shared$Q_off)
      osi_abs <- max(d, na.rm = TRUE)
      mref    <- mean(shared$Q_ref, na.rm = TRUE)
      tibble::tibble(
        width = r$width, offset = r$offset,
        OSI_abs = osi_abs,
        OSI_rel = if (is.finite(mref) && mref > 0) osi_abs / mref else NA_real_
      )
    }
  })
  out
}

# Rank-align curves across OFFSETS for a single WIDTH and compute:
# (i) overall OSI (max envelope spread across offsets), and
# (ii) per-offset OSI vs the rank-median curve (diagnostic).
#
# Let Q_{o,r} be pooled mean at rank r for offset o (ranks are by ascending interior bucket).
# Define spread_r = max_o Q_{o,r} - min_o Q_{o,r}.
# Overall OSI_abs = max_r spread_r; OSI_rel = OSI_abs / mean_r(mean_o Q_{o,r}).
# Per-offset OSI_abs(o) = max_r |Q_{o,r} - median_o(Q_{o,r})|; OSI_rel(o) uses SAME denominator.
ocp_compute_osi_ranked <- function(ens,
                                   width,
                                   min_frac   = 1.0,
                                   thr_green  = 0.05,
                                   thr_orange = 0.10) {
  stopifnot(is.list(ens), "results" %in% names(ens))
  # Collect rank-aligned curves for this width
  curves <- lapply(ens$results, function(r) {
    if (r$width != width) return(NULL)
    ib <- interior_mask(r$coverage_tbl, min_frac)
    df <- r$pooled %>%
      dplyr::mutate(H = as.numeric(horizon_bucket)) %>%
      dplyr::inner_join(ib, by = "horizon_bucket") %>%
      dplyr::filter(is_interior) %>%
      dplyr::arrange(H) %>%
      dplyr::transmute(offset = r$offset, rank = dplyr::row_number(), Q = Q_bar)
    if (nrow(df) >= 2) df else NULL
  })
  curves <- dplyr::bind_rows(curves)
  if (nrow(curves) == 0) {
    return(list(
      overall    = tibble::tibble(width = width, OSI_abs = NA_real_, OSI_rel = NA_real_,
                                  n_offsets_used = 0L, min_k = 0L, note = "no interior data"),
      per_offset = tibble::tibble(width = width, offset = integer(), OSI_abs = numeric(),
                                  OSI_rel = numeric(), status = character())
    ))
  }
  # Use only ranks common to all offsets (min_k)
  k_by_off <- curves %>% dplyr::group_by(offset) %>% dplyr::summarise(k = max(rank), .groups="drop")
  min_k    <- min(k_by_off$k)
  if (!is.finite(min_k) || min_k < 2) {
    return(list(
      overall    = tibble::tibble(width = width, OSI_abs = NA_real_, OSI_rel = NA_real_,
                                  n_offsets_used = dplyr::n_distinct(curves$offset),
                                  min_k = min_k, note = "min_k<2"),
      per_offset = tibble::tibble(width = width, offset = integer(), OSI_abs = numeric(),
                                  OSI_rel = numeric(), status = character())
    ))
  }
  C <- curves %>% dplyr::filter(rank <= min_k)
  
  # Overall envelope spread across offsets per rank
  spread <- C %>% dplyr::group_by(rank) %>%
    dplyr::summarise(minQ = min(Q, na.rm=TRUE),
                     maxQ = max(Q, na.rm=TRUE),
                     meanQ = mean(Q, na.rm=TRUE), .groups="drop")
  OSI_abs <- max(spread$maxQ - spread$minQ, na.rm=TRUE)
  denom   <- mean(spread$meanQ, na.rm=TRUE)
  OSI_rel <- ifelse(is.finite(denom) && denom > 0, OSI_abs / denom, NA_real_)
  
  # Per-offset deviation from the rank-median curve
  med_rank <- C %>% dplyr::group_by(rank) %>%
    dplyr::summarise(medQ = stats::median(Q, na.rm=TRUE), .groups="drop")
  per_off <- C %>% dplyr::left_join(med_rank, by = "rank") %>%
    dplyr::group_by(offset) %>%
    dplyr::summarise(
      OSI_abs = max(abs(Q - medQ), na.rm=TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      width  = width,
      OSI_rel = ifelse(is.finite(denom) && denom > 0, OSI_abs / denom, NA_real_),
      status = dplyr::case_when(
        OSI_rel <= thr_green  ~ "green",
        OSI_rel <= thr_orange ~ "orange",
        TRUE                  ~ "red"
      )
    ) %>% dplyr::select(width, offset, OSI_abs, OSI_rel, status) %>%
    dplyr::arrange(offset)
  
  list(
    overall    = tibble::tibble(
      width = width, OSI_abs = OSI_abs, OSI_rel = OSI_rel,
      n_offsets_used = dplyr::n_distinct(C$offset), min_k = min_k, note = NA_character_
    ),
    per_offset = per_off
  )
}


# Simple OSI status dot
# Plot OSI_rel per OFFSET for a given width, traffic-light colored,
# with horizontal threshold lines. Accepts the $per_offset tibble from ocp_compute_osi_ranked().
ocp_plot_osi_status <- function(osi_per_offset_df,
                                metric = c("OSIrel_q95","OSIrel_RMSE","OSIrel_Linf"),
                                thr_green, thr_orange,
                                out_path = NULL) {
  metric <- match.arg(metric)
  df <- osi_per_offset_df %>%
    dplyr::mutate(
      status = dplyr::case_when(
        .data[[metric]] <= thr_green ~ "green",
        .data[[metric]] <= thr_orange ~ "orange",
        TRUE ~ "red"
      )
    )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(offset), y = .data[[metric]], fill = status)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = thr_green, linetype = 2) +
    ggplot2::geom_hline(yintercept = thr_orange, linetype = 2) +
    ggplot2::scale_fill_manual(values = c(green = "#2ca02c", orange = "#ff7f0e", red = "#d62728")) +
    ggplot2::labs(x = "Offset", y = metric, title = "Offset sensitivity per offset (rank-aligned)") +
    ggplot2::theme_minimal()
  if (!is.null(out_path)) ocp_save_plot(p, out_path, width = 7, height = 4)
  p
}


#' metric_replicate_stability — Between vs within variance
#'
#' On interior buckets, compares replicate-to-replicate variability
#' (between variance of mean entries) to within-replicate variance
#' (first replicate’s var_hat). Reports their means and the ratio
#' \eqn{B/W}. Large ratios indicate unstable stochastic allocations.
metric_replicate_stability <- function(res, min_frac = 1.0) {
  ib <- interior_mask(res$coverage_tbl, min_frac)
  M <- length(res$reps)
  
  # Collect replicate series (use corrected per_bucket if present)
  reps_long <- purrr::imap_dfr(res$reps, function(rep, idx) {
    pb <- if ("per_bucket" %in% names(rep)) rep$per_bucket else rep$per_bucket_raw
    pb %>%
      dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
      dplyr::inner_join(ib, by = "horizon_bucket") %>%
      dplyr::filter(is_interior) %>%
      dplyr::transmute(horizon_bucket, mean_entries, rep = idx)
  })
  
  if (M < 2) {
    # Deterministic: define between variance as 0
    return(tibble::tibble(
      width = res$width, offset = res$offset,
      mean_var_between = 0, mean_var_within = {
        # take within from any replicate on interior
        pb1 <- if ("per_bucket" %in% names(res$reps[[1]])) res$reps[[1]]$per_bucket else res$reps[[1]]$per_bucket_raw
        jj  <- pb1 %>%
          dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
          dplyr::inner_join(ib, by = "horizon_bucket") %>%
          dplyr::filter(is_interior)
        mean(jj$var_hat, na.rm = TRUE)
      },
      ratio_B_over_W = 0, M_reps = M
    ))
  }
  
  # Common interior across replicates
  ok_bks <- reps_long %>% dplyr::count(horizon_bucket) %>%
    dplyr::filter(n == M) %>% dplyr::pull(horizon_bucket)
  reps_common <- reps_long %>% dplyr::filter(horizon_bucket %in% ok_bks)
  
  # Between variance per bucket (only if >=2 reps by construction)
  between <- reps_common %>% dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(var_between = stats::var(mean_entries), .groups="drop")
  
  # Within variance from the first replicate, restricted to common interior
  pb1 <- if ("per_bucket" %in% names(res$reps[[1]])) res$reps[[1]]$per_bucket else res$reps[[1]]$per_bucket_raw
  w_tbl <- pb1 %>%
    dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
    dplyr::filter(horizon_bucket %in% ok_bks) %>%
    dplyr::select(horizon_bucket, var_hat)
  
  jj <- dplyr::left_join(between, w_tbl, by = "horizon_bucket")
  tibble::tibble(
    width = res$width, offset = res$offset,
    mean_var_between = mean(jj$var_between, na.rm = TRUE),
    mean_var_within  = mean(jj$var_hat,     na.rm = TRUE),
    ratio_B_over_W   = mean(jj$var_between, na.rm = TRUE) / mean(jj$var_hat, na.rm = TRUE),
    M_reps = M
  )
}

#' metric_meta_dependence — How much meta fallback drives mass?
#'
#' On interior buckets, computes the mean and max fraction of total mass
#' contributed by meta fallback (\eqn{\sum \text{meta} / \sum \text{total}}).
#' High values indicate heavy reliance on meta intervals.
metric_meta_dependence <- function(res, min_frac = 1.0) {
  sb <- res$reps[[1]]$diagnostics$site_bucket
  if (!all(c("primary_share","meta_share") %in% names(sb))) {
    return(tibble(width = res$width, offset = res$offset,
                  mean_meta_share = NA_real_, max_meta_share = NA_real_))
  }
  ib <- interior_mask(res$coverage_tbl, min_frac)
  tmp <- sb %>% inner_join(ib, by="horizon_bucket") %>% filter(is_interior) %>%
    group_by(horizon_bucket) %>%
    summarise(total = sum(primary_share + meta_share, na.rm = TRUE),
              meta  = sum(meta_share, na.rm = TRUE), .groups="drop") %>%
    mutate(meta_prop = ifelse(total > 0, meta/total, NA_real_))
  tibble(width = res$width, offset = res$offset,
         mean_meta_share = mean(tmp$meta_prop, na.rm=TRUE),
         max_meta_share  = max(tmp$meta_prop,  na.rm=TRUE))
}

#' metric_loso_influence — Leave-one-site-out influence
#'
#' For each site i and bucket H, computes how much the bucket mean would
#' change if site i were removed:
#' \eqn{\mu = \text{mean}(Q)}, \eqn{n = \#\text{sites}},
#' \eqn{\mu_{-i} = (n\mu - Q_i)/\max(n-1,1)},
#' \eqn{\mathrm{infl}_i = |\mu_{-i} - \mu|}.
#' Reports max and mean influence across sites; large values indicate
#' undue leverage by individual sites.
metric_loso_influence <- function(res, min_frac = 1.0) {
  sb <- res$reps[[1]]$diagnostics$site_bucket
  ib <- interior_mask(res$coverage_tbl, min_frac)
  sb_int <- sb %>% inner_join(ib, by="horizon_bucket") %>% filter(is_interior)
  mu <- sb_int %>% group_by(horizon_bucket) %>%
    summarise(mu = mean(normalized_entries), n = n(), .groups="drop")
  jj <- sb_int %>% left_join(mu, by="horizon_bucket") %>%
    mutate(mean_minus_i = (n*mu - normalized_entries)/pmax(n-1, 1),
           infl = abs(mean_minus_i - mu))
  tibble(width = res$width, offset = res$offset,
         max_influence = max(jj$infl, na.rm=TRUE),
         mean_influence = mean(jj$infl, na.rm=TRUE))
}

#' ocp_compute_loso — Leave-One-Site-Out (LOSO) diagnostics
#'
#'Computes:
#'  - bucketwise LOSO (max |mu - mu_{-i}| across sites)
#'  - summary (max/mean LOSO abs; relative to mean pooled level and median SE)
#'  - per-site cumulative influence (sum across buckets)
#' Optionally bootstraps site-level control limits (off by default).
ocp_compute_loso <- function(res,
                             min_frac       = 1.0,
                             do_bootstrap   = FALSE,
                             B_boot         = 200L,
                             boot_seed      = ocp_get_global_seed()) {
  stopifnot(is.list(res), "reps" %in% names(res), length(res$reps) >= 1)
  sb  <- res$reps[[1]]$diagnostics$site_bucket
  ib  <- interior_mask(res$coverage_tbl, min_frac)
  sbI <- sb %>%
    dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
    dplyr::inner_join(ib, by = "horizon_bucket") %>%
    dplyr::filter(is_interior)
  
  if (nrow(sbI) == 0) {
    return(list(
      summary = tibble::tibble(
        width = res$width, offset = res$offset,
        max_infl_abs = NA_real_, mean_infl_abs = NA_real_,
        rel_to_mean  = NA_real_, rel_to_medSE  = NA_real_
      ),
      bucketwise = tibble::tibble(horizon_bucket = numeric(0), LOSO_abs = numeric(0), SE = numeric(0)),
      top_sites  = tibble::tibble(site_id = character(0), cum_infl = numeric(0), max_infl = numeric(0), n_buckets = integer(0))
    ))
  }
  
  # Bucketwise mean across sites (mu) and n
  mu_tbl <- sbI %>%
    dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(
      mu = mean(normalized_entries, na.rm = TRUE),
      n  = dplyr::n(),
      .groups = "drop"
    )
  
  # For each (site, bucket): mu_{-i} = (n*mu - x_i) / (n-1); influence = |mu_{-i} - mu|
  jj <- sbI %>%
    dplyr::left_join(mu_tbl, by = "horizon_bucket") %>%
    dplyr::mutate(
      mu_minus_i = (n * mu - normalized_entries) / pmax(n - 1, 1),
      infl       = abs(mu_minus_i - mu)
    )
  
  # Bucketwise LOSO = max influence across sites (per bucket)
  loso_bucket <- jj %>%
    dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(LOSO_abs = max(infl, na.rm = TRUE), .groups = "drop")
  
  # Join pooled SE for reference comparisons (if available)
  pooled_int <- res$pooled %>%
    dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
    dplyr::inner_join(ib, by = "horizon_bucket") %>%
    dplyr::filter(is_interior)
  loso_bucket <- loso_bucket %>%
    dplyr::left_join(pooled_int %>% dplyr::select(horizon_bucket, SE), by = "horizon_bucket")
  
  # Per-site cumulative influence (sum across buckets)
  infl_sites <- jj %>%
    dplyr::group_by(site_id) %>%
    dplyr::summarise(
      cum_infl  = sum(infl, na.rm = TRUE),
      max_infl  = max(infl, na.rm = TRUE),
      n_buckets = dplyr::n_distinct(horizon_bucket),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(cum_infl))
  
  # Absolute & relative summaries
  max_abs  <- max(loso_bucket$LOSO_abs, na.rm = TRUE)
  mean_abs <- mean(loso_bucket$LOSO_abs, na.rm = TRUE)
  mean_Q   <- mean(pooled_int$Q_bar, na.rm = TRUE)
  med_SE   <- stats::median(pooled_int$SE, na.rm = TRUE)
  
  summ <- tibble::tibble(
    width = res$width,
    offset = res$offset,
    max_infl_abs = max_abs,
    mean_infl_abs = mean_abs,
    rel_to_mean  = if (is.finite(mean_Q) && mean_Q > 0) max_abs / mean_Q else NA_real_,
    rel_to_medSE = if (is.finite(med_SE) && med_SE > 0) max_abs / med_SE else NA_real_
  )
  
  # Optional bootstrap control limit (bucketwise 95th of |mu_boot - mu|)
  if (isTRUE(do_bootstrap)) {
    if (!is.null(boot_seed)) set.seed(boot_seed)
    Hvals <- sort(unique(sbI$horizon_bucket))
    boot_max <- numeric(length(Hvals))
    names(boot_max) <- as.character(Hvals)
    
    # Precompute site lists per bucket
    site_by_H <- sbI %>% dplyr::group_by(horizon_bucket) %>%
      dplyr::summarise(sites = list(split(normalized_entries, site_id)), .groups = "drop")
    
    for (h in Hvals) {
      lst <- site_by_H$sites[[ match(h, site_by_H$horizon_bucket) ]]
      if (is.null(lst) || length(lst) < 2) next
      mu0 <- mu_tbl$mu[match(h, mu_tbl$horizon_bucket)]
      # Bootstrap mean across sites (resample sites with replacement)
      boots <- replicate(B_boot, {
        sids <- sample(seq_along(lst), replace = TRUE)
        vals <- vapply(sids, function(i) mean(lst[[i]], na.rm = TRUE), numeric(1))
        mean(vals, na.rm = TRUE)
      })
      boot_max[as.character(h)] <- stats::quantile(abs(boots - mu0), probs = 0.95, na.rm = TRUE)
    }
    loso_bucket$boot95 <- as.numeric(boot_max[as.character(loso_bucket$horizon_bucket)])
  }
  
  list(summary = summ, bucketwise = loso_bucket, top_sites = infl_sites)
}

#' title - short summary
#' 
#' explanation
#'
#'

# ------------------------
# Plots
# ------------------------

#' ocp_plot_loso_timeseries
#' LOSO time series (bucketwise)
#' Draws LOSO_abs across interior buckets; optional reference bands at 0.25×SE_ref, 0.5×SE_ref
ocp_plot_loso_timeseries <- function(res, loso_bucketwise,
                                     se_ref = NULL,
                                     out_path = NULL, print_plot = TRUE) {
  df <- loso_bucketwise %>%
    dplyr::mutate(H = as.numeric(horizon_bucket))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = H, y = LOSO_abs)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::labs(
      title = "LOSO influence (bucketwise max |µ - µ_{-i}|)",
      subtitle = sprintf("width=%s · offset=%s", res$width, res$offset),
      x = "Pseudo Horizon (bucket)", y = "LOSO (absolute)"
    ) +
    ggplot2::theme_minimal()
  
  if (is.finite(se_ref) && se_ref > 0) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0.25 * se_ref, linetype = 2, alpha = 0.5) +
      ggplot2::geom_hline(yintercept = 0.50 * se_ref, linetype = 2, alpha = 0.5)
  }
  if (!is.null(out_path)) ocp_save_plot(p, out_path, width = 8, height = 4.5)
  if (isTRUE(print_plot)) print(p)
  invisible(p)
}

#' ocp_plot_loso_summary 
#' 
#' LOSO summary bars (+ optional bootstrap control line)
ocp_plot_loso_summary <- function(res, loso_summary, boot_max95 = NA_real_,
                                  out_path = NULL, print_plot = TRUE) {
  df <- tidyr::pivot_longer(
    loso_summary %>% dplyr::select(max_infl_abs, mean_infl_abs),
    cols = everything(), names_to = "metric", values_to = "value"
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = metric, y = value, fill = metric)) +
    ggplot2::geom_col(width = 0.55, show.legend = FALSE) +
    ggplot2::labs(
      title = "LOSO summary",
      subtitle = sprintf("width=%s · offset=%s", res$width, res$offset),
      x = NULL, y = "Absolute influence"
    ) +
    ggplot2::theme_minimal()
  if (is.finite(boot_max95) && boot_max95 > 0) {
    p <- p + ggplot2::geom_hline(yintercept = boot_max95, linetype = 2)
  }
  if (!is.null(out_path)) ocp_save_plot(p, out_path, width = 5.2, height = 4.2)
  if (isTRUE(print_plot)) print(p)
  invisible(p)
}

# Bucketwise LOSO influence for a single (width, offset) result.
# Returns horizon_bucket, infl_mean, infl_max plus per-bucket n.
ocp_loso_bucketwise <- function(res, min_frac = 1.0) {
  sb <- res$reps[[1]]$diagnostics$site_bucket
  ib <- interior_mask(res$coverage_tbl, min_frac)
  sb_int <- sb %>% dplyr::inner_join(ib, by="horizon_bucket") %>% dplyr::filter(is_interior)
  if (!nrow(sb_int)) return(tibble::tibble())
  mu <- sb_int %>%
    dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(mu = mean(normalized_entries, na.rm = TRUE),
                     n  = dplyr::n(),
                     .groups="drop")
  jj <- sb_int %>% dplyr::left_join(mu, by = "horizon_bucket") %>%
    dplyr::mutate(mean_minus_i = (n*mu - normalized_entries)/pmax(n-1, 1),
                  infl = abs(mean_minus_i - mu))
  jj %>% dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(infl_mean = mean(infl, na.rm = TRUE),
                     infl_max  = max(infl, na.rm = TRUE),
                     n = dplyr::n(),
                     .groups="drop")
}

# Plot for first width+offset; traffic-light on points
ocp_plot_loso_bucketwise <- function(res, min_frac = 1.0,
                                     green_thr = 1e-3, red_thr = 5e-3,
                                     out_path = NULL, print_plot = TRUE) {
  df <- ocp_loso_bucketwise(res, min_frac)
  if (!nrow(df)) {
    p <- ggplot2::ggplot() + ggplot2::labs(title="LOSO bucketwise", subtitle="No interior", x="H", y="Influence") + ggplot2::theme_minimal()
    if (!is.null(out_path)) ocp_save_plot(p, out_path, 7.5, 4.5)
    if (print_plot) print(p); return(invisible(p))
  }
  df <- df %>%
    dplyr::mutate(H = as.numeric(horizon_bucket),
                  status = ocp_status(infl_mean, green_thr, red_thr),
                  col = dplyr::recode(status, green="#2ca02c", orange="#ff7f0e", red="#d62728", .default="#7f7f7f"))
  p <- ggplot2::ggplot(df, ggplot2::aes(H, infl_mean)) +
    ggplot2::geom_line(alpha = 0.4) +
    ggplot2::geom_point(ggplot2::aes(colour = status)) +
    ggplot2::scale_colour_manual(values = c(green="#2ca02c", orange="#ff7f0e", red="#d62728", grey="#7f7f7f"), guide = "none") +
    ggplot2::geom_hline(yintercept = c(green_thr, red_thr), linetype = 2, alpha = 0.4) +
    ggplot2::labs(title = "LOSO influence by bucket (first width+offset)",
                  subtitle = sprintf("width=%s · offset=%s", res$width, res$offset),
                  x = "Pseudo Horizon (bucket)", y = "Mean LOSO influence") +
    ggplot2::theme_minimal()
  if (!is.null(out_path)) ocp_save_plot(p, out_path, 8.2, 4.8)
  if (print_plot) print(p)
  invisible(p)
}

# Median LOSO across offsets at fixed width (rank aligned)
ocp_median_loso_across_offsets <- function(ens, width, min_frac = 1.0) {
  Rs <- purrr::keep(ens$results, ~ .x$width == width)
  if (!length(Rs)) return(tibble::tibble())
  series <- purrr::map(Rs, function(r) {
    df <- ocp_loso_bucketwise(r, min_frac)
    if (!nrow(df)) return(df)
    df %>% dplyr::mutate(H = as.numeric(horizon_bucket)) %>%
      dplyr::arrange(H) %>%
      dplyr::transmute(idx = dplyr::row_number(), H = H, L = infl_mean, offset = r$offset)
  })
  series <- purrr::keep(series, ~ nrow(.x) > 0)
  if (!length(series)) return(tibble::tibble())
  K <- min(vapply(series, nrow, integer(1)))
  long <- purrr::map_dfr(series, ~ dplyr::filter(.x, idx <= K))
  long %>% dplyr::group_by(idx) %>%
    dplyr::summarise(H_center = stats::median(H, na.rm = TRUE),
                     L_median = stats::median(L, na.rm = TRUE),
                     L_q25    = stats::quantile(L, 0.25, na.rm = TRUE, names = FALSE),
                     L_q75    = stats::quantile(L, 0.75, na.rm = TRUE, names = FALSE),
                     n_offsets = dplyr::n(),
                     .groups = "drop") %>%
    dplyr::mutate(width = width)
}

ocp_plot_median_loso_across_offsets <- function(ens, width, min_frac = 1.0,
                                                green_thr = 1e-3, red_thr = 5e-3,
                                                out_path = NULL, print_plot = TRUE) {
  med <- ocp_median_loso_across_offsets(ens, width, min_frac)
  if (!nrow(med)) {
    p <- ggplot2::ggplot() + ggplot2::labs(title="Median LOSO across offsets", subtitle="No data", x="H", y="Influence") + ggplot2::theme_minimal()
    if (!is.null(out_path)) ocp_save_plot(p, out_path, 8, 4.5)
    if (print_plot) print(p); return(invisible(p))
  }
  p <- ggplot2::ggplot(med, ggplot2::aes(H_center, L_median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = L_q25, ymax = L_q75), alpha = 0.20) +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(colour = ocp_status(L_median, green_thr, red_thr))) +
    ggplot2::scale_colour_manual(values = c(green="#2ca02c", orange="#ff7f0e", red="#d62728", grey="#7f7f7f"), guide = "none") +
    ggplot2::geom_hline(yintercept = c(green_thr, red_thr), linetype = 2, alpha = 0.4) +
    ggplot2::labs(title = "Median LOSO influence across offsets (fixed width)",
                  subtitle = sprintf("width=%s · IQR ribbon across offsets", width),
                  x = "Pseudo Horizon (median center)", y = "Mean LOSO influence") +
    ggplot2::theme_minimal()
  if (!is.null(out_path)) ocp_save_plot(p, out_path, 8.5, 4.8)
  if (print_plot) print(p)
  invisible(p)
}

# Top-10 LOSO sites (first width+offset). Adds site_name if available.
ocp_top_loso_sites <- function(res, meta = NULL, min_frac = 1.0) {
  sb <- res$reps[[1]]$diagnostics$site_bucket
  ib <- interior_mask(res$coverage_tbl, min_frac)
  sb_int <- sb %>% dplyr::inner_join(ib, by="horizon_bucket") %>% dplyr::filter(is_interior)
  if (!nrow(sb_int)) return(tibble::tibble())
  mu <- sb_int %>% dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(mu = mean(normalized_entries, na.rm = TRUE),
                     n  = dplyr::n(), .groups="drop")
  jj <- sb_int %>% dplyr::left_join(mu, by="horizon_bucket") %>%
    dplyr::mutate(mean_minus_i = (n*mu - normalized_entries)/pmax(n-1, 1),
                  infl = abs(mean_minus_i - mu))
  site_infl <- jj %>% dplyr::group_by(site_id) %>%
    dplyr::summarise(influence = mean(infl, na.rm = TRUE), .groups="drop") %>%
    dplyr::arrange(dplyr::desc(influence)) %>% dplyr::slice_head(n = 10)
  if (!is.null(meta) && all(c("site_id","site_name") %in% names(meta))) {
    site_infl <- site_infl %>% dplyr::left_join(meta %>% dplyr::select(site_id, site_name), by="site_id")
  } else {
    site_infl <- site_infl %>% dplyr::mutate(site_name = NA_character_)
  }
  site_infl %>% dplyr::select(site_id, site_name, influence)
}

ocp_plot_top_loso_sites <- function(res, meta = NULL, min_frac = 1.0,
                                    green_thr = 1e-3, red_thr = 5e-3,
                                    out_path = NULL, print_plot = TRUE) {
  df <- ocp_top_loso_sites(res, meta, min_frac)
  if (!nrow(df)) {
    p <- ggplot2::ggplot() + ggplot2::labs(title="Top-10 LOSO sites", subtitle="No data", x=NULL, y="Influence") + ggplot2::theme_minimal()
    if (!is.null(out_path)) ocp_save_plot(p, out_path, 6.5, 4.5)
    if (print_plot) print(p); return(invisible(p))
  }
  df <- df %>%
    dplyr::mutate(label = dplyr::coalesce(site_name, as.character(site_id)),
                  status = ocp_status(influence, green_thr, red_thr))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(label, influence), y = influence, fill = status)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(green="#2ca02c", orange="#ff7f0e", red="#d62728", grey="#7f7f7f"), guide = "none") +
    ggplot2::geom_hline(yintercept = c(green_thr, red_thr), linetype = 2, alpha = 0.4) +
    ggplot2::labs(title = "Top-10 sites by LOSO influence (first width+offset)",
                  x = NULL, y = "Mean LOSO influence across interior buckets") +
    ggplot2::theme_minimal()
  if (!is.null(out_path)) ocp_save_plot(p, out_path, 7.2, 5.4)
  if (print_plot) print(p)
  invisible(p)
}


#' ocp_plot_loso_top_sites
#' 
#' Top-k sites by cumulative LOSO-influence
ocp_plot_loso_top_sites <- function(res, top_sites, k = 10L,
                                    out_path = NULL, print_plot = TRUE) {
  df <- top_sites %>%
    dplyr::slice_head(n = k) %>%
    dplyr::mutate(site_id = as.factor(site_id)) %>%
    dplyr::arrange(cum_infl)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = cum_infl, y = site_id)) +
    ggplot2::geom_col(fill = "grey40") +
    ggplot2::labs(
      title = sprintf("Top-%d sites by cumulative LOSO influence", k),
      subtitle = sprintf("width=%s · offset=%s", res$width, res$offset),
      x = "Cumulative influence (sum over buckets)", y = "site_id"
    ) +
    ggplot2::theme_minimal()
  if (!is.null(out_path)) ocp_save_plot(p, out_path, width = 7.2, height = 5.0)
  if (isTRUE(print_plot)) print(p)
  invisible(p)
}


#' ocp_plot_loso_vs_osi
#' 
#' LOSO vs OSI (for a given width): bars per offset
#' Expect: loso_summary_tbl has (width, offset, max_infl_abs) rows for that width
#'        osi_tbl has (width, offset, OSI_abs)
ocp_plot_loso_vs_osi <- function(width, loso_summary_tbl, osi_tbl,
                                 out_path = NULL, print_plot = TRUE) {
  df <- dplyr::inner_join(
    loso_summary_tbl %>% dplyr::filter(width == !!width) %>%
      dplyr::select(offset, max_infl_abs),
    osi_tbl %>% dplyr::filter(width == !!width) %>%
      dplyr::select(offset, OSI_abs),
    by = "offset"
  ) %>%
    tidyr::pivot_longer(c(max_infl_abs, OSI_abs), names_to = "metric", values_to = "value")
  
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(offset), y = value, fill = metric)) +
    ggplot2::geom_col(position = "dodge", width = 0.6) +
    ggplot2::labs(
      title = sprintf("LOSO (max) vs. OSI by offset (width = %s)", width),
      x = "Offset", y = "Absolute value", fill = NULL
    ) +
    ggplot2::theme_minimal()
  if (!is.null(out_path)) ocp_save_plot(p, out_path, width = 8, height = 4.5)
  if (isTRUE(print_plot)) print(p)
  invisible(p)
}


#' plot_offset_ribbon — Visualize offset sensitivity
#'
#' Builds, for a fixed width, the union of interior buckets across
#' offsets and plots the envelope ribbon [min, max] of Q̄(H) across
#' offsets with the mean line. A thin ribbon indicates robustness.
plot_offset_ribbon <- function(ens, width, min_frac = 1.0) {
  Rs <- purrr::keep(ens$results, ~ .x$width == width)
  if (!length(Rs)) {
    return(ggplot2::ggplot() + ggplot2::labs(title="Offset ribbon", subtitle="No results", x="PH", y="Q") + ggplot2::theme_minimal())
  }
  series <- purrr::map(Rs, function(r) {
    ib <- interior_mask(r$coverage_tbl, min_frac)
    r$pooled %>%
      dplyr::mutate(H = as.numeric(horizon_bucket)) %>%
      dplyr::inner_join(ib, by = "horizon_bucket") %>%
      dplyr::filter(is_interior) %>%
      dplyr::arrange(H) %>%
      dplyr::transmute(idx = dplyr::row_number(), H, Q = Q_bar, offset = r$offset)
  })
  series <- purrr::keep(series, ~ nrow(.x) > 0)
  if (!length(series)) {
    return(ggplot2::ggplot() + ggplot2::labs(
      title = "Offset sensitivity ribbon",
      subtitle = paste("width =", width, "· no interior buckets after filtering"),
      x = "Pseudo Horizon (median center)", y = "Normalized occupancy"
    ) + ggplot2::theme_minimal())
  }
  K <- min(vapply(series, nrow, integer(1)))
  long <- purrr::map_dfr(series, ~ dplyr::filter(.x, idx <= K))
  env <- long %>%
    dplyr::group_by(idx) %>%
    dplyr::summarise(H = stats::median(H), Q_min = min(Q), Q_max = max(Q), Q_mean = mean(Q), .groups="drop")
  ggplot2::ggplot(env, ggplot2::aes(H, Q_mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Q_min, ymax = Q_max), alpha = 0.20) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Offset sensitivity ribbon (rank-aligned)",
                  subtitle = paste("width =", width, "· interior ranked & aligned by order"),
                  x = "Pseudo Horizon (median center)", y = "Normalized occupancy") +
    ggplot2::theme_minimal()
}

#' plot_coverage — Coverage fraction per bucket
#'
#' Bar plot of \eqn{c(b)} showing which buckets are fully inside the
#' global PHU domain and which are edge-affected.
plot_coverage <- function(res) {
  # width for bars: compute a pleasant width automatically
  hb <- sort(as.numeric(res$coverage_tbl$horizon_bucket))
  bw <- if (length(hb) >= 2) min(diff(hb)) * 0.9 else 1
  ggplot(res$coverage_tbl, aes(x = as.numeric(horizon_bucket), y = cover_frac)) +
    geom_col(width = bw) +
    geom_hline(yintercept = 1, linetype = 2) +
    labs(title = "Bucket coverage fraction",
         subtitle = paste("width =", res$width, "· offset =", res$offset),
         x = "Pseudo Horizon (bucket)", y = "Coverage (0–1)") +
    theme_minimal()
}

#' plot_meta_share — Primary vs meta contributions over H
#'
#' Stacked columns of total mass per bucket split by assign_source
#' ("primary","meta") on interior buckets (by coverage threshold).
plot_meta_share <- function(res, min_frac = 1.0) {
  sb <- res$reps[[1]]$diagnostics$site_bucket
  ib <- interior_mask(res$coverage_tbl, min_frac)
  sb_int <- sb %>% dplyr::inner_join(ib, by = "horizon_bucket") %>%
    dplyr::filter(is_interior)
  
  has_cols <- all(c("primary_share","meta_share") %in% names(sb_int))
  
  if (!has_cols) {
    df <- sb_int %>%
      dplyr::group_by(horizon_bucket) %>%
      dplyr::summarise(total = sum(entries_per_period, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(source = "total", mass = total)
    
    return(
      ggplot2::ggplot(df, ggplot2::aes(as.numeric(horizon_bucket), mass, fill = source)) +
        ggplot2::geom_col() +
        ggplot2::labs(
          title = "Primary vs. meta contributions (interior)",
          subtitle = paste0("width = ", res$width, " · offset = ", res$offset,
                            " · note: primary/meta unavailable; showing total mass"),
          x = "Pseudo Horizon (bucket)", y = "Total mass"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::guides(fill = "none")
    )
  }
  
  df <- sb_int %>%
    dplyr::group_by(horizon_bucket) %>%
    dplyr::summarise(
      primary = sum(primary_share, na.rm = TRUE),
      meta    = sum(meta_share,    na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(c(primary, meta), names_to = "source", values_to = "mass")
  
  ggplot2::ggplot(df, ggplot2::aes(as.numeric(horizon_bucket), mass, fill = source)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::labs(
      title = "Primary vs. meta contributions (interior)",
      subtitle = paste("width =", res$width, "· offset =", res$offset),
      x = "Pseudo Horizon (bucket)", y = "Total mass"
    ) +
    ggplot2::theme_minimal()
}

#' ocp_build_minimal_package — One-click calibration + analysis bundle
#'
#' Runs a small calibration (synthetic uniform) and the real dataset
#' analysis end-to-end: ensemble curves, coverage/meta plots, and a set
#' of QC metrics, writing both figures and CSVs to `out_dir`.
#' Serves as a reproducible driver for typical study runs.
ocp_build_minimal_package <- function(
    # ---- Inputs ----
    harm_tab, meta, indiv = NULL, agg = NULL,
    out_dir             = "ocp_review",
    save_plots          = TRUE,
    save_csv            = TRUE,
    run_calibration     = TRUE,
    
    # ---- Calibration fixture knobs (kept small/fast) ----
    calib_sites         = 5L,
    calib_phases        = 20L,
    calib_total_per_site= 200L,
    calib_phase_len     = 10L,
    
    # ===================== Pass-through of run_stochastic_ensemble parameters =====================
    M                   = 5L,
    widths              = 10L,
    offsets             = 0,
    sampling            = c("deterministic_max","stochastic","fractional"),
    base_seed           = ocp_get_global_seed(),
    primary_profile     = weight_profile("uniform"),
    meta_profile        = founding_rise_taper(),
    meta_weight         = 0.25,
    meta_cap_fraction   = 0.30,
    normalization       = c("size","duration","size_duration","none"),
    region_col          = "cemetery_region",
    B_blocks            = 200L,
    use_t_df            = TRUE,
    edge_correction     = c("none","trim","coverage"),
    min_coverage_frac   = 1.0,
    # OSI thresholds
    thr_green           = 0.05,
    thr_orange          = 0.10,
    # --- Regional facets ---
    do_regional         = TRUE,
    regional_loess_span = 0.75,
    # --- Spatial slices ---
    do_spatial                = TRUE,
    spatial_selected_buckets  = c(30, 60, 90, 120),
    spatial_use_normalized    = FALSE,
    spatial_bins              = 7,
    spatial_bbox              = c(7, 15, 51, 58),
    spatial_point_color       = "red3",
    spatial_point_alpha       = 0.6,
    spatial_anim              = FALSE,
    spatial_anim_bucket_bin   = 10
) {
  # ---------------- tiny helpers ----------------
  stamp <- function(msg, t0) {
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)
    message(sprintf("[ocp_build_minimal_package] %s (%.2fs)", msg, dt))
  }
  osi_color <- function(x, g, o) {
    if (!is.finite(x)) return("#888888")
    if (x <= g) "#2ca02c" else if (x <= o) "#ff7f0e" else "#d62728"
  }
  
  sampling        <- match.arg(sampling)
  normalization   <- match.arg(normalization)
  edge_correction <- match.arg(edge_correction)
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  t0 <- Sys.time()
  
  # ---------------- Calibration (optional & fast) ----------------
  calib <- NULL
  if (isTRUE(run_calibration)) {
    fx <- ocp_make_uniform_fixture(
      n_sites = calib_sites, n_phases = calib_phases,
      per_site_total = calib_total_per_site,
      phase_length_phu = calib_phase_len,
      start_phu = 0
    )
    
    calib <- run_stochastic_ensemble(
      M = 1,
      mort_indiv = fx$indiv, mort_agg = fx$agg,
      harm_tab   = fx$harm_tab, meta_data = fx$meta,
      widths              = calib_phase_len,
      offsets             = offsets,
      sampling            = "fractional",
      base_seed           = base_seed,
      primary_profile     = weight_profile("uniform"),
      meta_profile        = weight_profile("uniform"),
      meta_weight         = 0,
      meta_cap_fraction   = NULL,
      normalization       = "size_duration",
      region_col          = "cemetery_region",
      B_blocks            = min(50L, B_blocks),
      use_t_df            = use_t_df,
      edge_correction     = edge_correction,
      min_coverage_frac   = min_coverage_frac,
      thr_green           = thr_green,
      thr_orange          = thr_orange,
      do_regional         = FALSE,
      do_spatial          = FALSE
    )
    
    if (length(calib$results)) {
      
      # --- Calibration metrics (flatness) ---
      cal_metrics <- calib$results %>%
        purrr::map_dfr(~ metric_flatness(.x, min_frac = min_coverage_frac))
      if (isTRUE(save_csv)) {
        ocp_write_csv(as.data.frame(cal_metrics),
                      file.path(out_dir, "calibration_metrics.csv"))
      }
      stamp("calibration metrics written", t0)
      
      # --- Calibration offset sensitivity (rank-aligned width = calib_phase_len) ---
      cal_osi <- metric_offset_sensitivity(
        calib, width = as.numeric(calib_phase_len), min_frac = min_coverage_frac
      )
      if (isTRUE(save_csv)) {
        ocp_write_csv(as.data.frame(cal_osi),
                      file.path(out_dir, "calibration_offset_sensitivity.csv"))
      }
      stamp("calibration OSI written", t0)
      
      
      if (isTRUE(save_plots)) {
        ocp_save_plot(
          plot_offset_ribbon(calib, width = calib_phase_len, min_frac = min_coverage_frac),
          file.path(out_dir, "calibration_offset_ribbon.png"), width = 8, height = 4.5
        )
        ocp_save_plot(
          plot_coverage(calib$results[[1]]),
          file.path(out_dir, "calibration_coverage_offset1.png"), width = 8, height = 4.5
        )
        # Meta-dependence for calibration (falls back to total-mass if primary/meta missing)
        ocp_save_plot(
          plot_meta_share(calib$results[[1]], min_frac = min_coverage_frac),
          file.path(out_dir, sprintf("calibration_meta_share_offset%d.png", as.integer(calib$results[[1]]$offset))),
          width = 8, height = 4.5
        )
      }
      stamp("calibration plots saved", t0)
    }
  }
  
  # ---------------- Real dataset ensemble ----------------
  ens <- run_stochastic_ensemble(
    M                   = M,
    mort_indiv          = indiv,
    mort_agg            = agg,
    harm_tab            = harm_tab,
    meta_data           = meta,
    widths              = widths,
    offsets             = offsets,
    sampling            = sampling,
    base_seed           = base_seed,
    primary_profile     = primary_profile,
    meta_profile        = meta_profile,
    meta_weight         = meta_weight,
    meta_cap_fraction   = meta_cap_fraction,
    normalization       = normalization,
    region_col          = region_col,
    B_blocks            = B_blocks,
    use_t_df            = use_t_df,
    edge_correction     = edge_correction,
    min_coverage_frac   = min_coverage_frac,
    thr_green           = thr_green,
    thr_orange          = thr_orange,
    do_regional         = do_regional,
    regional_loess_span = regional_loess_span,
    do_spatial                = do_spatial,
    spatial_selected_buckets  = spatial_selected_buckets,
    spatial_use_normalized    = spatial_use_normalized,
    spatial_bins              = spatial_bins,
    spatial_bbox              = spatial_bbox,
    spatial_point_color       = spatial_point_color,
    spatial_point_alpha       = spatial_point_alpha,
    spatial_anim              = spatial_anim,
    spatial_anim_bucket_bin   = spatial_anim_bucket_bin
  )
  stamp("ensemble finished", t0)
  
  if (!length(ens$results)) {
    warning("[ocp_build_minimal_package] No results; returning.")
    return(invisible(list(out_dir = out_dir, calibration = calib, ensemble = ens)))
  }
  
  # Representative (first width+offset)
  res_first     <- ens$results[[1]]
  primary_width <- as.numeric(ens$grid$width[[1]])
  
  # ---------------- Plots: main + ribbon (preserved) ----------------
  if (isTRUE(save_plots)) {
    if (!is.null(res_first$overall_plot)) {
      ocp_save_plot(
        res_first$overall_plot,
        file.path(out_dir, sprintf("main_curve_width%d_offset%d.png", as.integer(res_first$width), as.integer(res_first$offset)))
      )
    }
    ocp_save_plot(
      plot_offset_ribbon(ens, width = primary_width, min_frac = min_coverage_frac),
      file.path(out_dir, sprintf("offset_ribbon_width%d.png", as.integer(primary_width)))
    )
    stamp("main + ribbon plots saved", t0)
  }
  
  # ---------------- Median across offsets (preserved; with fallback) ----------------
  # Function may exist in your code; if not, we do a fallback.
  if (exists("ocp_median_curve_across_offsets", mode = "function")) {
    med_tbl <- ocp_median_curve_across_offsets(ens, width = primary_width, min_frac = min_coverage_frac)
  } else {
    # Fallback: use ens$diag$osi[[width]]$median_curve (rank) and map ranks to the interior horizon buckets of res_first
    wkey <- as.character(primary_width)
    med_tbl <- tibble::tibble()
    if (!is.null(ens$diag$osi[[wkey]]) && nrow(ens$diag$osi[[wkey]]$median_curve)) {
      mc <- ens$diag$osi[[wkey]]$median_curve # rank, Q
      ib <- interior_mask(res_first$coverage_tbl, min_frac = min_coverage_frac)
      ref <- res_first$pooled %>%
        dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
        dplyr::inner_join(ib, by = "horizon_bucket") %>%
        dplyr::filter(is_interior) %>%
        dplyr::arrange(horizon_bucket) %>%
        dplyr::mutate(rank = dplyr::row_number()) %>%
        dplyr::select(rank, horizon_bucket)
      med_tbl <- mc %>%
        dplyr::inner_join(ref, by = "rank") %>%
        dplyr::select(horizon_bucket, Q_median = Q)
    }
  }
  
  if (nrow(med_tbl)) {
    if (isTRUE(save_csv)) ocp_write_csv(as.data.frame(med_tbl), file.path(out_dir, sprintf("median_across_offsets_width%d.csv", as.integer(primary_width))))
    if (isTRUE(save_plots)) {
      if (exists("ocp_plot_median_across_offsets", mode = "function")) {
        ocp_plot_median_across_offsets(
          ens, width = primary_width, min_frac = min_coverage_frac,
          out_path = file.path(out_dir, sprintf("median_across_offsets_width%d.png", as.integer(primary_width))),
          print_plot = FALSE
        )
      } else {
        # Minimal fallback plot
        p_med <- ggplot2::ggplot(med_tbl, ggplot2::aes(as.numeric(horizon_bucket), Q_median)) +
          ggplot2::geom_line() + ggplot2::geom_point() +
          ggplot2::labs(title = "Median across offsets (rank-aligned, mapped to first (w,o))",
                        x = "Pseudo Horizon (bucket)", y = "Median normalized occupancy") +
          ggplot2::theme_minimal()
        ocp_save_plot(p_med, file.path(out_dir, sprintf("median_across_offsets_width%d.png", as.integer(primary_width))))
      }
    }
    stamp("median-across-offsets written", t0)
  }
  
  # ---------------- OSI ranked CSV + status plot (preserved; with fallback) ----------------
    # Fallback: use ens$diag$osi[[width]]
    wkey <- as.character(primary_width)
    os <- ens$diag$osi[[wkey]]
    overall <- os$overall
    per_off <- os$per_offset
    osi_df <- list(overall = overall, per_offset = per_off)
    
  if (isTRUE(save_csv)) {
    # Write both overall and per_offset if available
    if (is.list(osi_df) && !is.null(osi_df$overall)) {
      ocp_write_csv(as.data.frame(osi_df$overall), file.path(out_dir, sprintf("osi_ranked_width%d_overall.csv", as.integer(primary_width))))
      if (!is.null(osi_df$per_offset) && nrow(osi_df$per_offset)) {
        ocp_write_csv(as.data.frame(osi_df$per_offset), file.path(out_dir, sprintf("osi_ranked_width%d_per_offset.csv", as.integer(primary_width))))
      }
    } else {
      ocp_write_csv(as.data.frame(osi_df), file.path(out_dir, sprintf("osi_ranked_width%d.csv", as.integer(primary_width))))
    }
  if (isTRUE(save_plots)) {
    osi_w <- ens$diag$osi[[as.character(primary_width)]]
    ocp_plot_osi_status(
      osi_w$per_offset,
      metric = "OSIrel_RMSE",
      thr_green = thr_green,
      thr_orange = thr_orange,
      out_path = file.path(out_dir, sprintf("osi_ranked_width%d_per_offset.png", primary_width))
    )
  }
  stamp("OSI written", t0)
  }
  
  
  # ---------------- LOSO (first width+offset) bucketwise + top-10 (preserved) ----------------
  if (isTRUE(save_plots)) {
    if (exists("ocp_plot_loso_bucketwise", mode = "function")) {
      ocp_plot_loso_bucketwise(
        res_first, min_frac = min_coverage_frac,
        out_path = file.path(out_dir, sprintf("loso_bucketwise_width%d_offset%d.png",
                                              as.integer(res_first$width), as.integer(res_first$offset)))
      )
    }
    if (exists("ocp_plot_top_loso_sites", mode = "function")) {
      ocp_plot_top_loso_sites(
        res_first, meta = meta, min_frac = min_coverage_frac,
        out_path = file.path(out_dir, sprintf("loso_top10_width%d_offset%d.png",
                                              as.integer(res_first$width), as.integer(res_first$offset)))
      )
    }
    stamp("LOSO plots saved", t0)
  }
  
  # ---------------- LOSO median across offsets (fixed width) (preserved; with fallback) ----------------
  if (exists("ocp_median_loso_across_offsets", mode = "function")) {
    loso_med <- ocp_median_loso_across_offsets(ens, width = primary_width, min_frac = min_coverage_frac)
  } else {
    wkey <- as.character(primary_width)
    loso_med <- ens$diag$loso_median[[wkey]] %||% tibble::tibble()
  }
  if (nrow(loso_med)) {
    if (isTRUE(save_csv)) ocp_write_csv(as.data.frame(loso_med), file.path(out_dir, sprintf("loso_median_across_offsets_width%d.csv", as.integer(primary_width))))
    if (isTRUE(save_plots)) {
      if (exists("ocp_plot_median_loso_across_offsets", mode = "function")) {
        ocp_plot_median_loso_across_offsets(
          ens, width = primary_width, min_frac = min_coverage_frac,
          out_path = file.path(out_dir, sprintf("loso_median_across_offsets_width%d.png", as.integer(primary_width)))
        )
      } else {
        p_lmed <- ggplot2::ggplot(loso_med, ggplot2::aes(rank, infl_mean_med)) +
          ggplot2::geom_line() + ggplot2::geom_point() +
          ggplot2::labs(title = "LOSO median across offsets (by rank)",
                        x = "Ranked bucket index", y = "Median LOSO mean influence") +
          ggplot2::theme_minimal()
        ocp_save_plot(p_lmed, file.path(out_dir, sprintf("loso_median_across_offsets_width%d.png", as.integer(primary_width))))
      }
    }
    stamp("LOSO median saved", t0)
  }
  
  # ---------------- keep full CSVs for *all* widths/offsets (preserved) ----------------
  if (isTRUE(save_csv)) {
    flat_tbl <- ens$results %>% purrr::map_dfr(~ metric_flatness(.x, min_frac = min_coverage_frac))
    repstab  <- ens$results %>% purrr::map_dfr(~ metric_replicate_stability(.x, min_frac = min_coverage_frac))
    meta_dep <- ens$results %>% purrr::map_dfr(~ metric_meta_dependence(.x, min_frac = min_coverage_frac))
    loso_all <- ens$results %>% purrr::map_dfr(function(r) {
      if (!exists("ocp_loso_bucketwise", mode = "function")) return(tibble::tibble())
      bw <- ocp_loso_bucketwise(r, min_frac = min_coverage_frac)
      if (!nrow(bw)) return(tibble::tibble())
      bw %>% dplyr::mutate(width = r$width, offset = r$offset)
    })
    
    ocp_write_csv(as.data.frame(flat_tbl), file.path(out_dir, "all_flatness_by_offset.csv"))
    ocp_write_csv(as.data.frame(repstab),  file.path(out_dir, "all_replicate_stability.csv"))
    ocp_write_csv(as.data.frame(meta_dep), file.path(out_dir, "all_meta_dependence.csv"))
    ocp_write_csv(as.data.frame(loso_all), file.path(out_dir, "all_loso_bucketwise.csv"))
    stamp("all CSVs written", t0)
    
    loso_inf <- ens$results %>%
      purrr::map_dfr(~ metric_loso_influence(.x, min_frac = min_coverage_frac))
    ocp_write_csv(as.data.frame(loso_inf), file.path(out_dir, "all_loso_influence.csv"))
    stamp("LOSO influence CSV written", t0)
  }
  
  # ---------------- Support table for representative (first w,o) (preserved) ----------------
  if (isTRUE(save_csv)) {
    sb_first <- res_first$reps[[1]]$diagnostics$site_bucket
    sup_tbl <- sb_first %>%
      dplyr::group_by(horizon_bucket) %>%
      dplyr::summarise(
        n_sites      = dplyr::n_distinct(site_id),
        total_weight = sum(entries_per_period, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::left_join(res_first$coverage_tbl, by = "horizon_bucket")
    ocp_write_csv(as.data.frame(sup_tbl), file.path(out_dir, "support_by_bucket_first_width_offset.csv"))
    stamp("support table written", t0)
  }
  
  # ---------------- Meta-dependence plot for first (w,o) (preserved/new) ----------------
  if (isTRUE(save_plots)) {
    ocp_save_plot(
      plot_meta_share(res_first, min_frac = min_coverage_frac),
      file.path(out_dir, sprintf("meta_share_width%d_offset%d.png", as.integer(res_first$width), as.integer(res_first$offset))),
      width = 8, height = 4.5
    )
  }
  
  # ---------------- Manifest ----------------
  manifest <- tibble::tibble(
    timestamp = as.character(Sys.time()),
    out_dir   = out_dir,
    # calibration knobs
    run_calibration, calib_sites, calib_phases, calib_total_per_site, calib_phase_len,
    # ensemble knobs
    M, widths = paste(widths, collapse = ","), offsets = paste(offsets, collapse = ","),
    sampling, base_seed, meta_weight, meta_cap_fraction,
    normalization, region_col, B_blocks, use_t_df,
    edge_correction, min_coverage_frac,
    thr_green, thr_orange,
    do_regional, regional_loess_span,
    do_spatial,
    spatial_selected_buckets = paste(spatial_selected_buckets, collapse = ","),
    spatial_use_normalized, spatial_bins,
    spatial_bbox = paste(spatial_bbox, collapse = ","),
    spatial_point_color, spatial_point_alpha,
    spatial_anim, spatial_anim_bucket_bin
  )
  ocp_write_csv(as.data.frame(manifest), file.path(out_dir, "manifest.csv"))
  stamp("manifest written", t0)
  
  invisible(list(
    out_dir = out_dir,
    calibration = calib,
    ensemble = ens,
    figures = list.files(out_dir, full.names = TRUE)
  ))
}



#' plot_harmonization_table — Gantt-like phase plot
#'
#' Visualizes the harmonization table as horizontal bars (per system)
#' over PHU with optional fade-in/out error bars, alternating row bands,
#' vertical grid (width+offset), wrapped labels with adaptive sizing, and
#' optional richtext. Autosizes figure if width/height not provided.
plot_harmonization_table <- function(
    table,
    cols = list(
      system  = "system_name",
      start   = "horizon_start",
      end     = "horizon_end",
      fade_in = "fade_in_start",
      fade_out= "fade_out_end",
      label   = "phase_name"
    ),
    alternate_bands   = TRUE,
    band_alpha        = 0.06,
    show_grid         = TRUE,
    grid_width        = 10L,
    grid_offset       = 0,
    wrap_labels       = TRUE,
    wrap_phu_per_char = 1.8,
    adaptive_size     = TRUE,
    size_min          = 2.2,
    size_max          = 3.6,
    stagger_labels    = TRUE,
    y_stagger         = 0.22,
    angle_deg         = 0,
    use_richtext      = FALSE,
    richtext_fill     = NA,
    richtext_alpha    = 0.0,
    font_family       = NULL,
    out_path          = "graphs/harmonization_table.png",
    out_width         = NULL,
    out_height        = NULL,
    dpi               = 300,
    print_plot        = TRUE,
    debug             = FALSE
) {
  # ---- Column standardization ----
  if (!is.null(cols)) {
    ren <- list(
      system_name   = cols$system,
      horizon_start = cols$start,
      horizon_end   = cols$end,
      fade_in_start = cols$fade_in,
      fade_out_end  = cols$fade_out,
      phase_label   = cols$label
    )
    ren <- ren[!vapply(ren, is.null, logical(1))]
    missing <- setdiff(unname(unlist(ren)), names(table))
    if (length(missing)) stop("In `cols=`, not found: ", paste(missing, collapse = ", "))
    dat <- dplyr::rename(table, !!!rlang::set_names(rlang::syms(unname(unlist(ren))), names(ren)))
  } else {
    dat <- table
    if (!"phase_label" %in% names(dat)) {
      if ("phase_name" %in% names(dat))          dat <- dplyr::rename(dat, phase_label = phase_name)
      else if ("relative_chron" %in% names(dat)) dat <- dplyr::rename(dat, phase_label = relative_chron)
      else dat$phase_label <- NA_character_
    }
  }
  if (!"phase_label" %in% names(dat)) {
    lbl <- if ("phase_name" %in% names(dat)) dat$phase_name else
      if ("relative_chron" %in% names(dat)) dat$relative_chron else
        if ("phase_id" %in% names(dat)) as.character(dat$phase_id) else
          rep("", nrow(dat))
    dat$phase_label <- as.character(lbl)
  } else {
    dat$phase_label <- as.character(dat$phase_label)
  }
  
  # ---- Requirements & numeric coercion ----
  req <- c("system_name","horizon_start","horizon_end")
  miss <- setdiff(req, names(dat))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse=", "))
  
  numc <- intersect(c("horizon_start","horizon_end","fade_in_start","fade_out_end"), names(dat))
  dat[numc] <- lapply(dat[numc], function(x) suppressWarnings(as.numeric(x)))
  dat <- dplyr::filter(dat, is.finite(horizon_start) & is.finite(horizon_end))
  
  # ---- System ordering & per-phase fields ----
  sys_order <- dat %>%
    dplyr::group_by(system_name) %>%
    dplyr::summarise(earliest_start = min(horizon_start, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(earliest_start, system_name) %>%
    dplyr::mutate(system_numeric = dplyr::row_number())
  
  dat <- dat %>% dplyr::left_join(sys_order, by = "system_name") %>%
    dplyr::group_by(system_name) %>%
    dplyr::arrange(horizon_start, .by_group = TRUE) %>%
    dplyr::mutate(
      phase_idx = dplyr::row_number(),
      span_phu  = pmax(0, horizon_end - horizon_start),
      x_mid     = (horizon_start + horizon_end)/2
    ) %>% dplyr::ungroup()
  
  # ---- Label wrapping & sizing ----
  lab <- dat$phase_label
  if (wrap_labels) {
    char_budget <- pmax(4, floor(dat$span_phu / wrap_phu_per_char))
    lab <- purrr::map2_chr(dat$phase_label, char_budget,
                           ~ stringr::str_wrap(.x %||% "", width = as.integer(.y)))
  }
  dat$lab_wrapped <- lab
  
  if (adaptive_size) {
    span_max <- if (any(is.finite(dat$span_phu))) max(dat$span_phu, na.rm = TRUE) else 1
    base_sz  <- size_min + (pmin(dat$span_phu, span_max) / span_max) * (size_max - size_min)
    base_sz  <- base_sz - pmax(0, (nchar(gsub("\n","",dat$lab_wrapped)) - 18)) * 0.03
    dat$lab_size <- pmin(size_max, pmax(size_min, base_sz))
  } else {
    dat$lab_size <- (size_min + size_max) / 2
  }
  
  dat$y_lab <- dat$system_numeric
  if (stagger_labels) {
    dat$y_lab <- dat$system_numeric + ifelse(dat$phase_idx %% 2 == 0, y_stagger, -y_stagger)
  }
  
  thr <- max(1e-6, 0.2 * grid_width)
  dat <- dat %>%
    dplyr::group_by(system_name) %>%
    dplyr::arrange(horizon_start, .by_group = TRUE) %>%
    dplyr::mutate(
      mid_prev = dplyr::lag(x_mid),
      mid_next = dplyr::lead(x_mid),
      close_l  = is.finite(mid_prev) & (x_mid - mid_prev) < thr,
      close_r  = is.finite(mid_next) & (mid_next - x_mid) < thr,
      nudge0   = ifelse(close_l, +1, 0) + ifelse(close_r, -1, 0),
      pad      = pmax(0.15 * span_phu, pmin(2, 0.15 * span_phu)),
      max_nudge= pmax(0, (span_phu/3)),
      x_mid2   = pmin(horizon_end - pad, pmax(horizon_start + pad,
                                              x_mid + pmax(-max_nudge, pmin(max_nudge, nudge0 * thr/2))))
    ) %>% dplyr::ungroup()
  
  # ---- Grid & bands ----
  gs <- min(dat$horizon_start, na.rm = TRUE)
  ge <- max(dat$horizon_end,   na.rm = TRUE)
  b0 <- floor((gs - grid_offset)/grid_width) * grid_width + grid_offset
  b1 <- floor(((ge - 1e-9) - grid_offset)/grid_width) * grid_width + grid_offset
  grid_x <- seq(b0, b1, by = grid_width)
  
  band_df <- sys_order %>%
    dplyr::mutate(is_alt = (system_numeric %% 2 == 0),
                  y0 = system_numeric - 0.5,
                  y1 = system_numeric + 0.5)
  
  # ---- Font handling: coerce to single string; omit if NULL ----
  ff <- if (!is.null(font_family) && length(font_family) >= 1L) as.character(font_family[[1L]]) else NULL
  base_ff <- ff %||% ""  # theme is happy with ""
  
  # ---- Plot ----
  p <- ggplot2::ggplot() +
    ggplot2::theme_minimal(base_family = base_ff) +
    ggplot2::labs(x = "Pseudo Horizon (PHU)", y = "Systems",
                  title = "Harmonization table",
                  subtitle = sprintf("Grid width=%s · offset=%s", grid_width, grid_offset)) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(size = 11),
      axis.text.x  = ggplot2::element_text(size = 10),
      axis.title   = ggplot2::element_text(size = 13),
      panel.border = ggplot2::element_rect(colour = "grey80", fill = NA),
      plot.margin  = grid::unit(c(6, 8, 6, 8), "pt")
    ) +
    ggplot2::scale_y_continuous(
      breaks = sys_order$system_numeric,
      labels = sys_order$system_name,
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::coord_cartesian(xlim = c(gs, ge), clip = "off")
  
  if (alternate_bands && nrow(band_df)) {
    p <- p + ggplot2::geom_rect(
      data = dplyr::filter(band_df, is_alt),
      ggplot2::aes(xmin = gs, xmax = ge, ymin = y0, ymax = y1),
      inherit.aes = FALSE, fill = "grey40", alpha = band_alpha, colour = NA
    )
  }
  if (show_grid && length(grid_x)) {
    p <- p + ggplot2::geom_vline(xintercept = grid_x, colour = "grey85", linewidth = 0.3)
  }
  p <- p + ggplot2::geom_segment(
    data = dat,
    ggplot2::aes(y = system_numeric, yend = system_numeric,
                 x = horizon_start, xend = horizon_end),
    linewidth = 5, colour = "grey70"
  )
  
  sep <- dat %>%
    dplyr::group_by(system_name, system_numeric) %>%
    dplyr::arrange(horizon_start, .by_group = TRUE) %>%
    dplyr::mutate(next_start = dplyr::lead(horizon_start)) %>%
    dplyr::filter(!is.na(next_start)) %>%
    dplyr::select(system_numeric, sep_x = horizon_end)
  if (nrow(sep)) {
    p <- p + ggplot2::geom_segment(
      data = sep,
      ggplot2::aes(x = sep_x, xend = sep_x, y = system_numeric - 0.22, yend = system_numeric + 0.22),
      linewidth = 0.4, colour = "grey55"
    )
  }
  
  if ("fade_in_start" %in% names(dat)) {
    f_in <- dplyr::filter(dat, is.finite(fade_in_start))
    if (nrow(f_in)) {
      p <- p + ggplot2::geom_errorbar(
        data = f_in,
        ggplot2::aes(xmin = fade_in_start, xmax = pmin(horizon_start, horizon_end), y = system_numeric),
        width = 0.15
      )
    }
  }
  if ("fade_out_end" %in% names(dat)) {
    f_out <- dplyr::filter(dat, is.finite(fade_out_end))
    if (nrow(f_out)) {
      p <- p + ggplot2::geom_errorbar(
        data = f_out,
        ggplot2::aes(xmin = pmax(horizon_start, horizon_end), xmax = fade_out_end, y = system_numeric),
        width = 0.15
      )
    }
  }
  
  # ---- Label layer (omit 'family' when NULL) ----
  if (use_richtext && !requireNamespace("ggtext", quietly = TRUE)) {
    warning("ggtext not installed; falling back to plain geom_text().")
    use_richtext <- FALSE
  }
  if (use_richtext) {
    args <- list(
      data = dat,
      mapping = ggplot2::aes(x = x_mid2, y = y_lab, label = lab_wrapped),
      label.size = 0,
      fill = richtext_fill,
      alpha = richtext_alpha,
      size = dat$lab_size,
      angle = angle_deg,
      vjust = 0.5
    )
    if (!is.null(ff)) args$family <- ff
    p <- p + do.call(ggtext::geom_richtext, args)
  } else {
    args <- list(
      data = dat,
      mapping = ggplot2::aes(x = x_mid2, y = y_lab, label = lab_wrapped),
      size = dat$lab_size,
      angle = angle_deg,
      vjust = 0.5
    )
    if (!is.null(ff)) args$family <- ff
    p <- p + do.call(ggplot2::geom_text, args)
  }
  
  # ---- Auto size if not provided ----
  if (is.null(out_width) || is.null(out_height)) {
    n_sys  <- nrow(sys_order)
    ph_per <- dat %>% dplyr::count(system_name, name = "k") %>%
      dplyr::summarise(m = max(k)) %>% dplyr::pull(m)
    ph_per <- ifelse(is.finite(ph_per), ph_per, 10)
    out_width  <- out_width  %||% (8 + 0.30 * ph_per)
    out_height <- out_height %||% (2.8 + 0.55 * n_sys)
  }
  
  # ---- Save & print ----
  if (!is.null(out_path)) {
    dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
    ocp_save_plot(p, out_path, width = out_width, height = out_height, dpi = dpi)
    if (debug) message(sprintf("[plot_harmonization_table] Saved to %s (%.2f×%.2f in, %ddpi)",
                               out_path, out_width, out_height, dpi))
  }
  if (print_plot) print(p)
  invisible(p)
}

#' ocp_pool_region_curves — Pool region×bucket across replicates
#'
#' For each region and bucket, pools replicate means via Rubin’s rules
#' (Q̄, W, B, T, df, CI), optionally scaling by coverage if edge method
#' is "coverage". Requires `region_col` in `site_bucket`.
ocp_pool_region_curves <- function(res,
                                   region_col = "cemetery_region",
                                   min_coverage_frac = 1.0,
                                   use_t_df = TRUE) {
  stopifnot(is.list(res), "reps" %in% names(res), length(res$reps) >= 1)
  if (!region_col %in% names(res$reps[[1]]$diagnostics$site_bucket)) {
    stop(sprintf("Column '%s' not found in site_bucket.", region_col))
  }
  ib <- interior_mask(res$coverage_tbl, min_frac = min_coverage_frac)
  
  # Per-replicate region×bucket means + within-var (var/n) over sites in region
  reg_long <- purrr::imap_dfr(res$reps, function(rep, irep) {
    sb <- rep$diagnostics$site_bucket %>%
      dplyr::mutate(horizon_bucket = as.numeric(horizon_bucket)) %>%
      dplyr::inner_join(ib, by = "horizon_bucket") %>%
      dplyr::filter(is_interior)
    
    sb %>%
      dplyr::group_by(.data[[region_col]], horizon_bucket) %>%
      dplyr::summarise(
        mean_entries = mean(normalized_entries, na.rm = TRUE),
        n_sites      = dplyr::n_distinct(site_id),
        var_hat_reg  = {
          v <- stats::var(normalized_entries, na.rm = TRUE)
          if (!is.finite(v)) 0 else v / pmax(n_sites, 1)  # variance of the mean
        },
        .groups = "drop"
      ) %>%
      dplyr::transmute(
        region = as.character(.data[[region_col]]),
        horizon_bucket = as.numeric(horizon_bucket),
        mean_entries, var_hat_reg, .rep = irep
      )
  })
  
  if (nrow(reg_long) == 0) {
    return(tibble::tibble(
      region = character(), horizon_bucket = numeric(),
      M = integer(), Q_bar = numeric(), W = numeric(), B = numeric(),
      T = numeric(), df = numeric(), SE = numeric(), CI_lower = numeric(), CI_upper = numeric()
    ))
  }
  
  # Pool across replicates per region×bucket
  pooled <- reg_long %>%
    dplyr::group_by(region, horizon_bucket) %>%
    dplyr::summarise(
      M     = dplyr::n(),
      Q_bar = mean(mean_entries, na.rm = TRUE),
      W     = mean(var_hat_reg,  na.rm = TRUE),
      B     = stats::var(mean_entries, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      T  = W + (1 + 1/M) * B,
      df = dplyr::if_else(B > 0,
                          (M - 1) * (1 + W/((1 + 1/M) * B))^2,
                          Inf),
      SE = sqrt(T),
      crit = if (use_t_df) stats::qt(0.975, df = pmax(df, 1)) else 1.96,
      CI_lower = Q_bar - crit * SE,
      CI_upper = Q_bar + crit * SE
    )
  pooled
}

#' ocp_plot_region_curves — Faceted regional pooled trends
#'
#' Plots pooled region curves with points, CI error bars, and a LOESS
#' smoother; facets by region (free y-scale). Designed to accompany
#' the main pooled curve with regional detail.
ocp_plot_region_curves <- function(res,
                                   region_col = "cemetery_region",
                                   min_coverage_frac = 1.0,
                                   use_t_df = TRUE,
                                   loess_span = 0.75,
                                   out_path = NULL,
                                   print_plot = TRUE) {
  pooled_reg <- ocp_pool_region_curves(res, region_col, min_coverage_frac, use_t_df)
  
  if (nrow(pooled_reg) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::labs(
        title = "Regional trends",
        subtitle = "No interior buckets after coverage filtering.",
        x = "Pseudo Horizon (bucket)", y = "Normalized occupancy"
      ) + ggplot2::theme_minimal()
    if (!is.null(out_path)) ocp_save_plot(p, out_path, width = 9, height = 5)
    if (print_plot) print(p)
    return(invisible(p))
  }
  
  p <- ggplot2::ggplot(pooled_reg,
                       ggplot2::aes(x = as.numeric(horizon_bucket), y = Q_bar)) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, alpha = 0.7) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, span = loess_span) +
    ggplot2::facet_wrap(~ region, scales = "free_y") +
    ggplot2::labs(
      title = "Regional temporal trends (pooled across replicates)",
      subtitle = sprintf("width=%s · offset=%s · edge=%s (min=%.2f)",
                         res$width, res$offset, res$edge_corr %||% "none",
                         ifelse(is.na(res$min_cover), 0, res$min_cover)),
      x = "Pseudo Horizon (bucket)", y = "Normalized occupancy"
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(out_path)) ocp_save_plot(p, out_path, width = 9.5, height = 6)
  if (print_plot) print(p)
  invisible(p)
}

#' ocp_spatial_slices — Spatial density & point-size maps by bucket
#'
#' For selected horizon buckets, overlays 2D kernel density (normalized
#' density via `stat_density_2d(contour_var="ndensity")`) on a world
#' base map (sf), and draws site points sized by either
#' `entries_per_period` or `normalized_entries`. Can assemble a combined
#' grid with per-panel legends and export an animation sweeping through
#' horizon time (gganimate). Also writes a `site_bucket.csv` with
#' columns: site_id, site_name (if provided via `meta_lookup`),
#' entries_per_period, normalized_entries, coord_x, coord_y, cemetery_size.
ocp_spatial_slices <- function(
    res,
    meta_lookup         = NULL,                 # <— add: optional meta df with site_id, site_name
    selected_buckets    = c(30, 60, 90, 120),
    use_normalized      = FALSE,
    bins                = 7,
    bbox                = c(7, 15, 51, 58),
    point_color         = "red3",
    point_alpha         = 0.6,
    save_dir            = "graphs",
    filename_prefix     = "spatial_slice",
    width               = 7,
    height              = 6,
    dpi                 = 300,
    show_plots          = TRUE,
    make_combined       = TRUE,
    combined_filename   = "spatial_slices_grid.png",
    combined_title      = "Spatial density by horizon bucket",
    combined_ncol       = NULL,
    combined_width      = NULL,
    combined_height     = NULL,
    legend_rel_width    = 0.35,
    make_animation      = TRUE,
    anim_filename       = "spatial_density_anim.gif",
    anim_bucket_bin     = 10,
    anim_fps            = 12,
    anim_duration       = 8,
    # --- NEW (CSV export) ---
    export_site_bucket  = TRUE,                 # write site_bucket.csv?
    site_bucket_filename = "site_bucket.csv"    # filename inside save_dir
) {
  stopifnot(is.list(res), "reps" %in% names(res))
  sb <- res$reps[[1]]$diagnostics$site_bucket
  
  need <- c("horizon_bucket","entries_per_period","normalized_entries",
            "coord_x","coord_y","site_id","cemetery_size")
  miss <- setdiff(need, names(sb))
  if (length(miss)) stop("[ocp_spatial_slices] Missing columns in site_bucket: ",
                         paste(miss, collapse = ", "))
  
  # --- write site_bucket.csv (site_id, site_name, entries_per_period, normalized_entries, coord_x, coord_y, cemetery_size)
  if (isTRUE(export_site_bucket)) {
    sb_export <- sb %>%
      dplyr::select(site_id, entries_per_period, normalized_entries, coord_x, coord_y, cemetery_size) %>%
      dplyr::distinct()
    
    if (!is.null(meta_lookup) && all(c("site_id","site_name") %in% names(meta_lookup))) {
      sb_export <- sb_export %>%
        dplyr::left_join(meta_lookup %>% dplyr::select(site_id, site_name), by = "site_id") %>%
        dplyr::relocate(site_name, .after = site_id)
    } else {
      sb_export <- sb_export %>%
        dplyr::mutate(site_name = NA_character_) %>%
        dplyr::relocate(site_name, .after = site_id)
      warning("[ocp_spatial_slices] meta_lookup missing 'site_id'/'site_name'; writing site_name=NA.")
    }
    
    sb_export <- sb_export %>%
      dplyr::select(site_id, site_name, entries_per_period, normalized_entries, coord_x, coord_y, cemetery_size)
    
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    csv_path <- file.path(save_dir, site_bucket_filename)
    if (exists("ocp_write_csv", mode = "function")) {
      ocp_write_csv(sb_export, csv_path)
    } else {
      readr::write_csv(sb_export, csv_path)
    }
  }

  value_col <- if (isTRUE(use_normalized)) "normalized_entries" else "entries_per_period"
  
  slice_base <- sb %>%
    dplyr::filter(.data$horizon_bucket %in% selected_buckets) %>%
    dplyr::filter(is.finite(.data$coord_x), is.finite(.data$coord_y)) %>%
    dplyr::mutate(value = .data[[value_col]])
  
  if (!requireNamespace("rnaturalearth", quietly = TRUE) ||
      !requireNamespace("sf", quietly = TRUE)) {
    stop("Packages 'rnaturalearth' and 'sf' are required for map background.")
  }
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  aes_density <- ggplot2::aes(
    x = .data$coord_x, y = .data$coord_y,
    fill = after_stat(level), alpha = after_stat(level)
  )
  
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------- Static slices ----------
  base_plots <- list()
  for (hb in selected_buckets) {
    slice_data <- dplyr::filter(slice_base, .data$horizon_bucket == hb)
    if (nrow(slice_data) == 0) next
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = world, fill = "gray95", color = "gray80") +
      ggplot2::stat_density_2d(
        data = slice_data,
        mapping = aes_density,
        geom = "polygon",
        contour = TRUE,
        contour_var = "ndensity",
        bins = bins,
        # IMPORTANT: density only contributes to 'fill' legend, never 'size'
        show.legend = c(fill = TRUE, alpha = FALSE, size = FALSE)
      ) +
      ggplot2::scale_fill_viridis_c(option = "C") +
      ggplot2::scale_alpha(range = c(0.2, 0.7), guide = "none") +
      ggplot2::geom_point(
        data = slice_data,
        ggplot2::aes(x = .data$coord_x, y = .data$coord_y, size = .data$value),
        color = point_color, alpha = point_alpha,
        # point is the *only* contributor to size legend
        show.legend = c(size = TRUE)
      ) +
      ggplot2::scale_size_continuous(
        name = if (use_normalized) "Normalized occupancy" else "Expected burials",
        guide = ggplot2::guide_legend(
          override.aes = list(shape = 19, colour = point_color, fill = "white", alpha = 1, stroke = 0)
        )
      ) +
      ggplot2::labs(
        title = sprintf("Horizon Bucket: %s", hb),
        x = "Longitude", y = "Latitude", fill = "Density level"
      ) +
      ggplot2::coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE) +
      ggplot2::guides(
        alpha = "none",
        fill  = ggplot2::guide_colorbar(order = 1)
        # no guides(size=...) here; defined above in scale_size_continuous()
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position       = "right",
        legend.background     = ggplot2::element_rect(fill = "white", colour = NA),
        legend.box.background = ggplot2::element_rect(fill = "white", colour = NA),
        legend.key            = ggplot2::element_rect(fill = "white", colour = NA)
      )
    
    file_out <- file.path(save_dir, sprintf("%s_bucket_%s.png", filename_prefix, as.character(hb)))
    if (exists("ocp_save_plot", mode = "function")) {
      ocp_save_plot(p, file_out, width = width, height = height, dpi = dpi)
    } else {
      ggplot2::ggsave(file_out, p, width = width, height = height, dpi = dpi)
    }
    if (isTRUE(show_plots)) print(p)
    base_plots[[as.character(hb)]] <- p
  }
  
  # ---------- Combined grid with per-panel legends ----------
  combined_plot <- NULL
  if (make_combined && length(base_plots) > 1) {
    have_cow   <- requireNamespace("cowplot", quietly = TRUE)
    have_patch <- requireNamespace("patchwork", quietly = TRUE)
    
    if (!have_patch) {
      warning("[ocp_spatial_slices] 'patchwork' not installed; skipping combined grid.")
    } else if (!have_cow) {
      warning("[ocp_spatial_slices] 'cowplot' not installed; legends will remain inside panels.")
      np <- length(base_plots)
      if (is.null(combined_ncol)) combined_ncol <- if (np == 4) 2 else min(3, np)
      nrow <- ceiling(np / combined_ncol)
      combined_plot <- patchwork::wrap_plots(base_plots, ncol = combined_ncol) +
        patchwork::plot_annotation(title = combined_title)
      cw <- combined_width  %||% (width  * combined_ncol)
      ch <- combined_height %||% (height * nrow)
      file_grid <- file.path(save_dir, combined_filename)
      if (exists("ocp_save_plot", mode = "function")) {
        ocp_save_plot(combined_plot, file_grid, width = cw, height = ch, dpi = dpi)
      } else {
        ggplot2::ggsave(file_grid, combined_plot, width = cw, height = ch, dpi = dpi)
      }
      if (isTRUE(show_plots)) print(combined_plot)
    } else {
      with_leg <- lapply(base_plots, function(p) {
        leg   <- cowplot::get_legend(p)
        panel <- p + ggplot2::theme(legend.position = "none")
        cowplot::plot_grid(panel, leg, ncol = 2,
                           rel_widths = c(1, legend_rel_width),
                           align = "h", axis = "tb")
      })
      
      np <- length(with_leg)
      if (is.null(combined_ncol)) combined_ncol <- if (np == 4) 2 else min(3, np)
      nrow <- ceiling(np / combined_ncol)
      
      combined_plot <- patchwork::wrap_plots(with_leg, ncol = combined_ncol) +
        patchwork::plot_annotation(title = combined_title)
      
      cw <- combined_width  %||% (width  * combined_ncol * (1 + legend_rel_width * 0.35))
      ch <- combined_height %||% (height * nrow)
      
      file_grid <- file.path(save_dir, combined_filename)
      if (exists("ocp_save_plot", mode = "function")) {
        ocp_save_plot(combined_plot, file_grid, width = cw, height = ch, dpi = dpi)
      } else {
        ggplot2::ggsave(file_grid, combined_plot, width = cw, height = ch, dpi = dpi)
      }
      if (isTRUE(show_plots)) print(combined_plot)
    }
  }
  
  # ---------- Animation ----------
  anim_path <- NULL
  if (make_animation) {
    if (!requireNamespace("gganimate", quietly = TRUE)) {
      warning("[ocp_spatial_slices] 'gganimate' not installed; skipping animation.")
    } else {
      anim_base <- sb %>%
        dplyr::filter(is.finite(.data$coord_x), is.finite(.data$coord_y)) %>%
        dplyr::mutate(value = .data[[value_col]],
                      combined_bucket = floor(.data$horizon_bucket / anim_bucket_bin) * anim_bucket_bin) %>%
        dplyr::group_by(.data$combined_bucket) %>%
        dplyr::filter(dplyr::n() > 1) %>%
        dplyr::ungroup()
      
      if (nrow(anim_base) > 0) {
        p_anim <- ggplot2::ggplot(anim_base, ggplot2::aes(x = .data$coord_x, y = .data$coord_y)) +
          ggplot2::geom_sf(data = world, fill = "gray95", color = "gray80", inherit.aes = FALSE) +
          ggplot2::coord_sf(xlim = c(bbox[1], bbox[2]), ylim = c(bbox[3], bbox[4]), expand = FALSE) +
          ggplot2::stat_density_2d(
            ggplot2::aes(fill = after_stat(level), alpha = after_stat(level)),
            geom = "polygon", contour = TRUE, contour_var = "ndensity", bins = bins,
            show.legend = c(fill = TRUE, alpha = FALSE, size = FALSE)
          ) +
          ggplot2::scale_fill_viridis_c(option = "C") +
          ggplot2::scale_alpha(range = c(0.1, 0.9), guide = "none") +
          ggplot2::geom_point(
            ggplot2::aes(size = .data$value),
            color = point_color, alpha = point_alpha,
            show.legend = c(size = TRUE)
          ) +
          ggplot2::scale_size_continuous(
            name = if (use_normalized) "Normalized occupancy" else "Expected burials",
            guide = ggplot2::guide_legend(
              override.aes = list(shape = 19, colour = point_color, fill = "white", alpha = 1, stroke = 0)
            )
          ) +
          ggplot2::labs(
            title = "Horizon Bucket: {closest_state}",
            x = "Longitude", y = "Latitude", fill = "Density level"
          ) +
          ggplot2::guides(
            alpha = "none",
            fill  = ggplot2::guide_colorbar(order = 1)
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            legend.key           = ggplot2::element_rect(fill = "white", colour = NA),
            legend.background    = ggplot2::element_rect(fill = "white", colour = NA),
            legend.box.background= ggplot2::element_rect(fill = "white", colour = NA)
          ) +
          gganimate::transition_states(as.factor(.data$combined_bucket),
                                       transition_length = 15, state_length = 30) +
          gganimate::ease_aes("circular-in-out")
        
        px_w <- as.integer(width * dpi)
        px_h <- as.integer(height * dpi)
        anim <- gganimate::animate(
          p_anim,
          renderer = gganimate::gifski_renderer(),
          fps = anim_fps, duration = anim_duration,
          width = px_w, height = px_h, res = dpi, bg = "white"
        )
        anim_path <- file.path(save_dir, anim_filename)
        gganimate::anim_save(anim_path, animation = anim)
        if (isTRUE(show_plots)) print(p_anim)
      } else {
        warning("[ocp_spatial_slices] Not enough points to animate density; no frames saved.")
      }
    }
  }
  
  invisible(list(plots = base_plots, combined = combined_plot, animation_file = anim_path))
}

# ===============================
# END MODULE
# ===============================