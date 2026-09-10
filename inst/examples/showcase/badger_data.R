# Loading and reshaping the badger data. Shared by badger_friendly.R and
# badger_power.R, which differ only in how they build the MODEL.
#
# None of this is about the package: it reads the reference's CSVs and converts
# two conventions.
#   * the reference is individual-major (X[i, t]); ET is time-major (X[t, i])
#   * the reference encodes states sparsely (S=0, E=3, I=1, D=9); ET uses the
#     position in the state space, so S=1, E=2, I=3, D=4

BADGER_DATA_DIR <- Sys.getenv("BADGER_DATA_DIR", unset = file.path(
  path.expand("~"), ".julia", "dev", "EpidemicTrajectories", "badger_ref", "RData2"))

BROCK_CHANGEPOINT_RAW <- 80   # what the raw TestMat columns encode
BROCK_CHANGEPOINT     <- 101  # what we fix it to

# All of this is user code: reading the CSVs and converting the reference's
# conventions into ET's. Two differ and are converted here:
#   * the reference is individual-major (X[i, t]); ET is time-major (X[t, i])
#   * the reference encodes states sparsely (S=0, E=3, I=1, D=9); ET uses the
#     position in the state space, so S=1, E=2, I=3, D=4
load_badger_data <- function(dir = BADGER_DATA_DIR,
                             brock_changepoint = BROCK_CHANGEPOINT,
                             known_sex_only = TRUE) {
  rd <- function(f) as.matrix(utils::read.csv(file.path(dir, f)))

  dims <- utils::read.csv(file.path(dir, "dimensions.csv"))
  m <- dims$m[1]; maxt <- dims$maxt[1]
  n_groups <- dims$G[1]; n_tests <- dims$numTests[1]
  n_seasons <- dims$numSeasons[1]; n_nu_times <- dims$numNuTimes[1]

  Xinit_raw   <- rd("Xinit.csv")
  test_mat    <- rd("TestMat.csv")
  capt_hist   <- rd("CaptHist.csv")
  capt_effort <- rd("CaptEffort.csv")
  birth_time  <- as.vector(rd("birthTimes.csv"))
  start_p     <- as.vector(rd("startSamplingPeriod.csv"))
  end_p       <- as.vector(rd("endSamplingPeriod.csv"))
  nu_times    <- as.vector(rd("nuTimes.csv"))
  sex_raw     <- as.vector(rd("sex.csv"))
  K <- utils::read.csv(file.path(dir, "Kay.csv"))$K[1]
  k <- utils::read.csv(file.path(dir, "k.csv"))$k[1]

  # The reference filters to badgers of known sex; the base model has no sex
  # effects, but keeping the filter makes the dataset comparable with the sex
  # models.
  keep <- if (known_sex_only) which(sex_raw != 0) else seq_len(m)
  old_to_new <- integer(m)
  old_to_new[keep] <- seq_along(keep)

  Xinit_raw  <- Xinit_raw[keep, , drop = FALSE]
  capt_hist  <- capt_hist[keep, , drop = FALSE]
  birth_time <- birth_time[keep]; start_p <- start_p[keep]; end_p <- end_p[keep]

  kept_rows <- old_to_new[test_mat[, 2]] != 0
  test_mat  <- test_mat[kept_rows, , drop = FALSE]
  test_mat[, 2] <- old_to_new[test_mat[, 2]]
  m <- length(keep)

  ref_code <- c("0" = 1L, "3" = 2L, "1" = 3L, "9" = 4L)
  X_init <- matrix(1L, maxt, m)
  for (i in seq_len(m)) {
    codes <- Xinit_raw[i, ]
    X_init[, i] <- ifelse(codes == -10, 1L, ref_code[as.character(codes)])
  }

  # Group membership genuinely varies with t (badgers move); 0 means "absent".
  social_group <- matrix(0L, m, maxt)
  for (i in seq_len(m)) {
    rows <- which(test_mat[, 2] == i)
    if (!length(rows)) next
    times_i  <- test_mat[rows, 1]
    groups_i <- test_mat[rows, 3]
    g <- as.integer(groups_i[which.min(times_i)])
    for (t in max(1, birth_time[i]):maxt) {
      hit <- match(t, times_i)
      if (!is.na(hit)) g <- as.integer(groups_i[hit])
      social_group[i, t] <- g
    }
  }

  age <- matrix(-10L, m, maxt)
  for (i in seq_len(m)) {
    tt <- max(1, birth_time[i]):maxt
    age[i, tt] <- as.integer(tt - birth_time[i])
  }

  # THE REPEAT-TEST FIX. `n_neg` / `n_pos` record how many ROWS read 0 / 1 for
  # each (t, i, j), so a cell with repeat captures in one quarter contributes
  # each reading once -- as both reference implementations do. `tests` keeps the
  # last row's value only for the presence check.
  tests <- array(-1L, dim = c(maxt, m, n_tests))
  n_neg <- array(0L,  dim = c(maxt, m, n_tests))
  n_pos <- array(0L,  dim = c(maxt, m, n_tests))
  for (r in seq_len(nrow(test_mat))) {
    t <- test_mat[r, 1]; i <- test_mat[r, 2]
    if (t < 1 || t > maxt || i < 1 || i > m) next
    for (j in seq_len(n_tests)) {
      v <- test_mat[r, 3 + j]
      if (is.na(v) || v == -10) next
      tests[t, i, j] <- as.integer(v)
      if (v == 1) n_pos[t, i, j] <- n_pos[t, i, j] + 1L
      else        n_neg[t, i, j] <- n_neg[t, i, j] + 1L
    }
  }
  # The Brock changepoint. The raw columns are already correct at xi = 80, and
  # the reference only ever SWAPS the two Brock columns for tests lying between
  # the current and proposed changepoint -- so fixing it at 101 means applying
  # that same swap once, for the tests in [80, 101).
  if (brock_changepoint != BROCK_CHANGEPOINT_RAW) {
    win <- if (brock_changepoint > BROCK_CHANGEPOINT_RAW)
      BROCK_CHANGEPOINT_RAW:(brock_changepoint - 1)
    else brock_changepoint:(BROCK_CHANGEPOINT_RAW - 1)
    win <- win[win >= 1 & win <= maxt]
    swap <- function(A) { tmp <- A[win, , 1]; A[win, , 1] <- A[win, , 2]
                          A[win, , 2] <- tmp; A }
    tests <- swap(tests); n_neg <- swap(n_neg); n_pos <- swap(n_pos)
  }

  capture <- t(capt_hist)
  storage.mode(capture) <- "integer"

  last_capture_time <- vapply(seq_len(m), function(i) {
    w <- which(capture[, i] == 1L); if (length(w)) max(w) else 0L
  }, integer(1))
  first_capture_time <- vapply(seq_len(m), function(i) {
    w <- which(capture[, i] == 1L)
    if (length(w)) min(w) else as.integer(max(birth_time[i], 1))
  }, integer(1))

  cam <- rd("capturesAfterMonit.csv")
  for (r in seq_len(nrow(cam))) {
    oid <- cam[r, 1]
    if (oid < 1 || oid > length(old_to_new) || old_to_new[oid] == 0) next
    nid <- old_to_new[oid]
    last_capture_time[nid] <- max(last_capture_time[nid], cam[r, 2])
  }

  season <- integer(maxt); season[1] <- 1L
  for (t in 2:maxt) season[t] <- if (season[t - 1] < n_seasons) season[t - 1] + 1L else 1L

  list(n_individuals = m, n_timepoints = maxt, n_groups = n_groups,
       n_tests = n_tests, n_seasons = n_seasons, n_nu_times = n_nu_times,
       X_init = X_init, social_group = social_group, age = age,
       capture = capture, capt_effort = matrix(as.integer(capt_effort),
                                               nrow = nrow(capt_effort)),
       tests = tests, n_neg = n_neg, n_pos = n_pos,
       sampling_period = cbind(as.integer(start_p), as.integer(end_p)),
       birth_time = as.integer(birth_time),
       last_capture_time = last_capture_time,
       first_capture_time = first_capture_time,
       season = season, nu_times = as.integer(nu_times),
       K = as.numeric(K), k = as.integer(k))
}

raw <- load_badger_data()
cat(sprintf("Badger data: %d individuals x %d quarters, %d groups, %d tests\n",
            raw$n_individuals, raw$n_timepoints, raw$n_groups, raw$n_tests))

