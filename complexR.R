#' butterfly
#'
#' This function returns the sum of two numbers.
#' @param <a> first addend
#' @param <b> second addend
#' @keywords add
#' @export
#' @examples
#' butterfly(2,3)

butterfly <- function (a,b) {
  return(a + b)
}

#' cusp
#'
#' This function returns the sum of two numbers.
#' @param <df> data frame that contains the independent (i.e., control parameters) and dependent variables, and possibly the time index.
#' @keywords cusp
#' @export
#' @examples
#' cusp(df)

cusp <- function (df) {
  # Look at the existing cusp method for consideration of other
  #  possible parameters that need to be passed.
  list() -> param_list
  # Here's the approach to use:
  # Follow Butner's 2014 piece, except:
  # a) Figure out if this is a single longitudinal case, cross-sectional
  #    data or a mixed set
  # b) Do the PDF approach to establish a prior, for both the cusp and the
  #    linear fit.
  # c) Find a posterior (via bootstrapping?) using Approximate Bayesian
  #    Computation (e.g., package:abc). Use the posterior to create
  #    confidence intervals (credible intervals?) for the fit parameters
  #    and for the R^2 as well.
  # d) Return a list with the important information: First element in
  #    the list is a list of parameters for the cusp, while the second
  #    element in the list is a list of parameters for the linear fit.
  #    Each of the two sub-lists includes estimates & standard errors for
  #    the coefficients and the R^2. Anything else?
  return(param_list)
}

#' cpt.kurtosis
#'
#' This function returns the sum of two numbers.
#' @param <a> first addend
#' @param <b> second addend
#' @keywords add
#' @export
#' @examples
#' cpt.kurtosis(2,3)

cpt.kurtosis <- function (a,b) {
  return(a + b)
}

#' cpt.skew
#'
#' This function returns the sum of two numbers.
#' @param <a> first addend
#' @param <b> second addend
#' @keywords add
#' @export
#' @examples
#' cpt.skew(2,3)

cpt.skew <- function (a,b) {
  return(a + b)
}

#' fold
#'
#' This function returns the sum of two numbers.
#' @param <a> first addend
#' @param <b> second addend
#' @keywords add
#' @export
#' @examples
#' fold(2,3)

fold <- function (a,b) {
  return(a + b)
}

#' networkEntropy
#'
#' This function returns the sum of two numbers.
#' @param <a> first addend
#' @param <b> second addend
#' @keywords add
#' @export
#' @examples
#' networkEntropy(2,3)

networkEntropy <- function (a,b) {
  return(a + b)
}

#' orbde
#'
#' This function returns the sum of two numbers.
#' @param <data_seq> sequence of (categorical) data
#' @keywords orbital decomposition
#' @export
#' @examples
#' orbde(c("A","B","B","C","B","B","C","D","D","D","E","A","E","A","D","E","A","B","B","C"))

  orbde <- function(data_seq) {
      # Returns a table  with the final orbital decomposition analysis:
      #  - String length, C
      #  - Trace of the transition matrix for the string length C, M^C
      #  - Topological entropy, H_T (this is better than H_S for these purposes;
      #    use minimum of H_T)
      #  - Shannon entropy, H_S
      #  - D_L - dimensionality via Lyapunov
      #  - chi-square
      #  - degrees of freedom, df
      #  - p-value, p
      #  - Number of possible sub-strengths of length C, N* = length + 1 - C
      #  - phi-square (phi-sq = chi-sq/N*)
      # This all follows Guastello, Chapter 21, in Guastello & Gregson, 2011,
      #  _Nonlinear Dynamical Systems Analysis for the Behavioral Sciences Using
      #  Real Data_.
      #
      # 52 codes should be sufficient for any application; it is the limit on the
      #  ORBDE software as well.
      require(data.table)  # Necessary to allow parallelizing to work
      # data.table::shift() is used here.
      require(glue)
      
      # First, recast everything into single character codes:
      unique(data_seq) -> uni_seq
      length(uni_seq) -> n_codes
      if(n_codes > 52) {
        cat("Cannot do more than 52 codes\n") # See below.
        return(NULL)
      }
      c(LETTERS, letters)[1:n_codes] -> 
        uni_rep                         # c(LETTERS, letters) has 52 elements
      uni_seq -> names(uni_rep)
      uni_rep[data_seq] -> data_seq
      
      table(data_seq) -> freqs
      
      glue_collapse(data_seq) -> coll_seq
      # Begin processing data: For C = 1
      1 -> C
      length(data_seq) -> seq_len -> N_star
      
      shift(data_seq, 1) -> shifted_seq
      length(which(data_seq == shifted_seq)) -> recurs
      
      if(recurs > 0) {
        unique(data_seq[data_seq == shifted_seq]) -> repeats
        length(repeats) -> trMC
        if (trMC > 0) {
          log2(trMC) / C -> H_T
        } else {
          -Inf -> H_T
        }
        exp(H_T) -> D_L
      } else {
        0 -> trMC
        -Inf -> H_T
        0 -> D_L
      }
      n_codes - 1 -> dof # For the singlets, this is true;
      
      table(data_seq) -> freqs  # I think this is a repeated command,
      #  but everything breaks if I remove it.
      freqs / seq_len -> p_obs_keep  # This gets used for C > 1 calculations
      rep(seq_len / n_codes,
          n_codes) -> F_exp
      freqs / seq_len -> p_obs_tab
      - 1 * sum(p_obs_tab * log(p_obs_tab)) -> H_S
      sum(freqs * log(freqs / F_exp)) * 2 ->
        chi_sq
      dchisq(chi_sq, dof) -> p
      chi_sq / N_star -> phi_sq
      
      data.frame("C" = 1,
                 "trM" = trMC,
                 "Ht" = H_T,
                 "Dl" = D_L,
                 "chi^2" = chi_sq,
                 "df" = dof,
                 "N*" = N_star,
                 "Phi^2" = phi_sq,
                 "Hs" = H_S,
                 "p" = p) -> OD_tab
      
      # Processing for 1 < C < N_star; while trace(C^M) > 0
      seq_len -> N_star
      for(len in 1:(length(data_seq) / 2)) {    # This is the ORBDE default
        # Identify all recurrences of length C
        C + 1 -> C
        N_star - 1 -> N_star
        
        # I suspect that this next part is the part that takes time
        vector(mode = "character",
               length = N_star) -> current_seq
        C -> end_index
        # The old way: routine is >10 times slower with this:
        #    for (index in 1:N_star) {
        #      paste0(data_seq[index:end_index],
        #             collapse = "") -> current_seq[index]
        #      end_index + 1 -> end_index
        #    }
        
        # New way: routine is much faster than with the old way
        for(index in 1:N_star) {
          substr(coll_seq, index, end_index) -> current_seq[index]
          end_index + 1 -> end_index
        }
        
        shift(current_seq, C) -> shifted_seq
        
        # Can we do a loop here and skip the which? 
        # E.g.: if(current_seq == shifted_seq) {...}
        which(current_seq == shifted_seq) -> repeat_nums
        # I doubt it is faster.
        
        if (length(repeat_nums) > 0) {
          length(repeat_nums) -> recurs
          
          unique(current_seq[repeat_nums]) -> repeats  # Which codes are repeated?
          length(repeats) -> trMC                   # How many repeated codes
          #  are there?
          log2(trMC) / C -> H_T                  # Topological entropy
          exp(H_T) -> D_L                        # Lyapunov dimension
          
          # The number of unique repeated codes:
          table(current_seq) -> repeated_codes -> F_obs_tab
          length(repeated_codes[repeated_codes > 1]) -> dof
          # Shannon entropy. Must do this before collapsing codes:
          sum((F_obs_tab / N_star) * 
                (log(N_star/F_obs_tab))) -> H_S
          
          F_obs_tab[F_obs_tab > 1] -> F_rep_obs_tab  # Repeated codes
          
          # Frequency expected
          length(F_rep_obs_tab) -> n_rep
          rep(1, n_rep) -> F_exp           # For all the repeats and
          
          # Is there a faster way here?
          for (i in 1:n_rep) {
            for (j in 1:C) {
              substr(names(F_rep_obs_tab)[i], j, j) -> fn
              p_obs_keep[fn] * F_exp[i] -> F_exp[i]
            }
          }
          F_exp * N_star -> F_exp
          
          # Guastello eq. 21.6: 
          #  chi-squared = 2 * sum( F_obs * ln(F_obs/F_expected) )
          sum(F_rep_obs_tab * log(F_rep_obs_tab / F_exp)) * 2 -> chi_sq
          abs(N_star - sum(F_rep_obs_tab)) -> singles
          # Need to do something for all the non-recurrent sequences, as they
          #  contribute here too.
          if(singles > 0.01 ) {        # 0.01 is capricious, but probably good
            chi_sq + 2 *singles * log(singles / (N_star - sum(F_exp))) -> chi_sq
          }
          
          # Now for p-value and phi-squared (Guastello eq 21.7):
          dchisq(chi_sq, dof) -> p
          chi_sq / N_star -> phi_sq
          
          if(!is.na(trMC)) {
            if(trMC != 0) {
              data.frame("C" = C,
                         "trM" = trMC,
                         "Ht" = H_T,
                         "Dl" = D_L,
                         "chi^2" = chi_sq,
                         "df" = dof,
                         "N*" = N_star,
                         "Phi^2" = phi_sq,
                         "Hs" = H_S,
                         "p" = p) -> OD_line
              rbind(OD_tab, OD_line) -> OD_tab
            }
          }
        } 
      }
      
      if(OD_tab$trM[1] == 0) {
        OD_tab[-1,] -> OD_tab
      }
      
      return(OD_tab)
    }

#' swallowtail
#'
#' This function returns the sum of two numbers.
#' @param <a> first addend
#' @param <b> second addend
#' @keywords add
#' @export
#' @examples
#' swallowtail(2,3)

swallowtail <- function (a,b) {
  return(a + b)
}

#' transitionNetwork
#'
#' This function returns a list containing the a transition matrix constructed from a time series and also a matrix of bootstrapped significance for the frequencies.
#' @param <ts> a vector (character or integer) containing the sequence of states, ordered temporally.
#' @param <alpha> significance level to use.
#' @keywords transition, network
#' @export
#' @examples
#' transitionNetwork(c("A", "A", "B", "B", "C", "D", "A")) -> tran_l

transitionNetwork <- function (ts, alpha = 0.05) {
  1e3 -> R           # Default number of replications
  sort(unique(ts)) -> nodes
  length(ts) -> ts_len

  if(ts_len < 2) {
    print("Time series is too short to create a transition matrix.")
    return(NULL)
  }

  if(length(nodes) < 2) {
    print("There is only one state. Transitions are irrelevant.")
    return(NULL)
  }

  if(length(nodes) > 50) {
    print("There are more than 50 distinct codes.\n
          Plots and interpretations will be messy.\n
          You have been warned.")
  }

  # Create the transition matrix:
  matrix(0,
         ncol = length(nodes),
         nrow = length(nodes),
         dimnames = list(nodes,
                         nodes)) -> trans
  for(node in 1:(ts_len - 1)) {
    trans[ts[node], ts[node + 1]] + 1 ->
      trans[ts[node], ts[node + 1]]
  }

  # Bootstrap significance of the transitions
  # There are multiple approaches to bootstrapping. We will generate
  #  all the permutations of the time series, and use those to
  #  generate the distribution of transition matrices.
  table(ts) -> obs_node_freq
  length(obs_node_freq) -> num_nodes

  max(R, 50 * 1/alpha) -> R

  # Create a large number of transition matrices
  # No seed is set here. The user will be responsible for setting
  #  the seed.
  replicate(R,
            sample(ts,
                   ts_len,
                   replace = TRUE)) -> replicates
  array(0,
        c(num_nodes, num_nodes, R),
        dimnames = list(nodes,
                        nodes,
                        1:R)) -> rep_array

  # This is probably faster in C++, especially for a large number
  #  of replications.
  for(index in 1:nrow(replicates)) {
    replicates[index,] -> ts
    # Replace this inner for-loop with sapply, or purrr:: something?
    for(node in 1:(ts_len - 1)) {
      rep_array[ts[node], ts[node + 1], index] + 1 ->
        rep_array[ts[node], ts[node + 1], index]
    }
  }

  # Create the base significance matrix
  matrix("black",
         ncol = num_nodes,
         nrow = num_nodes,
         dimnames = list(nodes,
                         nodes)) -> sig

  # Confidence intervals
  alpha -> min_quant
  1 - min_quant -> max_quant

  for (i in 1:num_nodes) {
    for (j in 1:num_nodes) {
      if (trans[i, j] < quantile(rep_array[i, j, ], min_quant)) {
        "blue" -> sig[i, j]
      } else {
        if (trans[i, j] > quantile(rep_array[i, j, ], max_quant)) {
          "gold" -> sig[i, j]
        }
      }
    }
  }

  return(list(trans, sig))
}

#' windowedEntropy
#'
#' This function returns the sum of two numbers.
#' @param <data_seq> vector of sequenced data. Data may be character or numeric.
#' @param <parts> number of partitions to use. For categorical data, will default to the number of unique codes used.
#' @param <win_max> Maximum possible length of sliding window.
#' @param <alpha> Probability to use for significance.
#' @keywords entropy, sliding window
#' @export
#' @examples
#' windowedEntropy()

windowedEntropy <- function (data_seq,
                             parts = 16,
                             lag_max = 50,
                             alpha = 0.05) {
  require(tseriesChaos)



  ami_window_size <- function(data_seq,
                              n_codes = 32,
                              lag_max = 50) {
    # The next finds the absolute minimum. We want the first local minimum.
    # Need to set arbitrary partitions for continuous data.

    # This finds the first local minimum; do we want the local or the first absolute minimum?
    #  Fraser & Swinney (1986) suggest the former.
    min(n_codes, length(unique(data_seq))) -> n_codes
    min(lag_max, length(data_seq)) -> lag_max

    1:n_codes -> nums
    unique(data_seq) -> names(nums)

    rep(0, length = length(data_seq)) -> coded_data
    for(i in 1:length(data_seq)) {
      nums[which(names(nums) == data_seq[i])] -> coded_data[i]
    }

    # For testing:
    #mutual(coded_data, partitions = n_codes,
    #       lag.max = lag_max, plot = TRUE) -> wwss

    #diff(mutual(coded_data, partitions = n_codes,
    #            lag.max = lag_max, plot = TRUE)) -> wwss
    #plot(wwss)

    if(length(which(diff(mutual(coded_data,
                                partitions = n_codes,
                                lag.max = lag_max,
                                plot = FALSE)) > 0)) > 0)
    {
      min(which(diff(mutual(coded_data,
                            partitions = n_codes,
                            lag.max = lag_max,
                            plot = FALSE)) > 0)) -> wind_sz
      return(wind_sz)
    } else {
      cat("No AMI minimum found!\n")
      return(lag_max)
    }
  }






  RJB1 <- function(data_seq, ws) {
    # Compute the windowed entropy sequence
    # Probably should vectorize this in some way.

    # The probabilities are from the entire sequence
    require(data.table)
    table(data_seq) / sum(table(data_seq)) -> probs      # Probabilities

    # Create a sequence of probabilities
    rep(0, length = length(data_seq)) -> prob_seq
    for(i in 1:length(data_seq)) {
      probs[which(names(probs) == data_seq[i])] -> prob_seq[i]
    }

    rep(0, length = length(data_seq) - ws + 1) -> ent_seq
    for(i in 1:length(ent_seq)) {
      -sum(prob_seq[i:(i+ws-1)] *
             log(prob_seq[i:(i+ws-1)], 2)) ->
        ent_seq[i]
    }
    return(ent_seq)
  }








  census <- function(data_vec, codes) {
    rep(0, length(codes)) -> census_vec
    codes -> names(census_vec)

    for(i in seq_along(data_vec)) {
      census_vec[data_vec[i]] + 1 -> census_vec[data_vec[i]]
    }
    return(census_vec)
  }

  win_ent <- function(data_seq,
                      parts = 16,
                      lag_max = 50,
                      alpha = 0.05) {
    # Follows Wiltshire, Fiore, & Butner (2017)
    # Returns a data frame (or NULL, if failed) of two time series:
    #  1. A (numeric) windowed entropy
    #  2. A marker of maxima (TRUE for a maxima, FALSE for not)
    #
    if(mode(data_seq) == "integer") {
      as.character(data_seq) -> data_seq
    }

    if(mode(data_seq) == "numeric") {
      max(data_seq) -> max_data
      min(data_seq) -> min_data
      as.character(1:parts) -> codes
      ceiling( parts *
                 (data_seq - min_data) /
                 (max_data - min_data)) -> partition
      1 -> partition[which(partition == 0)]  # partition goes 0:parts, not 1:parts
      codes[partition] -> data_seq
    }

    if(mode(data_seq) != "character") {
      return(NULL)
    }

    require(data.table)
    unique(data_seq) -> codes
    length(codes) -> n_codes

    # Because some sequences may be short
    n_codes -> parts
    min(lag_max, length(data_seq)-n_codes) -> lag_max

    table(data_seq) -> code_counts
    1:n_codes -> nums
    codes -> names(nums)

    code_counts / sum(code_counts) -> code_probs

    # The appropriate lag, according to Fraser & Swinney, 1986, is at
    #  the first minimum of mutual information.
    if(n_codes > 1) {
      which.min(mutual(unname(nums[data_seq]),
                       partitions = n_codes,       # Arbitrary partitions needed                                                     #  for continuous data; any
                       #  partition >= n_codes works
                       #  for categorical data.
                       lag.max = lag_max,
                       plot = FALSE)) -> wind_sz
    } else {
      return(NULL)
    }

    # Calculate the windowed entropy, ent[]:
    # This is different than what Wiltshire, Fiore, & Butner (2017) do. However,
    #  I think that they overestimate the issue (and they probably would think
    #  that I underestimate the issue.) They use the natural log, and only use
    #  the probabilities within the window; it is the latter that I find problematic.
    # However, there are also problems with this approach, so look at this more carefully.
    #
    # They also do some smoothing; that may be helpful
    length(data_seq) -> seq_len
    rep(0.0, length = (seq_len - wind_sz + 1)) -> ent

    # Old way
    #  for(i in 1:(length(ent))) {
    #    for(j in 0:(wind_sz-1)) {
    #      ent[i] - (code_probs[data_seq[i+j]] *
    #                  log2(code_probs[data_seq[i+j]])) -> ent[i]
    #    }
    #  }

    # A faster way, it appears, than the old way
    for(i in 1:(length(ent))) {
      ent[i] - sum(code_probs[data_seq[i:(i+wind_sz-1)]] *
                     log2(code_probs[data_seq[i:(i+wind_sz-1)]])) -> ent[i]
    }

    # Find the appropriate maxima:
    # Pass 1 - find all maxima:
    vector("logical", length(ent)) -> maxes  # All entries are FALSE by default

    if(length(ent) > 2) {
      for (i in 2:(length(ent)-1)) {
        if(ent[i] >= ent[i-1] &&
           ent[i] >= ent[i+1]) {
          TRUE -> maxes[i]
        }
      }
    } else {
      return(NULL)
    }

    # Pass 2 - which maxima indicate code probability changes from before to after a
    #  maximum located in Pass 1
    1 -> i_old
    length(data_seq) -> n_max

    for(i in 1:length(maxes)) {
      if(maxes[i] == TRUE) {
        0 -> cont_tab
        if((i+1) < (n_max - n_codes - wind_sz)) {   # Skip if too close to the end
          census(data_seq[i_old:i], codes) -> before
          census(data_seq[(i+1):n_max], codes) -> after
          rbind(before,after) -> cont_tab

          # The fisher.test() is appropriate for the next one, but it sometimes
          #  runs out of memory. The chi-squared sometimes fails for really
          #  small cell values (and zeroes); hence the is.na().
          #  The chisq.test() is also very biased, and it will generate warnings if any
          #  cell is < 5. Hence, use the fisher.test() is any cell size < 10.
          #  Probably should just use the chisq.test() if all(cont_tab > 100); will do
          #  that on the 3rd pass.
          if (sum(cont_tab) > 0) {
            if (any(cont_tab < 10)) {
              fisher.test(cont_tab) -> p_temp
            } else {
              chisq.test(cont_tab) -> p_temp
            }
            if (is.na(p_temp$p.value)) {
              FALSE -> maxes[i]
            } else {
              # The suppressWarnings seems to be necessary for about
              #  a quarter of the vector entries. I don't know why yet:
              if(p_temp$p.value > 0.10) { # 0.1 is arbitrary
                FALSE -> maxes[i]  # No change
              } else {
                i -> i_old         # Only include info from this sub-segment
              }
            }
          } else {
            FALSE -> maxes[i]
          }
        }
      }
    }
    # In principle, this refinement should be iterated until it converges. However,
    #  as we know this approach is biased, we won't try to refine something that
    #  will still be a bit biased.

    # However, let's do a 3rd pass:
    #  Fisher test for probably all, but subset by subset, with updating.

    which(maxes == TRUE) -> dividing_points
    length(dividing_points) -> num_possibles
    if(num_possibles > 0) {

      1 -> begin_a
      length(data_seq) -> end_b
      dividing_points[1] -> end_a -> begin_b
      census(data_seq[begin_a:end_a], codes) -> before
      census(data_seq[begin_b:end_b], codes) -> after
      rbind(before, after) -> cont_tab

      if(any(cont_tab < 100)) {
        fisher.test(cont_tab)$p.value -> p_value
      } else {
        chisq.test(cont_tab)$p.value -> p_value
      }

      if(p_value > alpha) {
        FALSE -> maxes[dividing_points[1]]
        1 -> dividing_points[1]
      }

      if (num_possibles > 1) {
        TRUE -> continue
        2 -> index
        while (continue) {
          dividing_points[index - 1] -> begin_a
          dividing_points[index] -> begin_b -> end_a
          census(data_seq[begin_a:end_a], codes) -> before
          census(data_seq[begin_b:end_b], codes) -> after
          rbind(before, after) -> cont_tab
          if (any(cont_tab < 100)) {
            if(length(which(colSums(cont_tab)>0))>1) {
              fisher.test(cont_tab,
                          simulate.p.value = TRUE)$p.value -> p_value
            } else {
              0 -> p_value
            }
          } else {
            chisq.test(cont_tab)$p.value -> p_value
          }
          if (p_value > alpha) {
            FALSE -> maxes[dividing_points[index]]
            dividing_points[-index] -> dividing_points
            index - 1 -> index
          }
          index + 1 -> index
          if (index > length(dividing_points)) {
            FALSE -> continue
          }
        }
      }

      #  if (length(which(maxes == TRUE)) > 1) {
      #    length(data_seq) -> end_b
      #    max(which(maxes == TRUE)) -> begin_b -> end_a
      #    if (length(which(maxes == TRUE)) > 2) {
      #      which(maxes == TRUE) -> temp_maxes
      #      temp_maxes[length(temp_maxes) - 1] -> begin_a
      #    } else {
      #      1 -> begin_a
      #    }
      #  } else {
      #    1 -> begin_a
      #    length(data_seq) -> end_b
      #    dividing_points[num_possibilities] -> begin_b -> end_a
      #  }

      #  unname(table(data_seq[begin_a:end_a])) -> before
      #  unname(table(data_seq[begin_b:end_b])) -> after
      #  rbind(before, after) -> cont_tab
      #  if (any(cont_tab < 100)) {
      #    fisher.test(cont_tab)$p.value -> p_value
      #  } else {
      #    chisq.test(cont_tab)$p.value -> p_value
      #  }
      #  if (p_value > alpha) {
      #    FALSE -> maxes[dividing_points[num_possibiles]]
      #  }
    }

    data.frame("Windowed_Entropy" = ent,    # Entropy
               "Maxima" = maxes) -> ent_df  # True if a local maximum of entropy
    return(ent_df)
  }






  return(a + b)
}
