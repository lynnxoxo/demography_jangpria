# Load all required libraries once
library(tidyverse)
library(ggpattern)
library(rlang)
library(paletteer)  # For color palettes used in plotting
library(patchwork)  
library(sf)  
library(fuzzyjoin)
library(gdata)


####################
# Helper Functions
####################

# Helper function for renaming columns to ensure consistent naming
rename_columns <- function(data, col_mapping) {
  # Remove mappings where the new name is NA
  valid_mappings <- col_mapping[!is.na(col_mapping)]
  
  if (length(valid_mappings) < length(col_mapping)) {
    warning("Column mapping contains NA values. These columns will not be renamed.")
  }
  
  # Check if old names exist in data
  missing_cols <- setdiff(valid_mappings, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns to rename do not exist in the data:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Correctly create a named vector with new_name = old_name
  # setNames(values, names)
  rename_map <- setNames(valid_mappings, names(valid_mappings))
  
  # Perform renaming using dplyr::rename explicitly to avoid masking
  data %>% dplyr::rename(!!!rename_map)
}

# Helper function to ensure numeric conversion
ensure_numeric <- function(data, cols) {
  existing_cols <- intersect(cols, names(data))
  
  if (length(existing_cols) < length(cols)) {
    warning("Some columns specified for numeric conversion do not exist in the data.")
  }
  
  data %>%
    mutate(across(all_of(existing_cols), ~ as.numeric(.)))
}


####################
# Core Functions
####################

# Harmonization Table Creation Function
create_harmonization_table <- function(relative_chron_data,
                                       system_name_col,
                                       relative_chron_col,
                                       horizon_start_col = NULL,
                                       horizon_end_col = NULL,
                                       fade_in_start_col = NULL,
                                       fade_out_end_col = NULL) {
  # Select mandatory and optional columns
  cols_to_select <- c(system_name_col, relative_chron_col)
  optional_cols <- c(horizon_start_col, horizon_end_col, fade_in_start_col, fade_out_end_col)
  optional_cols <- optional_cols[!sapply(optional_cols, is.null)]
  cols_to_select <- c(cols_to_select, optional_cols)
  
  harmonization_table <- relative_chron_data %>%
    select(all_of(cols_to_select))
  
  # Calculate period_duration if both horizon_start_col and horizon_end_col are provided
  if (!is.null(horizon_start_col) & !is.null(horizon_end_col)) {
    harmonization_table <- harmonization_table %>%
      mutate(period_duration = .data[[horizon_end_col]] - .data[[horizon_start_col]])
  } else {
    harmonization_table <- harmonization_table %>%
      mutate(period_duration = NA_real_)
  }
  
  return(harmonization_table)
}

plot_harmonization_table <- function(harmonization_table,
                                    system_name_col,
                                    relative_chron_col,
                                    horizon_start_col = NULL,
                                    horizon_end_col = NULL,
                                    fade_in_start_col = NULL,
                                    fade_out_end_col = NULL) {
  
  # Create a mapping from new names to existing names
  rename_map <- c(
    system_name = system_name_col,
    relative_chron = relative_chron_col
  )
  
  if (!is.null(horizon_start_col)) {
    rename_map["horizon_start"] <- horizon_start_col
  }
  if (!is.null(horizon_end_col)) {
    rename_map["horizon_end"] <- horizon_end_col
  }
  if (!is.null(fade_in_start_col)) {
    rename_map["fade_in_start"] <- fade_in_start_col
  }
  if (!is.null(fade_out_end_col)) {
    rename_map["fade_out_end"] <- fade_out_end_col
  }
  
  # Rename columns correctly using the helper function
  harmonization_table <- rename_columns(harmonization_table, rename_map)
  
  # Ensure numeric for existing columns
  numeric_cols <- c("horizon_start", "horizon_end", "fade_in_start", "fade_out_end")
  numeric_cols <- intersect(numeric_cols, colnames(harmonization_table))
  harmonization_table <- ensure_numeric(harmonization_table, numeric_cols)
  
  # Arrange by earliest horizon_start
  if ("horizon_start" %in% colnames(harmonization_table)) {
    harmonization_table <- harmonization_table %>%
      group_by(system_name) %>%
      mutate(earliest_start = min(horizon_start, na.rm = TRUE)) %>%
      ungroup() %>%
      arrange(earliest_start, system_name) %>%
      select(-earliest_start)
  }
  
  # Assign system_numeric for plotting
  harmonization_table <- harmonization_table %>%
    mutate(system_numeric = as.numeric(factor(system_name, levels = rev(unique(system_name)))))
  
  # Initialize ggplot
  p <- ggplot() +
    theme_minimal() +
    labs(#title = "Harmonization Table Visualization Sorted by Earliest Horizon Start Date",
         x = "Global Continuous Horizon",
         y = "Systems") +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "none",
          panel.border = element_rect(colour = "grey", fill = NA, size = 1))
  

  
  # Add geom_segment for horizon_start to horizon_end
  if ("horizon_start" %in% colnames(harmonization_table) && 
      "horizon_end" %in% colnames(harmonization_table)) {
    p <- p + geom_segment(
      data = harmonization_table %>% filter(!is.na(horizon_start) & !is.na(horizon_end)),
      aes(y = system_numeric, yend = system_numeric,
          x = horizon_start, xend = horizon_end),
      size = 5, 
      color = "lightgrey"
    )
  }

  
  # Adjust scale_y_continuous to display system names
  p <- p + scale_y_continuous(
    breaks = unique(harmonization_table$system_numeric),
    labels = unique(harmonization_table$system_name)
  )
  
  # Add vertical separators between periods within each system
  if ("horizon_start" %in% colnames(harmonization_table) && 
      "horizon_end" %in% colnames(harmonization_table)) {
    
    # Create a separator data frame
    separator_data <- harmonization_table %>%
      group_by(system_name, system_numeric) %>%
      arrange(horizon_start) %>%
      mutate(next_horizon_start = lead(horizon_start)) %>%
      filter(!is.na(next_horizon_start)) %>%
      select(system_numeric, separator_x = horizon_end) %>%
      ungroup()
    
    # Define the height of the separators (slightly larger than geom_segment)
    separator_height <- 0.2  # Adjust as needed based on system_numeric scale
    
    # Add geom_segment for vertical separators
    p <- p + geom_segment(
      data = separator_data,
      aes(x = separator_x, xend = separator_x,
          y = system_numeric - separator_height,
          yend = system_numeric + separator_height),
      size = 0.5, color = "darkgrey"
    )
  }
  
  # Add geom_errorbar for fade_in_start
  if ("fade_in_start" %in% colnames(harmonization_table) && 
      "horizon_start" %in% colnames(harmonization_table)) {
    p <- p + geom_errorbar(
      data = harmonization_table %>% filter(!is.na(fade_in_start)),
      aes(x = fade_in_start, xmin = fade_in_start,
          xmax = ifelse(!is.na(horizon_start), horizon_start, horizon_end),
          y = system_numeric),
      width = 0.15, color = "black"
    )
  }
  
  # Add geom_errorbar for fade_out_end
  if ("fade_out_end" %in% colnames(harmonization_table) && 
      "horizon_end" %in% colnames(harmonization_table)) {
    p <- p + geom_errorbar(
      data = harmonization_table %>% filter(!is.na(fade_out_end)),
      aes(x = fade_out_end,
          xmin = ifelse(!is.na(horizon_end), horizon_end, horizon_start),
          xmax = fade_out_end,
          y = system_numeric),
      width = 0.15, color = "black"
    )
  }    
  
  # Add geom_text for relative_chron
  if ("relative_chron" %in% colnames(harmonization_table) && 
      "horizon_start" %in% colnames(harmonization_table) && 
      "horizon_end" %in% colnames(harmonization_table)) {
    p <- p + geom_text(
      data = harmonization_table,
      aes(y = system_numeric, x = (horizon_start + horizon_end) / 2, label = relative_chron),
      vjust = -0.5, size = 3, color = "black"
    )
  }
  
  # Print and return the plot
  print(p)
  
  return(p)
}


# Individual Burial Data Preparation Function with Unique id_burial
prepare_mortuary_data_indiv <- function(jangdat_indiv, 
                                        harmonization_table = harmonization_table,
                                        site_name = "site_name", 
                                        system_name = "system_name", 
                                        id_burial_col = "id_burial",  
                                        period_col = "period", 
                                        age_col = "age", 
                                        sex_gender_col = "sex_gender") {
  
  # 
  # Step 1: Trim Whitespace from Specified Columns
  # 
  
  # Define the columns to trim
  columns_to_trim <- c(site_name, system_name, period_col, age_col, sex_gender_col)
  
  # Trim whitespace in the specified columns
  jangdat_cleaned <- jangdat_indiv %>%
    mutate(across(all_of(columns_to_trim), ~ str_trim(as.character(.)))) 
  
  #
  # Step 2: Add 'i' Prefix to Burial IDs
  # 
  
  # Ensure the burial ID column is of type character
  jangdat_cleaned <- jangdat_cleaned %>%
    mutate(
      !!id_burial_col := paste0("i", substr(.data[[site_name]], 1, 4), as.character(.data[[id_burial_col]]))
    )
  
  # 
  # Step 3: Return the Cleaned Dataframe
  # 
  
  return(jangdat_cleaned)
}

# Aggregated Burial Data Preparation Function
prepare_mortuary_data_aggreg <- function(data,
                                         harmonization_table,
                                         site_name = "site_name",
                                         system_name = "system_name",
                                         period = "period",
                                         amount = "amount",
                                         sex_gender = "sex_gender",
                                         age = "age") {
  
  # Step 1: Rename columns for consistency using the helper function rename_columns
  col_mapping <- c(
    site_name = site_name,
    system_name = system_name,
    period = period,
    amount = amount,
    sex_gender = sex_gender,
    age = age
  )
  
  
  # Use rename_columns to rename columns in the data
  data <- data %>%
    rename_columns(col_mapping)
  
  # Debug: Check the column names after renaming
  message("[DEBUG] Column names after renaming:")
  print(names(data))
  
  
  # Step 2: Handle missing values in the 'amount' column

  data <- data %>%
    filter(!is.na(amount))
  

  if (!"amount" %in% colnames(data)) {
    stop("The 'amount' column is missing after renaming. Please check the column renaming step.")
  }
  
  # Step 3: Expand burial counts into individual rows based on 'amount'
  data_expanded <- data %>%
    uncount(amount)  # Use tidyr::uncount to expand based on count in 'amount' column
  
  # Step 4: Assign unique IDs for each burial entry
  data_final <- data_expanded %>%
    mutate(id_burial = paste0("agg", substr(.data[[site_name]], 1, 4), row_number()))
  

  # Debug: Print the first few rows of the prepared dataset
  message("[DEBUG] Prepared aggregated data after harmonization:")
  print(head(data_final))
  
  return(data_final)
}

# Load necessary libraries
library(dplyr)
library(stringdist)
library(tidyr)
library(purrr)


#' Harmonize Chronologies #
#'
#' This function harmonizes individual mortuary data by aligning it with a harmonization table
#' based on both period and system_name. It ensures that each record is accurately mapped
#' to standardized chronologies, preserving data integrity and handling unmatched records appropriately.
#'
#' @param mortuary_data1 A dataframe containing the primary set of mortuary records.
#' @param mortuary_data2 An optional dataframe containing additional mortuary records to be combined with mortuary_data1.
#' @param harmonization_table A dataframe containing standardized mappings for period and system_name.
#' @param site_name_col The column name in mortuary_data representing cemetery names.
#' @param period_col The column name in mortuary_data representing period information (e.g., horizon_bucket).
#' @param horizon_start_col The column name in harmonization_table representing the start of the horizon.
#' @param horizon_end_col The column name in harmonization_table representing the end of the horizon.
#' @param harmonized_system_name_col The column name in harmonization_table for standardized system names.
#' @param harmonized_relative_chron_col The column name in harmonization_table for standardized relative chronologies.
#' @param horizon_bracket_size The size of each horizon bracket. Default is 10.
#'
#' @return A dataframe with harmonized chronologies, ready for further analysis.
#' @export

harmonize_chronologies <- function(mortuary_data1,
                                   mortuary_data2 = NULL,
                                   harmonization_table,
                                   site_name_col,
                                   period_col,
                                   horizon_start_col,
                                   horizon_end_col,
                                   harmonized_system_name_col,
                                   harmonized_relative_chron_col,
                                   horizon_bracket_size = 10) {
  
  # 
  # Step 1: Rename Columns in Harmonization Table
  # 
  
  # Create a named vector for renaming
  rename_map_harm <- c(
    harmonized_system_name_col,
    harmonized_relative_chron_col,
    horizon_start_col,
    horizon_end_col
  )
  names(rename_map_harm) <- c("system_name", "relative_chron", "horizon_start", "horizon_end")
  
  # Rename harmonization table columns
  harmonization_table <- harmonization_table %>%
    rename_columns(rename_map_harm) %>%
    ensure_numeric(c("horizon_start", "horizon_end"))
  
  #
  # Step 2: Rename Columns in Mortuary Data1
  #
  
  rename_map_mortuary <- c(
    site_name_col,
    period_col
  )
  names(rename_map_mortuary) <- c("site_name", "period")
  
  # Rename mortuary_data1 columns
  mortuary_data1 <- mortuary_data1 %>%
    rename_columns(rename_map_mortuary)
  
  # Check if mortuary_data2 is provided
  if(!is.null(mortuary_data2)) {
    # Identify columns in mortuary_data1 and mortuary_data2
    cols_data1 <- names(mortuary_data1)
    cols_data2 <- names(mortuary_data2)
    
    # Check if mortuary_data2 has all columns from mortuary_data1
    missing_in_2 <- setdiff(cols_data1, cols_data2)
    if(length(missing_in_2) > 0) {
      stop("mortuary_data2 is missing the following columns required by mortuary_data1: ",
           paste(missing_in_2, collapse = ", "))
    }
    
    # mortuary_data2 may have extra columns not in mortuary_data1, so we subset
    # Also reorder mortuary_data2 columns to match mortuary_data1 exactly
    mortuary_data2 <- mortuary_data2[, cols_data1, drop = FALSE]
  }  
  
  # 
  # Step 3: Combine Mortuary Data if Mortuary Data2 is Provided
  # 
  
  if (!is.null(mortuary_data2)) {
    
    # Rename columns in mortuary_data2 for consistency
    mortuary_data2 <- mortuary_data2 %>%
      rename_columns(rename_map_mortuary)
    
    # Combine the two datasets by appending them
    combined_data <- bind_rows(mortuary_data1, mortuary_data2)
    
    # Assign combined data to mortuary_data1 to continue with harmonization
    mortuary_data1 <- combined_data
    
    # Debug: Print the first few rows of combined data
    message("[DEBUG] Combined mortuary_data1 and mortuary_data2:")
    print(head(mortuary_data1))
  }
  
  # 
  # Step 4: Ensure Unique Combinations in Harmonization Table
  #
  
  # Check for duplicate combinations of system_name and relative_chron
  duplicates_harm <- harmonization_table %>%
    group_by(system_name, relative_chron) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if(nrow(duplicates_harm) > 0){
    warning("There are duplicate combinations of system_name and relative_chron in harmonization_table. These will be resolved by keeping the first occurrence.")
  }
  
  # Make harmonization_table_unique by distinct combinations of system_name and relative_chron
  harmonization_table_unique <- harmonization_table %>%
    distinct(system_name, relative_chron, .keep_all = TRUE)
  
  # 
  # Step 5: Merge Mortuary Data with Harmonization Table
  # 
  
  harmonized_data <- mortuary_data1 %>%
    regex_inner_join(
      harmonization_table_unique %>% select(system_name, relative_chron, horizon_start, horizon_end),
      by = c("system_name" = "system_name", "period" = "relative_chron"), ignore_case=TRUE
    )
  
  harmonized_data <- harmonized_data %>%
    regex_inner_join(
      jangdat_meta %>% select(site_name, relative_chronology_system, cemetery_start_ia, cemetery_end_ia),
      by = c("site_name" = "site_name"), ignore_case=TRUE
    )
  
  harmonized_data <- harmonized_data %>%
    mutate(site_name = coalesce(site_name.x, site_name.y)) %>%
    select(-site_name.x, -site_name.y)
  
  harmonized_data <- harmonized_data %>%
    mutate(system_name = coalesce(system_name.x, system_name.y)) %>%
    select(-system_name.x, -system_name.y)
  
  
  # 
  # Step 6: Check for Duplicates Based on 'id_burial'
  # 
  
  if("id_burial" %in% names(harmonized_data)){
    duplicates <- harmonized_data %>%
      group_by(id_burial) %>%
      filter(n() > 1) %>%
      ungroup()
    
    if (nrow(duplicates) > 0) {
      warning("There are still duplicate entries after ensuring unique 'relative_chron' and 'system_name'. Consider reviewing the harmonization data.")
    }
  }


  # 
  # Step 7: Compute Horizon Brackets
  # 
  
  #Takes a look at the lowest and highest occuring horizon values and 
  #multiplies them by the provided bracketing size
  
  harmonized_data <- harmonized_data %>%
    mutate(
      horizon_bracket_start = floor(horizon_start / horizon_bracket_size) * horizon_bracket_size,
      horizon_bracket_end = floor(horizon_end / horizon_bracket_size) * horizon_bracket_size
    )
  
  # 
  # Step 8: Generate Horizon Buckets
  # 
  
  harmonized_data <- harmonized_data %>%
    mutate(
      horizon_buckets = map2(
        horizon_bracket_start,
        horizon_bracket_end,
        ~ if(!is.na(.x) && !is.na(.y) && .y >= .x) {
          seq(.x, .y, by = horizon_bracket_size)
        } else {
          numeric(0) # Return an empty numeric vector instead of NULL
        }
      )
    )
  
  # 
  # Step 9: Assign Horizon Bucket by Sampling
  # 
  
  #set.seed(42)  # For testing only, otherwise set seed outside of function (e.g. in the Rmd)
  
  harmonized_data <- harmonized_data %>%
    mutate(horizon_buckets = map(horizon_buckets, ~ if(is.null(.x)) numeric(0) else .x)) %>%
    mutate(
      num_buckets = map_int(horizon_buckets, length),
      horizon_bucket = map2_dbl(
        horizon_buckets,
        num_buckets,
        ~ if (.y > 0) resample(.x, size = 1) else NA_real_
      ),
      horizon_bucket_start = horizon_bucket,
      horizon_bucket_end = horizon_bucket + horizon_bracket_size - 1
    )
  
  harmonized_data <- harmonized_data %>%
    group_by(site_name) %>%
    mutate(
      is_indiff = (period == "indiff"),
      # Gather all horizon_buckets from known (non-indiff) burials in this cemetery
      existing_buckets = list(unique(horizon_bucket[!is_indiff & !is.na(horizon_bucket)])),
      num_existing_buckets = length(existing_buckets[[1]]),

      

      # For indiff entries, if we have existing buckets, pick one randomly
      horizon_bucket = ifelse(
        is_indiff & num_existing_buckets > 0,
        resample(existing_buckets[[1]], 1),
        horizon_bucket
      ),

      # Set num_buckets for indiff entries to total number of these buckets
      num_buckets = ifelse(
        is_indiff & !is.na(horizon_bucket),
        num_existing_buckets,
        num_buckets
      ),

      # Update horizon_bucket_start and horizon_bucket_end for indiff entries now assigned a bucket
      horizon_bucket_start = ifelse(
        is_indiff & !is.na(horizon_bucket),
        horizon_bucket,
        horizon_bucket_start
      ),
      horizon_bucket_end = ifelse(
        is_indiff & !is.na(horizon_bucket),
        horizon_bucket + horizon_bracket_size - 1,
        horizon_bucket_end
      )
    ) %>%
    ungroup() %>%
    select(-is_indiff, -existing_buckets, -num_existing_buckets,
           -horizon_bracket_start, -horizon_bracket_end) 
  
  harmonized_data <- harmonized_data %>%
    group_by(site_name) %>%
    mutate(
      is_still_indiff = (period == "indiff" & is.na(horizon_bucket)),
      # For cemeteries with no assigned buckets left for indiff, use cemetery_start_ia & cemetery_end_ia
      # Filter harmonization_table for the system used by this cemetery:
      # Doing this in a rowwise manner (one row per cemetery) by summarizing first.
    ) %>%
    ungroup()
  
  cem_meta_ranges <- harmonized_data %>%
    filter(is_still_indiff) %>%  # Only consider cemeteries with still indiff burials
    distinct(site_name, relative_chronology_system, cemetery_start_ia, cemetery_end_ia) %>%
    group_by(site_name, relative_chronology_system) %>%
    summarise(
      cem_start_val = floor(mean(cemetery_start_ia, na.rm=TRUE) / horizon_bracket_size)*horizon_bracket_size,
      cem_end_val = floor(mean(cemetery_end_ia, na.rm=TRUE) / horizon_bracket_size)*horizon_bracket_size,
      .groups = "drop"
    )
  
  get_meta_buckets <- function(sys_name, start_val, end_val, horizon_bracket_size, harmonization_table) {
    if(is.na(start_val) || is.na(end_val) || start_val > end_val) {
      return(numeric(0))
    }
    
    # Filter harmonization_table for that system
    relevant_periods <- harmonization_table %>%
      filter(system_name == sys_name)
    
    # Identify horizon intervals overlapping [start_val, end_val]
    # Overlap means horizon_end >= start_val and horizon_start <= end_val
    overlapping <- relevant_periods %>%
      filter(horizon_end >= start_val & horizon_start <= end_val)
    
    if(nrow(overlapping) == 0) {
      return(numeric(0))
    }
    
    # For each overlapping row, generate buckets and combine them
    all_buckets <- overlapping %>%
      mutate(
        bucket_start = floor(horizon_start / horizon_bracket_size)*horizon_bracket_size,
        bucket_end = floor(horizon_end / horizon_bracket_size)*horizon_bracket_size
      ) %>%
      rowwise() %>%
      mutate(
        these_buckets = if (!is.na(bucket_start) && !is.na(bucket_end) && bucket_end >= bucket_start) {
          seq(bucket_start, bucket_end, by = horizon_bracket_size)
        } else {
          numeric(0)
        }
      ) %>%
      ungroup() %>%
      pull(these_buckets)
    
    # all_buckets is now a list of numeric vectors, flatten it:
    meta_buckets <- unique(unlist(all_buckets))
    
    # Filter meta_buckets to only those within [start_val, end_val]
    meta_buckets <- meta_buckets[meta_buckets >= start_val & meta_buckets <= end_val]
    
    meta_buckets
  }
  
  #Apply this logic
  cem_meta_ranges <- cem_meta_ranges %>%
    rowwise() %>%
    mutate(
      meta_buckets = list(get_meta_buckets(relative_chronology_system, cem_start_val, cem_end_val, horizon_bracket_size, harmonization_table)),
      num_meta_buckets = length(meta_buckets)
    ) %>%
    ungroup()
  
  # Join meta_buckets info back to harmonized_data
  harmonized_data <- harmonized_data %>%
    left_join(cem_meta_ranges %>% select(site_name, meta_buckets, num_meta_buckets), by = "site_name")
  
  # Assign horizon_bucket to still indiff graves if num_meta_buckets > 0
  harmonized_data <- harmonized_data %>%
    group_by(site_name) %>%
    mutate(
      # We only assign if still no horizon_bucket assigned and is_still_indiff
      # We use a round-robin assignment just like before
      # first identify those rows:
      indiff_seq = ifelse(is_still_indiff, row_number(), NA_integer_),
      horizon_bucket = ifelse(
        is_still_indiff & num_meta_buckets > 0,
        {
          chosen_bucket <- meta_buckets[[1]][ ((indiff_seq-1) %% num_meta_buckets) + 1 ]
          chosen_bucket
        },
        horizon_bucket
      ),
      num_buckets = ifelse(is_still_indiff & !is.na(horizon_bucket), num_meta_buckets, num_buckets),
      horizon_bucket_start = ifelse(is_still_indiff & !is.na(horizon_bucket), horizon_bucket, horizon_bucket_start),
      horizon_bucket_end = ifelse(is_still_indiff & !is.na(horizon_bucket), horizon_bucket + horizon_bracket_size - 1, horizon_bucket_end)
    ) %>%
    ungroup() %>%
    select(-is_still_indiff, -meta_buckets, -num_meta_buckets, -indiff_seq)

  harmonized_data <- harmonized_data %>%
    mutate(probability = ifelse(num_buckets > 0, (1 / num_buckets) * 100, NA_real_))  
  
  # Return final result
  return(harmonized_data)

  
  # 
  # Step 10: Filter Out Rows with Missing Horizon Values
  # 
  
  initial_row_count <- nrow(harmonized_data)
  debug_filter <- harmonized_data %>%
    filter(is.na(horizon_bucket_start) & is.na(horizon_bucket_end)) 
    print(debug_filter)
  harmonized_data <- harmonized_data %>%
    filter(!is.na(horizon_bucket_start) & !is.na(horizon_bucket_end))
  rows_dropped <- initial_row_count - nrow(harmonized_data)
  message("[DEBUG] Number of rows dropped due to NA horizon values: ", rows_dropped)
  
}

analyze_occupancy <- function(harmonized_data, 
                              meta_data,
                              site_name_col = "site_name", 
                              horizon_bucket_col = "horizon_bucket",
                              cemetery_size_col = "cemetery_size",
                              age_col = "age",
                              sex_col = "sex_gender",
                              cemetery_region_col = "cemetery_region",
                              coord_x_col = "coord_x",
                              coord_y_col = "coord_y",
                              fuzzy_match = TRUE, 
                              max_fuzzy_distance = 0.2, #for fuzzy matching on table join
                              change_threshold = 2,
                              distance_threshold = 50000, #for density plotting in m
                              horizon_bracket_size = 10,
                              beast_min_length = 5 #This did not end up getting used but is left in for compatibility
) {
  # Load required libraries
  library(dplyr)
  library(stringdist)
  library(ggplot2)
  library(viridis)
  library(sf)
  library(changepoint)
  library(Rbeast)
  library(dbscan)
  library(leaflet)
  library(spdep)
  library(purrr)
  library(tidyr)
  library(gganimate) # for animation
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(scales)
  
  #
  # Step 1: Standardize Site Names and Merge Meta Data
  #
  
  combined_data <- harmonized_data %>%
    left_join(meta_data, by = site_name_col)
  
  # Handle unmatched site_names with fuzzy matching if needed
  unmatched <- combined_data %>%
    filter(is.na(.data[[coord_x_col]]) | is.na(.data[[coord_y_col]])) %>%
    select(all_of(site_name_col)) %>%
    distinct()
  
  if(nrow(unmatched) > 0){
    warning("Unmatched site_names detected: ", paste(unmatched[[site_name_col]], collapse = ", "))
    
    if(fuzzy_match){
      message("[INFO] Performing fuzzy matching for unmatched site_names.")
      meta_site_names <- meta_data[[site_name_col]]
      
      find_closest_match <- function(site, meta_sites, max_dist) {
        distances <- stringdist(site, meta_sites, method = "jw")
        min_dist <- min(distances, na.rm = TRUE)
        closest_site <- meta_sites[which.min(distances)]
        if(min_dist < max_dist) closest_site else NA_character_
      }
      
      matched_names <- sapply(unmatched[[site_name_col]], find_closest_match, 
                              meta_sites = meta_site_names, max_dist = max_fuzzy_distance)
      
      harmonized_data <- harmonized_data %>%
        mutate(!!site_name_col := ifelse(
          .data[[site_name_col]] %in% unmatched[[site_name_col]] & !is.na(matched_names[.data[[site_name_col]]]),
          matched_names[.data[[site_name_col]]],
          .data[[site_name_col]]
        ))
      
      combined_data <- harmonized_data %>%
        left_join(meta_data, by = site_name_col)
      
      unmatched_after <- combined_data %>%
        filter(is.na(.data[[coord_x_col]]) | is.na(.data[[coord_y_col]])) %>%
        select(all_of(site_name_col)) %>%
        distinct()
      
      if(nrow(unmatched_after) > 0){
        warning("After fuzzy matching, remaining unmatched site_names: ", 
                paste(unmatched_after[[site_name_col]], collapse = ", "))
      } else {
        message("[INFO] All site_names matched after fuzzy matching.")
      }
    }
  }
  
  #
  # Step 2: Compute Entries Per Period and Variance from Probability (p=1/num_buckets)
  #
  
  ##### Probability per burial: p = 1/num_buckets if num_buckets>0, else p=0.
  ##### expected_entry for a bucket = sum of p for all burials in that bucket
  ##### variance for a bucket = sum over burials of p*(1-p)
  
  occupancy_summary <- combined_data %>%
    group_by(!!sym(site_name_col), !!sym(horizon_bucket_col)) %>%
    summarise(
      entries_per_period = sum(ifelse(num_buckets > 0, 1/num_buckets, 0), na.rm = TRUE),
      variance_entries = sum(ifelse(num_buckets > 0, (1/num_buckets)*(1 - (1/num_buckets)), 0), na.rm = TRUE),
      .groups = "drop"
    )
  
  #
  # Step 3: Temporal Trends (Normalized)
  #
  
  ##### We currently have entries_per_period as sum of p over all burials. We must also incorporate cemetery_size.
  ##### Each burial belongs to one cemetery; if we want normalized occupancy = (sum p)/cemetery_size?
  ##### Solution: re-join cemetery_size into occupancy_summary:
  
  occupancy_summary <- occupancy_summary %>%
    left_join(
      combined_data %>%
        distinct(!!sym(site_name_col), !!sym(cemetery_size_col)),
      by = site_name_col
    ) %>%
    mutate(
      normalized_entries = entries_per_period / .data[[cemetery_size_col]]
    )
  
  #####Calculate mean trend across all cemeteries with CI:
  ##### Combine variances: we have variance_entries for each cemetery-bucket. 

  global_trend <- occupancy_summary %>%
    group_by(!!sym(horizon_bucket_col)) %>% #grouped by horizon bucket
    summarise(
      mean_entries = mean(normalized_entries, na.rm = TRUE), #mean of entries per horizon bucket
      n_cem = n(), #n of cemeteries
      combined_variance = sum(variance_entries/(cemetery_size^2), na.rm = TRUE)/(n_cem^2), 
      se = sqrt(combined_variance), #standard deviation
      ci_lower = mean_entries - 1.96*se, #lower end of 2nd standard deviation
      ci_upper = mean_entries + 1.96*se, #upper end of 2nd standard deviation
      .groups = "drop"
    )
  
  p_global_trend <- ggplot(global_trend, aes(x = as.numeric(as.character(!!sym(horizon_bucket_col))), y = mean_entries)) +
    geom_point(color = "black", size = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "grey60") +
    geom_smooth(method = "loess", se = TRUE, span = 0.75, color = "blue", fill = "lightblue") +
    labs(
      #title = "Global Temporal Trends in Cemetery Occupancy (Normalized)",
      #subtitle = "Using probability-based entries and variance from assignment uncertainty",
      x = "Horizon Bucket",
      y = "Normalized Occupancy"
    ) +
    theme_minimal()
  
  # Region-specific trends
  region_trend <- occupancy_summary %>%
    left_join(
      combined_data %>% distinct(!!sym(site_name_col), !!sym(cemetery_region_col)),
      by = site_name_col
    ) %>%
    group_by(!!sym(cemetery_region_col), !!sym(horizon_bucket_col)) %>%
    summarise(
      mean_entries = mean(normalized_entries, na.rm = TRUE),
      n_cem = n(),
      combined_variance = sum(variance_entries/(cemetery_size^2), na.rm = TRUE)/(n_cem^2),
      se = sqrt(combined_variance),
      ci_lower = mean_entries - 1.96*se,
      ci_upper = mean_entries + 1.96*se,
      .groups = "drop"
    )
  
  p_region_trend <- ggplot(region_trend, aes(x = as.numeric(as.character(!!sym(horizon_bucket_col))), y = mean_entries)) +
    geom_point(color = "black", size = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "grey60") +
    geom_smooth(method = "loess", se = TRUE, span = 0.75, color = "darkred", fill = "orange", alpha = 2/10) +
    facet_wrap(vars(!!sym(cemetery_region_col)), scales = "free_y") +
    labs(
      #title = "Regional Temporal Trends in Cemetery Occupancy (Normalized)",
      x = "Horizon Bucket",
      y = "Normalized Occupancy"
    ) +
    theme_minimal()
  
  #
  # Step 4: Spatial Size-Density Over Time (Not Normalized)
  #
  
  ##### pick a subset of horizon buckets - for convenience hardcoded here
  
  # unique_buckets <- sort(unique(harmonized_data[[horizon_bucket_col]]))
  # if(length(unique_buckets) > 4) {
  # selected_buckets <- unique_buckets[round(seq(4, length(unique_buckets), length.out = 6))]
  # } else {
  # selected_buckets <- unique_buckets
  #}
  
  selected_buckets <- c(30, 60, 90, 120)
  
  spatial_slices <- occupancy_summary %>%
    filter(!!sym(horizon_bucket_col) %in% selected_buckets) %>%
    left_join(
      combined_data %>% distinct(site_name, coord_x = .data[[coord_x_col]], coord_y = .data[[coord_y_col]]),
      by = site_name_col
    ) %>%
    filter(!is.na(coord_x) & !is.na(coord_y))
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  
  plots_list <- list()
  
  for (hb in selected_buckets) {
    slice_data <- spatial_slices %>%
      filter(horizon_bucket == hb)
    
    p_slice <- ggplot() +
      geom_sf(data = world, fill = "gray95", color = "gray80") +
      stat_density2d(
        data = slice_data,
        aes(x = coord_x, y = coord_y, fill = after_stat(level), alpha = after_stat(level)),
        geom = "polygon",
        contour = TRUE
      ) +
      scale_fill_viridis_c(option = "C") +
      scale_alpha(range = c(0.2, 0.7), guide = "none") +
      geom_point(
        data = slice_data,
        aes(x = coord_x, y = coord_y, size = entries_per_period, show.legend = FALSE),
        color = "red3", alpha = 0.6
      ) +
      labs(
        title = paste("Horizon Bucket:", hb),
        x = "Longitude",
        y = "Latitude",
        size = "Expected Burials",
        fill = "Density Level"
      ) +
      # Restrict the map to a tight bounding box over Northern Germany and Jutland
      coord_sf(xlim = c(7, 15), ylim = c(51, 58), expand = FALSE)  +
      scale_size_continuous(guide = guide_legend(
        override.aes = list(
          fill = "white",     
          alpha = 1,         
          shape = 19         
        )
      )) +
      theme(
        legend.background = element_rect(fill = "white", color = NA),
        legend.key = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA)
      )
    
    
    plots_list[[as.character(hb)]] <- p_slice
  }
  
  #Combine the four maps into a 2x2 grid
  final_4_map_layout <- (plots_list[[1]] | plots_list[[2]]) /
    (plots_list[[3]] | plots_list[[4]])

  # --- PART B: All horizon buckets + animation
  
  all_buckets <- sort(unique(occupancy_summary[[horizon_bucket_col]]))
  
  spatial_slices_all <- occupancy_summary %>%
    filter(!!sym(horizon_bucket_col) %in% all_buckets) %>%
    left_join(
      combined_data %>% distinct(site_name, coord_x = .data[[coord_x_col]], coord_y = .data[[coord_y_col]]),
      by = site_name_col
    ) %>%
    filter(!is.na(coord_x) & !is.na(coord_y))
  
  spatial_slices_combined <- spatial_slices_all %>%
    mutate(
      combined_bucket = floor(horizon_bucket / 10) * 10
    )
  
  spatial_slices_combined <- spatial_slices_combined %>%
    group_by(horizon_bucket) %>%
    filter(n() > 1) %>%      # keep only buckets with actual data
    ungroup()
  
  spatial_slices_combined <- spatial_slices_combined %>%
    mutate(
      coord_x = as.numeric(coord_x),
      coord_y = as.numeric(coord_y),
    )


  
  # Use gganimate's transition through horizon buckets
  p_animation <- ggplot(spatial_slices_combined, aes(x = coord_x, y = coord_y)) +
    geom_sf(data = world, fill = "gray95", color = "gray80",inherit.aes = FALSE) +
    coord_sf(xlim = c(7.5, 14.5), ylim = c(51.5, 57.5)) +
    stat_density2d(
      aes(fill = after_stat(level), alpha = after_stat(level)),
      geom = "polygon",
      contour = TRUE,
      contour_var = "ndensity",
      #h = c(2, 2),
      bins = 7
    ) +
    scale_fill_viridis_c(option = "C") +
    scale_alpha(range = c(0.1, 0.9), guide = "none") +
    geom_point(aes(size = entries_per_period), color = "red3", alpha = 0.6) +
    scale_size_continuous(guide = guide_legend(
      override.aes = list(fill = "white", alpha = 1, shape = 19)
    )) +

    labs(
      title = "Horizon Bucket: {closest_state}",
      x = "Longitude",
      y = "Latitude",
      size = "Expected Burials",
      fill = "Density Level"
    ) +
    theme(
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.box.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "aliceblue")
    ) +
    transition_states(as.factor(combined_bucket), 
                      transition_length = 15, 
                      state_length = 30) +
    ease_aes("circular-in-out")
  
  # Optionally, outside this function you could do, for example:
  #   animate(p_animation, nframes = 100, fps = 5)
  #   anim_save(\"horizon_animation.gif\", animation = last_animation())
  
  return(list(
    p_global_trend       = p_global_trend,
    p_region_trend       = p_region_trend,
    final_4_map_layout   = final_4_map_layout,
    p_animation          = p_animation  # the gganimate object
  ))
}
