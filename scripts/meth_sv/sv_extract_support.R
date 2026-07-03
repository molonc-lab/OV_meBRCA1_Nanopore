require(stringr)
require(tidyverse)
require(data.table)
require(dplyr)

process_cigar_positions <- function(strand, ref_start, cigar) {
  
  # Parse the CIGAR string into operation lengths and types
  parsed <- str_match_all(cigar, "(\\d+)([MIDNSHP=X])")[[1]]
  lengths <- as.numeric(parsed[, 2])  # Operation lengths
  ops <- parsed[, 3]  # Operations
  
  # Initialize variables for reference and read positions
  ref_start_temp <- ref_start  # Temporary reference start
  ref_end <- ref_start_temp  # Reference mapped end
  read_start <- 1  # Read start position
  read_length <- 0  # Total read-consuming length

  # Reverse the CIGAR string for negative strand to make it like positive strand
  if (strand == "-") {
    lengths <- rev(lengths)
    ops <- rev(ops)
  }

  # Process the CIGAR string and calculate positions (same logic for both strands after reversal)
  for (i in seq_along(ops)) {
    # Skip soft clip for reference calculation
    if (ops[i] == "M" | ops[i] == "D" | ops[i] == "N" | ops[i] == "=" | ops[i] == "X") {
      ref_end <- ref_end + lengths[i]  # Update reference position
    }
    if (ops[i] == "M" | ops[i] == "I" | ops[i] == "=" | ops[i] == "X") {
      read_length <- read_length + lengths[i]  # Count read-consuming operations
    }
  }

  # Adjust the read start for leading soft clipping (for both strands)
  if (ops[1] == "S") {
    read_start <- read_start + lengths[1]
  }

  # Compute read end
  read_end <- read_start + read_length - 1

  # Adjust ref_end to be inclusive (we subtract 1 because reference ends are inclusive)
  ref_end <- ref_end - 1

  # Return the final positions
  return(paste(ref_start_temp, ref_end, read_start, read_end, collapse = ","))
}


deduce_reference_position <- function(cigar_vector, read_start_vector, read_position_vector, strand_vector) {
  result <- vector("numeric", length(cigar_vector))  # Preallocate result vector
  
  for (i in seq_along(cigar_vector)) {
    cigar <- cigar_vector[i]
    read_start <- read_start_vector[i]
    read_position <- read_position_vector[i]
    strand <- strand_vector[i]
    
    # Parse the CIGAR string into operations
    ops <- str_extract_all(cigar, "\\d+[MIDNSHP=X]")[[1]]
    
    if (strand == "+") {
      # For positive strand: Traverse CIGAR string sequentially
      ref_position <- read_start  # Start at the mapped reference position
      read_pos_counter <- 0      # Read position starts at 1
      
      for (op in ops) {
        length <- as.numeric(str_extract(op, "\\d+"))
        type <- str_extract(op, "[MIDNSHP=X]")
        
        if (type == "M") {
          # Match: Advances both read and reference positions
          if (read_pos_counter + length > read_position) {
            result[i] <- ref_position + (read_position - read_pos_counter)-1
            break
          }
          ref_position <- ref_position + length
          read_pos_counter <- read_pos_counter + length
        } else if (type == "I") {
          # Insertion: Advances read position only
          if (read_pos_counter + length > read_position) {
            result[i] <- NA
            break
          }
          read_pos_counter <- read_pos_counter + length
        } else if (type == "D") {
          # Deletion: Advances reference position only
          ref_position <- ref_position + length
        } else if (type == "S") {
          # Soft clipping: Advances read position only
          if (read_pos_counter + length > read_position) {
            result[i] <- NA
            break
          }
          read_pos_counter <- read_pos_counter + length
        }
      }
    } else if (strand == "-") {
      # For negative strand: Count backward from read_start
      ref_position <- read_start  # Start at the mapped reference position
      read_pos_counter <- 0       # Read position starts at 1
      
      # Reverse the CIGAR operations for backward traversal
      ops <- rev(ops)
      
      for (op in ops) {
        length <- as.numeric(str_extract(op, "\\d+"))
        type <- str_extract(op, "[MIDNSHP=X]")
        
        if (type == "M") {
          # Match: Advances both read and reference positions
          if (read_pos_counter + length > read_position) {
            result[i] <- ref_position - (read_position+1 - read_pos_counter)
            break
          }
          ref_position <- ref_position - length
          read_pos_counter <- read_pos_counter + length
        } else if (type == "I") {
          # Insertion: Advances read position only
          if (read_pos_counter + length > read_position) {
            result[i] <- NA
            break
          }
          read_pos_counter <- read_pos_counter + length
        } else if (type == "D") {
          # Deletion: Advances reference position only
          ref_position <- ref_position - length
        } else if (type == "S") {
          # Soft clipping: Advances read position only
          if (read_pos_counter + length > read_position) {
            result[i] <- NA
            break
          }
          read_pos_counter <- read_pos_counter + length
        }
      }
    }
  }
  
  return(result)
}



deduce_ref_pos_ins_shift <- function(cigar_vector, read_start_vector, read_position_vector, strand_vector, insertion_positions_vector) {
  result <- vector("numeric", length(cigar_vector))  # Preallocate result vector

  for (i in seq_along(cigar_vector)) {
    cigar <- cigar_vector[i]
    read_start <- read_start_vector[i]
    read_position <- read_position_vector[i]
    strand <- strand_vector[i]
    insertion_position <- insertion_positions_vector[i]

    # Parse the CIGAR string into operations
    ops <- str_extract_all(cigar, "\\d+[MIDNSHP=X]")[[1]]

    if (strand == "+") {
      # For positive strand: Traverse CIGAR string sequentially
      ref_position <- read_start  # Start at the mapped reference position
      read_pos_counter <- 0      # Read position starts at 1

      for (op in ops) {
        length <- as.numeric(str_extract(op, "\\d+"))
        type <- str_extract(op, "[MIDNSHP=X]")

        if (type == "M") {
          # Match: Advances both read and reference positions
          if (read_pos_counter + length >= read_position) {
            result[i] <- ref_position + (read_position - read_pos_counter) - 1
            break
          }
          ref_position <- ref_position + length
          read_pos_counter <- read_pos_counter + length
        } else if (type == "I") {
          # Insertion: Advances read position only
          if (read_pos_counter + length >= read_position) {
            # If the insertion overlaps the target position and matches the provided insertion position
            if (ref_position == insertion_position) {
              result[i] <- ref_position + (read_position - read_pos_counter) - 1
            } else {
              result[i] <- ref_position + (read_position - read_pos_counter) - 1
            }
            break
          }
          if (ref_position == insertion_position) {
            ref_position <- ref_position + length  # Increment reference position for insertion
          }
          read_pos_counter <- read_pos_counter + length
        } else if (type == "D") {
          # Deletion: Advances reference position only
          ref_position <- ref_position + length
        } else if (type == "S") {
          # Soft clipping: Advances read position only
          if (read_pos_counter + length >= read_position) {
            result[i] <- NA
            break
          }
          read_pos_counter <- read_pos_counter + length
        }
      }
    } else if (strand == "-") {
      # For negative strand: Count backward from read_start
      ref_position <- read_start  # Start at the mapped reference position
      read_pos_counter <- 0       # Read position starts at 1

      # Reverse the CIGAR operations for backward traversal
      ops <- rev(ops)

      for (op in ops) {
        length <- as.numeric(str_extract(op, "\\d+"))
        type <- str_extract(op, "[MIDNSHP=X]")

        if (type == "M") {
          # Match: Advances both read and reference positions
          if (read_pos_counter + length >= read_position) {
            result[i] <- ref_position - (read_position - read_pos_counter)
            break
          }
          ref_position <- ref_position - length
          read_pos_counter <- read_pos_counter + length
        } else if (type == "I") {
          # Insertion: Advances read position only
          if (read_pos_counter + length >= read_position) {
            # If the insertion overlaps the target position and matches the provided insertion position
            if (ref_position == insertion_position) {
              result[i] <- ref_position - (read_position - read_pos_counter)
            } else {
              result[i] <- ref_position - (read_position - read_pos_counter)
            }
            break
          }
          if (ref_position == insertion_position) {
            ref_position <- ref_position + length  # Increment reference position for insertion
          }
          read_pos_counter <- read_pos_counter + length
        } else if (type == "D") {
          # Deletion: Advances reference position only
          ref_position <- ref_position - length
        } else if (type == "S") {
          # Soft clipping: Advances read position only
          if (read_pos_counter + length >= read_position) {
            result[i] <- NA
            break
          }
          read_pos_counter <- read_pos_counter + length
        }
      }
    }
  }

  return(result)
}



