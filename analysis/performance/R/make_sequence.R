#' Generate a sequence of numbers
#'
#' This function creates a sequence of numbers from a starting value to a maximum value,
#' with an percentage increment between each number.
#'
#' @param start The starting value of the sequence.
#' @param max_val The maximum value for the sequence.
#' @param increment The percentage increment between each number in the sequence. Default is 0.10 (10%).
#' @param integer A logical value indicating whether the sequence values should be rounded to integers. Default is \code{FALSE}.
#' @return A numeric vector containing the generated sequence.
#' @examples
#' generate_sequence(1, 10)
#' generate_sequence(1, 10, increment = 0.20)
#' generate_sequence(1, 100, increment = 0.15, integer = TRUE)
#'
#'@keywords internal
#'
make_sequence <- function(start,
                              max_val,
                              increment = 0.10,
                              integer = FALSE) {
  sequence <- numeric()
  current <- start

  # Special handling for sequence starting at 0
  if (start == 0) {
    sequence <- c(sequence, 0)
    current <- increment
  }

  while (current <= max_val) {
    if (integer) {
      sequence <- c(sequence, ceiling(current))
      current <- ceiling(current) * (1 + increment)
    } else {
      sequence <- c(sequence, current)
      current <- current * (1 + increment)
    }
  }

  # # Add max_val if not already included
  # if(tail(sequence, 1) < max_val) {
  #   sequence <- c(sequence, max_val)
  # }

  return(sequence)
}
