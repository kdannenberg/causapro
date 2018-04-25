# small utilities

sink.reset <- function() {
  if (sink.number() > 0) {
    for (i in 1:sink.number()) {
      sink()
    }
  }
}

# behaves like paste, but does not add seps when an element is "",
# in other words, removes the "" before pasting
pastes <- function(..., sep = " ", collapse = NULL) {
  return(do.call("paste", c(list(...)[list(...) != ""], sep = sep, collapse = collapse)))
}

plus_minus_x <- function(vector, offset) {
  vector <- unlist(lapply(vector, function(n) return(seq(n - offset, n + offset))))
  return(sort(unique(vector)))
}

#' Remove the rows or columns in a matrix that only contain zeros
#'
#' @param matrix the matrix
#' @param rows_or_columns if set to 1, zero rows are removed, otherwise zero columns
remove_zero_rows_or_columns <- function(matrix, rows_or_columns) {
  if (rows_or_columns == 1) {
    return(matrix[rowSums(abs(matrix)) != 0, ])
  } else {
    return(matrix[,colSums(abs(matrix)) != 0])
  }
  # return(matrix[apply(matrix != 0, rows_or_columns, any), , drop = TRUE])
}

get_index <- function(position, data) {
  return(which(colnames(data) == as.character(position)))
}


#' This function is used e.g. for plotting connected components in pymol.
#' Normally, the connected components are sorted alphabetically.
#' Since in proteins the false impression that neighbouring residues are in the
#' same connected components could arise when they have subsequent and thus
#' indistinguishable rainbow colors, the components can be mixed (mix argment given).
#' Then, mix_offset determines how many other positions there are (at least) between two
#' originally neighbouring components

#' @param mix if mix = "every_mix_offset_th", the mixed version is obtained by appending to the
#' new list every i'th element from the original list, where i = mix_offset + 1,
#' continuing with incremeted offset every time the end of the list is reached
reorder_list_of_lists <- function(list, ordering, mix_mode = "mix_offset_between_original_neighbours", mix_offset = floor(sqrt(length(list))),
                                  sort_mode, sort_descending = TRUE) {
  if (ordering == "mix") {
    if (mix_offset > 0) {
      mixed_connected_components <- list()
      if (mix_mode == "every_mix_offset_th") {
        for (rest in seq(0, mix_offset)) {
          for (multipicity in seq(0, floor(length(list)/mix_offset))) {
            index_in_conn_comp <- multipicity * mix_offset + rest + 1
            if (index_in_conn_comp <= length(list)) {
              mixed_connected_components[[length(mixed_connected_components) + 1]] <- list[[index_in_conn_comp]]
              # mixed_connected_components <- append(mixed_connected_components, connected_components[[index_in_conn_comp]])
            }
          }## for(i in 1:length(list)) {
          ##              for(j in 1:length(connected_components)) {
          ##  if(list[i] %in% connected_components[[j]])
          ## names(list)[i] = colors[j]
          ##              }
          ##              }
          ##           for (i in 1:length(connected_components)) {
          ## for (j in 1:length(connected_components[[i]])) {

          ##           nAttrs$fillcolor[connected_components[[i]][j]] = colors[i]
          ## }
          ## }
        }
      } else {
        # there are at least mix_offest + 1 blcoks (+ one block with the rest of the division).
        # To obtain the permutation, first, the first elements of all the blocks are taken,
        # then the second ones and so on.
        block_size = floor(length(list) / mix_offset)
        for (offset in 0:(block_size-1)) {
          for (block in seq(0, ceiling(length(list) / block_size))) {
            index_in_conn_comp <- block * block_size + offset + 1
            if (index_in_conn_comp <= length(list)) {
              mixed_connected_components[[length(mixed_connected_components) + 1]] <- list[[index_in_conn_comp]]
              # mixed_connected_components <- append(mixed_connected_components, connected_components[[index_in_conn_comp]])
            }
          }
        }
      }
      list <- mixed_connected_components
    }
  } else if (grepl("sort", ordering)) {
    if (sort_mode == "length") {
      if (sort_descending) {
        sorted_connected_components <- list[order(vapply(list, length, 1L), decreasing = TRUE)]
      } else {
        sorted_connected_components <- list[order(vapply(list, length, 1L), decreasing = FALSE)]
      }
    }
    list <- sorted_connected_components
  }
  return(list)
}

