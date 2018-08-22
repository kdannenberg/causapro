library(graph)
library(dagitty)
library(pcalg)

## gets an alignment as a matrix of characters
## returns a binary matrix of the same size with the ones corresponding to
## the columnwise most frequent proteins
alignment_to_binary_matrix <- function(alignment) {
  ## n is the number of columns
  n <- dim(alignment)[2]
  maj <- function(t) {
    return(names(which.max(table(t))))
  }
  ## storing the columnwise most frequent proteins in a vector
  most_frequent_aminoacids <- apply(alignment, 2, maj)
  ## return binary matrix
  return(t(1*apply(alignment, 1, `==`, most_frequent_aminoacids)))
}

set_parameters <- function(FUN, parameters) {
  # return(function(...) return(do.call(FUN, parameters, ...)))
  # return(function(...) {return(do.call(FUN, list(parameters, list(...))))})  # previously: return(do.call(FUN, c(parameters, list(...)))) (problematic when arguemts are matrices)
  return(function(...) {return(do.call(FUN, c(parameters, list(...))))})
}

function_set_parameters <- set_parameters

# for_coloring -> output hierarchical (list with different sorts of interesting positions as vectors), otherwise one vector
# the for_coloring result of this function can be converted to the other one by conv_for_coloring
# (which is what happens internally if for_coloring is not set), so this result is more useful in general
# std-Reihenfolge: grün-gelb-rot-blau
interesting_positions <- function(protein, position_numbering = "crystal", for_coloring = FALSE, coloring = "", colors = "", counts) {
  if (!is.null(names(protein))) {
    # if ("int_pos" %in% names(protein$general)) {
    #   if (!for_coloring) {
    #     return(protein$general$int_pos)
    #   } else {
    #     return(conv_for_coloring(protein$general$int_pos))
    #     # TODO: function for converting between for_coloring and other
    #   }
    # }
    if ("caption" %in% names(protein$summary)) {
      if (grepl("bin_approx", protein$summary$caption)) {
        position_numbering <- "alignment"
      }
      protein <- strsplit(protein$summary$caption, "_")[[1]][1]
    } else if ("outpath" %in% names(protein$summary)) {
      if (grepl("bin_approx", protein$summary$outpath)) {
        position_numbering <- "alignment"
      }
      protein <- strsplit(protein$summary$outpath, "/")[[1]][2]
    }
  }
  if (is.null(coloring) || coloring == "none") {
    list <- list()
  } else if (grepl("pymol", tolower(coloring))) {
    # something like:
    # rainbow_colors <- rainbow(length(which(names(node_clustering) == "")))
    # names(node_clustering)[names(node_clustering) == ""] <- rainbow_colors
    # colors <- names(node_clustering)
  } else {
    interesting_pos <- NULL     ## list?!
    if ((protein == "pdz") || (protein == "PDZ")) {
      if (missing(position_numbering)) {
        position_numbering <- "crystal"
      }
      if (grepl("crystal", position_numbering)) {
        if (coloring == "es") {
          ligand = c(318, 322, 323, 324, 326, 327, 328, 330, 331, 332, 334, 337, 348, 352, 366, 373, 380)
          cluster = c(304, 305, 306, 309, 312, 354, 355, 357, 391, 395)
          list <- list(ligand = ligand, cluster = cluster)
          names(list) <- c("#3b73b1","#b0c7df")
        } else if (tolower(coloring) == "dds" || tolower(coloring) == "s") {
          high = c(322, 323, 325, 327, 329, 330, 336, 347, 351, 353, 359, 362, 363, 364, 372, 375, 376, 379, 386, 388)
          list <- list(high = high)
          names(list) <- c("#FFD700")
        } else {
          # numbering crystal structure
          main = c(372)
          high = c(322, 325, 329, 330, 347, 353, 362, 376, 380, 386) # interesting ones
          low = c(328, 340, 371, 385)
          rather_high = c(336, 341, 363, 365) # 352 (missing)
          # if (for_coloring && grepl("all", coloring)) {
          if (grepl("all", coloring)) {
            list <- list(main = main, high = high, low = low, rather_high = rather_high)
          } else {
            list <- list(main = main, high = high)
          }
          names(list) <- c("#69A019", "#FFD700", "#CC0000", "#FF9933")[1:length(list)]
        }
      } else if (position_numbering == "alignment") {
        # numbering alignment
        main = c(98)
        high = c(24, 28, 32, 33, 65, 72, 81, 102, 106, 119) # interesting ones
        list <- list(main = main, high = high)
        names(list) <- c("#69A019", "#FFD700", "#CC0000", "#FF9933")[1:length(list)]
      }
    } else if (protein == "GTB") {
      bind_donor <- c(c(121, 123, 126, 213, 346, 352), c(301, 302, 188, 211))       # source: first group:pdb, second group: GTAB-Paper von Friedemann
      bind_acceptor <- c(c(233, 245, 303, 326, 348), 266)  # source: first group: pdb, 266: GTA/B-Paper von Friedemann
      bind_donor_indirect_H2O <- c(124, 125, 351)
      bind_donor_indirct_Mn <- c(211, 213)
      close_to_bind_donor <- c() # plus_minus_x(bind_donor, 5)
      close_to_bind_donor <- c() # plus_minus_x(bind_donor, 5)
      list <- list("#FFD700" = bind_donor,             # yellow
                   "#69A019" = bind_acceptor,          # green
                   "#FF9933" = c(close_to_bind_donor, bind_donor_indirect_H2O, bind_donor_indirct_Mn),    # orange
                   "#6B8E23" = close_to_bind_donor)    # olive
    } else if (tolower(protein) == "pmo") {
      bind_ligand <- c(444, 443, 344, 345, 374, 442)
      others <- c(333, 347, 396, 348, 391, 392, 375)
      bile_acid <- c(486, 505, 506, 507, 508, 509)
      list <- list("#69A019" = bind_ligand,            # green
                   "#FFD700" = others,                 # yellow
                   "#CC0000" = bile_acid)              # red
    } else if (tolower(protein) == "pdi") {
        bind_ligand <- in_both_monomers(c(444, 443, 344, 345, 374, 442))
        others <-  in_both_monomers(c(347, 396, 348, 391, 392, 375))
        bile_acid <-  in_both_monomers(c(509))
        list <- list("#69A019" = bind_ligand,            # green
                     "#FFD700" = others,                 # yellow
                     "#CC0000" = bile_acid)              # red
    } else if (protein == "p38g") {
      if (coloring == "FS4") {
        blue_ <- c(53, 161, 215, 116, 137, 159, 154, 120, 130, 167, 219, 283, 125, 220, 209, 216, 212, 291, 287)
        green_left <- c(16, 26, 112, 170, 58, 337, 75, 87, 323, 78, 89, 109, 343) # left cluster in Fig. S4
        green_right <- c(23, 55, 30, 119, 33, 41, 141, 169, 45, 134)              # right cluster in Fig. S4
        yellow_left <- c(66, 182, 187, 186, 86, 77, 81, 149, 144, 174)
        yellow_right <- c(90, 357, 367, 150, 293, 361, 285, 333, 300)
        red_left <- c(197, 198, 201, 253, 268, 250, 262, 265)
        red_right <- c(225, 294, 238, 276, 239, 241, 277, 288, 292)

        list <- list("#69A019" = c(green_left, green_right),     # green
                     "#FFD700" = c(yellow_left, yellow_right),   # yellow
                     "#CC0000" = c(red_left, red_right),         # red
                     "#1874CD" = blue_)                          # blue

      } else #if (grepl("FS3", coloring)) {
        # if (grepl("FS3", coloring) && grepl("mix", coloring) && !(grepl("manual", coloring))) {
        #   # if (missing(counts)) {
        #   counts <- read.csv2("../Data/FigS3.csv", row.names = 1, check.names=FALSE, skip = 1)
        #   counts <- as.matrix(counts, ncol = 4)
        #   # }
        #   if (!is.na(as.numeric(colors))) {
        #     round_categories <- as.numeric(colors)
        #   } else {
        #     round_categories <- 1
        #   }
        #   list <- classify_nodes(counts, round_categories = round_categories, mix = TRUE, colors = colnames(counts))
        #   # fillcolor <- colors_for_nodes(clustering = node_clusters, colors = names(node_clusters))
        # } else
          if (grepl("FS3", coloring) && grepl("mix", coloring) && grepl("manual", coloring)) { # mixed manually (simple)
          # Mischungen nur, wenn mind. ein Achtel
          red_ <- c(250, 262, 265)
          red_little_blue <- c(241, 294, 276, 238) # < 3/4 rot
          red_some_blue <- c(198, 225, 239, 253, 268) # > 3/4 rot
          red_with_blue <- c(red_little_blue, red_some_blue)
          red_half_blue <- c(201, 197) # ca. 1/2 rot
          red_blue <- c(red_with_blue, red_half_blue)
          blue_ <- c(23, 55, 86, 187, 300, 333, 285, 186, 174, 144, 357, 66, 90, 150, 367, 361, 293, 292, 288, 277, 287, 220, 291, 216, 212, 45, 33, 283, 130, 209, 134, 120, 125)
          green_ <- c(26, 75, 30, 343)
          green_little_blue <- c(81, 87, 78, 16)
          green_some_blue <- c(149, 337, 77, 58, 141)
          green_with_blue <- c(green_little_blue, green_some_blue)
          green_half_blue <- c(182, 323)
          green_blue <- c(green_with_blue, green_half_blue)
          green_little_yellow <- c(89, 109, 170)
          yellow_half_green <- c(116, 112)
          yellow_some_green <- c(53, 161)
          yellow_green <- c(green_little_yellow, yellow_half_green, yellow_some_green)
          yellow_little_blue <- c(119, 169)
          yellow_some_blue <- c(159, 219, 137, 41)
          yellow_with_blue <- c(yellow_some_blue, yellow_little_blue)
          yellow_half_blue <- c(167)
          yellow_blue <- c(yellow_with_blue, yellow_half_blue)
          yellow_ <- c(154, 215)
          if (grepl("simple", coloring)) {
            yellow_ <- c(yellow_, yellow_blue, yellow_green)
            green_ <- c(green_, green_little_yellow, green_blue)
            red_ <- c(red_, red_blue)
            list <- list("#69A019" = sort(green_),       # green
                         "#FFD700" = sort(yellow_),      # yellow
                         "#CC0000" = sort(red_),         # red
                         "#1874CD" = sort(blue_))        # blue
          } else {
            list <- list("#69A019" = sort(green_),       # green
                         "#FFD700" = sort(yellow_),      # yellow
                         "#CC0000" = sort(red_),         # red
                         "#1874CD" = sort(blue_),        # blue
                         "#8B008B" = sort(red_blue),     # lilac
                         "#008B8B" = sort(green_blue),   # bluegreen
                         "#6B8E23" = sort(yellow_blue),  # olive
                         "#C0FF3E" = sort(yellow_green)) # brigthgreen
          }
      } else if (coloring == "modules") {
        module_list <- list()
        module_list[[1]] <- c(16, 26, 30, 33, 41, 53, 58, 75, 77, 78, 81, 87, 89, 109, 112, 116, 119, 137, 141, 144, 149, 161, 169, 170, 182, 215, 323, 337, 343)
        module_list[[2]] <- c(125, 130, 150, 197, 198, 201, 212, 220, 225, 238, 239, 241, 250, 253, 262, 265, 268, 276, 277, 287, 288, 291, 292, 293, 294, 300)
        module_list[[3]] <- c(16, 26, 30, 33, 41, 53, 78, 81, 89, 109, 112, 116, 119, 120, 130, 134, 137, 141, 154, 159, 161, 167, 169, 170, 215, 219, 220, 283, 343)
        module_list[[4]] <- c(26, 58, 66, 75, 77, 78, 81, 86, 87, 89, 90, 109, 141, 144, 149, 169, 170, 174, 182, 186, 187, 323, 343, 357)
        module_list[[5]] <- c(150, 174, 197, 198, 201, 225, 238, 239, 241, 250, 253, 262, 265, 268, 276, 288, 292, 293, 294)
        module_list[[6]] <- c(33, 41, 53, 55, 58, 77, 78, 81, 87, 89, 90, 109, 112, 116, 119, 130, 137, 141, 144, 149, 159, 161, 167, 169, 170, 174, 187, 201, 212, 215, 220, 283, 323, 333)
        module_list[[7]] <- c(30, 41, 58, 75, 77, 78, 81, 86, 87, 89, 109, 119, 141, 144, 149, 170, 174, 220, 323, 343, 357)
        module_list[[8]] <- c(16, 23, 26, 30, 33, 41, 45, 53, 55, 58, 66, 75, 77, 78, 81, 87, 89, 109, 112, 116, 119, 134, 137, 141, 149, 150, 159, 161, 169, 170, 283, 323, 337, 343, 357)
        module_list[[9]] <- c(58, 66, 77, 78, 81, 86, 89, 90, 144, 149, 174, 182, 186, 187, 239, 300, 333, 357, 367)
        module_list[[10]] <- c(116, 119, 120, 130, 137, 141, 144, 154, 159, 161, 167, 169, 212, 215, 219, 283, 292, 293)
        module_list[[11]] <- c(90, 150, 197, 209, 220, 238, 276, 283, 288, 292, 293, 357, 361, 367)
        module_list[[12]] <- c(45, 66, 77, 86, 90, 141, 144, 174, 182, 186, 187, 293, 323, 357, 361, 367)
        module_list[[13]] <- c(26, 58, 75, 81, 87, 149, 182, 343)
        module_list[[14]] <- c(33, 120, 125, 150, 169, 220, 277, 287, 288, 294, 323)
        module_list[[15]] <- c(137, 209, 216, 288, 291, 293)
        module_list[[16]] <- c(78, 87)
        module_list[[17]] <- c(58, 116, 141)
        module_list[[18]] <- c(89, 149)
        module_list[[19]] <- c(77, 109, 170)
        module_list[[20]] <- c(78, 323)
        module_list[[21]] <- c(141, 323)
        module_list[[22]] <- c(137, 170)
        module_list[[23]] <- c(58, 77)
        module_list[[24]] <- c(78, 174)
        module_list[[25]] <- c(141, 144, 149)
        module_list[[26]] <- c(141, 174)
        module_list[[27]] <- c(285)

        list <- module_list

        nAttrs <- list()
        nAttrs$fillcolor <- fillcolor
      } else { # different coloring for p38g # Figure S3 automatically mixed
        # return(list())
        if (missing(counts)) {
          counts <- read.csv2("Data/FigS3.csv", row.names = 1, check.names=FALSE, skip = 1)
          counts <- as.matrix(counts, ncol = 4)

          if (!is.na(as.numeric(colors)) && !is.null(colors)) {
            round_categories <- as.numeric(colors)
          } else {
            round_categories <- 1
          }

          list <- classify_nodes(counts, round_categories = round_categories, mix = TRUE, base_colors = colnames(counts))
        }
      }
    } else {
      list <- list()
    }
    if (!is.null(interesting_pos)) {  # !is.null(list) ?!
      warning("No interesting positions known")
      list <- list()
    }
  }

  if (for_coloring) {
    return(list)
    # n <- interesting_positions(protein = "PDZ")
    # list <- lapply(unique(names(n)), function(color) {return(unname(n[which(names(n) == color)]))})
    # names(list) <- unique(names(n))
  } else {
    # do not rename duplicate names
    return(conv_for_coloring(list))
    # return(unlist(list))
  }
}

in_both_monomers <- function(positions) {
  return(c(paste0(positions, "_1"), paste0(positions, "_2")))
}

## it seems much easier to convert from the for_coloring version to the other one,
## so for simplicity's sake I do it that way for now
## that means that the standardized way to pass int_pos would be the for_coloring version
conv_for_coloring <- function(int_pos) {
  # int and n are undefined in this function. I therefore copied the code from interesting_positions and
  # replaced ti by this function there
  # list <- lapply(unique(names(int)), function(color) {return(unname(int[which(names(int) == color)]))})
  # names(list) <- unique(names(n))
  # return(list)
  return(setNames(unlist(int_pos, use.names = FALSE), rep(names(int_pos), lengths(int_pos))))
}



# deprecated. Use
#
# parameters_to_info_file(parameters_for_info_file, outpath)
# # pc_func <- function(outpath) {return(pc_fun(outpath))}
# pc_func <- function_set_parameters(pc_fun, parameters = list(outpath = outpath))
# loaded_object_ok_fun <- function(pc) {return(length(pc@graph@nodes) != dim(data)[2])}
# compute_if_not_existent(filename = paste(outpath, "-pc.RData", sep = ""), FUN = pc_func, obj_name = "pc",
#                         compute_anew = compute_pc_anew, loaded_object_ok_fun = loaded_object_ok_fun)
#
# instead of get_pc(pc_fun, outpath, compute_pc_anew, parameters_for_info_file, data = data).
#
# Compute pc if necessary
get_pc <- function(pc_fun, outpath, compute_pc_anew, parameters_for_info, data) {
  parameters_to_info_file(parameters_for_info, outpath)
  if ((file.exists(paste(outpath, "-pc.RData", sep = ""))) && !(compute_pc_anew)) {
    filename <- paste(outpath, "-pc.RData", sep = "")
    load(filename)
    if (!exists("pc")) {
      warning("The file did not contain an object of name 'pc'!")
    } else if (length(pc@graph@nodes) != dim(data)[2]) {
      warning("The pc-object in the file did not correspond to the data!")
    } else {
      print(paste("pcAlgo-object loaded from ", filename, ".", sep = ""))
      print_pc_results_to_info_file(paste(outpath, sep = ""), pc)
      return(pc)
    }
    # file.copy(paste(outpath, "-info-pc.txt", sep = ""), paste(outpath, "-info.txt", sep = ""), overwrite = TRUE)
  }  else {
    warning("pc object not found.")
    warning(outpath)
    directories <- strsplit(outpath, "/")
    output_dir <- paste(directories[[1]][1:(length(directories[[1]])-1)], collapse = "/", sep = "/")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings = TRUE, recursive = TRUE, mode = "0777")
      print("Directory created.")
    }

  pc <- pc_fun(outpath)

  save(pc, file = paste(outpath, "-pc.RData", sep = ""))
  print_pc_results_to_info_file(paste(outpath, sep = ""), pc)
  }
  return(pc)
}


scale_effects <- function(effects, rank = FALSE, amplification_factor = FALSE, neg_effects = "pos") {  neg_effects = "abs"

  amplify_with_factor <- function(effects, element_that_should_be_scaled_to = 1,
                                  value_the_element_should_be_scaled_to = 0.9, cut_values_at = 1) {
    sorted_effects <- sort(effects, decreasing = TRUE)
    if (sorted_effects[element_that_should_be_scaled_to] == 0) {
      return(effects)
    }
    factor = value_the_element_should_be_scaled_to / sorted_effects[element_that_should_be_scaled_to]
    # factor <- 0.9 / sorted_effects_pos[2] # permutated_position nicht skalieren, alle anderen so, dass der Zweitgrößte bei 0.9 ist
    effects <- effects * factor
    if (is.numeric(cut_values_at)) {
      effects[,1][effects[,1] > cut_values_at] <- cut_values_at
    }
    # pos_with_colors[,2][pos_with_colors[,2] > 1] <- 1 # sollte nur eine Position sein, falls factor != 1, dann keine
    return(effects)
  }

  effects_na <- as.matrix(effects[,1][is.na(effects[,1])])
  effects <- as.matrix(effects[,1][!is.na(effects[,1])])

  if (!length(effects) == 0) {
    if (neg_effects == "discard") {
      effects[,1][effects[,1] < 0] <- 0
    } else if (neg_effects == "abs") {
      effects <- abs(effects)
    }
    if (rank) {
      if (neg_effects == "sep") {
        effects_pos <- as.matrix(effects[,1][effects[,1] >= 0])
        effects_neg <- as.matrix(effects[,1][effects[,1] < 0])
        effects_neg <- -effects_neg
        effects_pos <- cbind(apply(effects_pos, 2, rank))
        effects_neg <- cbind(apply(effects_neg, 2, rank))

        effects_pos <- effects_pos - min(effects_pos) + 1
        effects_neg <- effects_neg - min(effects_neg) + 1
        effects_pos <- effects_pos / max(effects_pos)
        effects_neg <- effects_neg / max(effects_neg)

        effects_neg <- -effects_neg
        effects <- rbind(effects_pos, effects_neg)
        # effects_ <- effects_[order(rownames(effects_)), , drop = FALSE]
        # effects <- effects_
      } else {
        effects <- cbind(apply(effects, 2, rank))

        effects <- effects - min(effects) + 1
        effects <- effects / max(effects)
        #effects <- effects/dim(effects)[1]
      }
      amplification_exponent <- 1               # will man das? I think so
    } else {
      if (neg_effects != "sep") {
        min_eff <- sort(effects)[1]
        if (min_eff < 0) {
          effects = (effects - min_eff) / 2
        }
        if (amplification_factor) {
          effects <- amplify_with_factor(effects)
        }
      } else {
        effects_pos <- as.matrix(effects[,1][effects[,1] >= 0])
        # effects_pos <- as.matrix(effects[,1][effects[,1] >= 0 & !is.na(effects[,1])])
        effects_neg <- as.matrix(effects[,1][effects[,1] < 0])
        # effects_neg <- as.matrix(effects[,1][effects[,1] < 0 & !is.na(effects[,1])])
        # effects_neg <- cbind(effects_neg[!is.na(effects_neg),]) # NAs nicht doppelt haben -- in neg rausschmeißen!
                                                                # TODO: in dan anderen Fällen ggf. auch !!!!!
        effects_neg <- -effects_neg


        if (amplification_factor) {
          effects_pos <- amplify_with_factor(effects_pos)
        }

        if (amplification_factor && dim(effects_neg)[1] > 1) {
          effects_neg <- amplify_with_factor(effects_neg)
        }

        effects_neg <- -effects_neg
        effects <- rbind(effects_pos, effects_neg)
        # effects_ <- effects_[order(rownames(effects_)), , drop = FALSE]
        # effects <- effects_
      }
    }
  }

  effects <- rbind(effects, effects_na)

  effects <- effects[order(as.numeric(rownames(effects))), , drop = FALSE]
  # effects <- effects

  return(effects)
}







