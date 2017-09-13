par(mfrow = c(1,3))

cons_S <- get_conservation("S", protein = "PDZ")
cons_G <- get_conservation("G", protein = "PDZ")

data <- read_data(get_data_description(protein = "PDZ", type_of_data = "DDS", subtype_of_data = ""))
vars <- apply(data, 2, var)

# TODO scale effects auch für vectoren
# vars <- scale_effects(effects = vars)
vars <- vars / max(vars)
colors_by_var <- color_by_effect(vars, int_pos)


# barplot(vars, col = colors_by_var, main = "Variance")



# cons_G <- cons_G / max(cons_G)
colors_by_cons <- color_by_effect(cons_G, int_pos)

barplot(cons_G, col = colors_by_cons, main = "Conservation DG")

 
# cons_S <- cons_S / max(cons_S)
colors_by_cons <- color_by_effect(cons_S, int_pos)

barplot(cons_S, col = colors_by_cons, main = "Conservation DS")

diff <- cons_S - cons_G
barplot(diff)