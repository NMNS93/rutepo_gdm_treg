# Helper functions for VDJ analysis

# Functions ----

#' Plot relative abundance per patient
plt_homeo <- function(cl, x) {
  pats <- names(cl) %>%
    purrr::keep(~ stringr::str_detect(.x, paste0("[0-9]+_", x)))
  clonalHomeostasis(cl[pats],
                    cloneCall = "gene",
                    cloneTypes = c(
                      Rare = 1e-04,
                      Small = 0.001,
                      Medium = 0.01,
                      Large = 0.1,
                      Hyperexpanded = 1
                    )
  ) +
    labs(title = x) + coord_flip()
}