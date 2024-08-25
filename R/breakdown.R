# Break down the item meta data with a data frame format into a list format
# This function was developed by referring to the output format of "as.irt.pars()"
# function in the plink (Weeks, 2010) R package.
breakdown <- function(x) {
  # extract the id information
  id <- x$id

  # extract the score category information
  cats <- x$cats

  # extract the unique model names
  model <- x$model
  uni.mod <- sort(unique(model))

  # extract the item parameter matrix
  pars <- data.matrix(x[, 4:ncol(x)])

  # classify the items into each of the models
  item <- purrr::map(.x = uni.mod, ~ {
    which(model == .x)
  })
  names(item) <- uni.mod

  # total number of items
  nitem <- nrow(x)

  # create a list by combining all extracted information
  list(
    id = id, cats = cats, pars = pars, model = model,
    item = item, nitem = nitem
  )
}
