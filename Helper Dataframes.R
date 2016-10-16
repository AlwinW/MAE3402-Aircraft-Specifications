#----------------------------
#--- Dataframes
#============================

## Repeat Rows of a Dataframe ======================================================================
RepeatRows <- function(data, each = 1) {
  # Determine each data type
  if (is.numeric(each) & length(each) == 1)
    each = each
  else if (is.data.frame(each))
    each = nrow(each)
  else
    each = length(each)
  # Expand the data and return it
  data <- data[rep(row.names(data), each = each), 1:length(data)]
  rownames(data) <- NULL
  return(as.data.frame(data))
}