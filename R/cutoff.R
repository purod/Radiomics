cutoff = makeFilter(
  name = "cutoff",
  desc = "determine the optical number of features included in the final model",
  pkg = character(0),
  supported.tasks = c("classif", "regr", "surv"),
  supported.features = c("numerics", "factors", "ordered"),
  fun = function(task, nselect, decreasing = TRUE, ...) {
    
    feature.name = getTaskFeatureNames(task)
    feature.value = c(length(feature.name):1)
    names(feature.value) = feature.name
  
    feature.value
  }

)
