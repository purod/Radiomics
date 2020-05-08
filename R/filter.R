makeFilter(
  name = "uni.cor.filter",
  desc = "Calculates scores according to alphabetical order of features",
  pkg = "mlr",
  supported.tasks = c("classif", "regr", "surv"),
  supported.features = c("numerics", "factors", "ordered"),
  fun = function(task, nselect, decreasing = TRUE, ...) {
    # univariate analysis
    #uv = generateFilterValuesData(task,imp.learner=makeLearner("surv.coxph"),method="univariate.model.score")
    #uv.list = uv$data$name[order(uv$data$univariate.model.score,decreasing = T)]
    
    uni_cox = multiple_uni_cox(time = getTaskData(task,target.extra = T)$target[,1],
                               status = ifelse(getTaskData(task, target.extra = T)$target[,2],1,0),
                               getTaskData(task,target.extra = T)$data)
    res = t(as.data.frame(uni_cox[1:841], check.names = FALSE)) 
    UniPval = res[,6]
    UniPval1 = UniPval[UniPval<0.1] #75
    storage.mode(UniPval1)="numeric"
    uv.list = names(UniPval1)[order(UniPval1,decreasing = F)]
    
    sort.data = getTaskData(task,features = uv.list, target.extra = T)$data
    sort.size = nrow(sort.data)
    # keep one feature among correlated features
    FT.name = colnames(RemoveCor(sort.data,0.8)) 
    select.feature = FT.name
    # output  
    FT.value = c(length(select.feature):1,rep(0,getTaskNFeats(task)-length(select.feature)))
    names(FT.value) = c(select.feature,setdiff(getTaskFeatureNames(task),select.feature))
    print("I am here.")
    FT.value
  }
  
)
