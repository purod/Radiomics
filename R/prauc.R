### 1. Define a function that calculates the misclassification rate
prauc.fun = function(task, model, pred, feats, extra.args) {
# the function doesn't work for predict.type="both"
  pr.roc <- pr.curve(scores.class0 = pred$data[,match(paste0("prob.",pred$task.desc$positive),colnames(pred$data))], weights.class0 = ifelse(pred$data[,2]==pred$task.desc$positive,1,0), curve = F)
  pr.roc$auc.davis.goadrich
}
### Generate the Measure object
prauc = makeMeasure(
id = "prauc", name = "area under precision-recall curve",
properties = c("classif", "classif.multi", "req.pred",
"req.truth"),
minimize = F, best = 1, worst = 0,
fun = prauc.fun
)




