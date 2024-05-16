get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

#s = "lambda.1se" or "lambda.min"
get_coef_dataframe <- function(cv_model_fit, s = "lambda.1se")
{
  temp <- as.data.frame(as.matrix(coef(cv_model_fit, s = "lambda.1se")))
  temp$gene <- row.names(temp)
  temp <- temp[-which(temp$gene == "(Intercept)"),]
  temp$minuslog10 <- -log10(abs(temp$s1))
  temp <- temp[order(temp$minuslog10),]
  row.names(temp) <- 1:length(temp[,1])
  
  return(temp)
}