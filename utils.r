#setwd("~/Dropbox/USP/Doutorado/ForceDecPaper")

require(nonlinearTseries)
require(TSDecomposition)
require(tseriesChaos)
require(fNonlinear)
require(parallel)
require(wavelets)
require(foreach)
require(c3net)
require(doMC)
require(Rssa)
require(FNN)

standardize <- function(x){
  if(length(x) == 1){
    return(x)
  } else if(sd(x) == 0){
    return(x-mean(x))
  } else {
    return((x-mean(x))/(sd(x)))
  }
}
normalize <- function(values, a = -1, b = 1){
  max = max(values)
  min = min(values)
  if((max - min) == 0) return(max)
  return( a + (((values - min)*(b - a))/(max - min)) )
}
toTimeSpace = function(attractor, m, d){
  n = nrow(attractor)
  ts = attractor[,1]
  for(i in 2:m)
    ts = c(ts, attractor[(n-d+1):n,i])
  return(ts)
}
md.dist <- function(d1, d2){
  d = d1 - d2
  return(sqrt(sum(diag(d%*%t(d)))))
}
mddl <-function(obs, pred){
  return(TSDecomposition::mddl(pred, obs))
}
mda <- function(obs, pred, m, d){
  at1 = tseriesChaos::embedd(pred, m=m, d=d)
  at2 = tseriesChaos::embedd(obs,  m=m, d=d)
  n   = nrow(at1)
  d   = nrow(at2) - nrow(at1)
  if(d == 0) return(mean(md.dist(at1,at2)))
  else return(min(unlist(lapply(0:d, function(x) mean(md.dist(at1, at2[(1+x):(n+x),]))))))
}
mae.md <- function(obs, pred, m, d){
  at1 = tseriesChaos::embedd(pred, m=m, d=d)
  at2 = tseriesChaos::embedd(obs,  m=m, d=d)
  n   = nrow(at1)
  d   = nrow(at2) - nrow(at1)

  if(d == 0) return(mean(abs(at1 - at2)))
  else return(min(unlist(lapply(0:d, function(x) mean(abs(at1 - at2[(1+x):(n+x),]))))))
}
rmse.md.aux <- function(d1, d2){
  d = d1 - d2
  return(sqrt(mean(diag(d%*%t(d)))))
}
rmse.md <- function(obs, pred, m, d){
  at1 = tseriesChaos::embedd(pred, m=m, d=d)
  at2 = tseriesChaos::embedd(obs,  m=m, d=d)
  n   = nrow(at1)
  d   = nrow(at2) - nrow(at1)

  if(d == 0) return(mean(md.dist(at1, at2)))
  else return(min(unlist(lapply(0:d, function(x) rmse.md.aux(at1, at2[(1+x):(n+x),])))))
}
mae <- function(obs, pred){
  n = length(pred)
  d = length(obs) - length(pred)
  if(d == 0)return(mean(abs(obs-pred)))
  else return(min(unlist(lapply(0:d,function(x) mean(abs(obs[(1+x):(n+x)]-pred))))))
}
rmse <- function(obs, pred){
  n = length(pred)
  d = length(obs) - length(pred)
  if(d == 0) return(sqrt(mean((obs-pred)^2)))
  else return(min(unlist(lapply(0:d,function(x) sqrt(mean((obs[(1+x):(n+x)]-pred)^2))))))
}
evaluateMetrics <- function(obs, pred, m, d){
  return(c(mddl(obs, pred),
           mda(obs, pred, m, d),
           mae.md(obs, pred, m, d),
           rmse.md(obs, pred, m, d),
           mae(obs, pred),
           rmse(obs, pred)
  ))
}
evaluateResult <- function(obs, resultSeries, params, techName, testId){
  resultTable = data.frame(testId   = numeric(0), tech    = character(0),
                           paramIdx = numeric(0), param   = character(0),
                           mddl     = numeric(0), mda     = numeric(0),
                           mae_md   = numeric(0), rmse_md = numeric(0),
                           mae      = numeric(0), rmse    = numeric(0),
                           dist     = numeric(0))

  validTestIdx = which(abs(rowSums(resultSeries)) > 0)
  resultSeries = matrix(resultSeries[validTestIdx,], ncol=ncol(resultSeries))

  m = unique(params$m)
  d = unique(params$d)

  er = apply(resultSeries, 1, function(pred) evaluateMetrics(obs, pred, m, d))
  er[which(er < 10^-12)] = 0

  resultTable[1:length(validTestIdx),5:10] = t(er)
  resultTable$dist = sqrt(normalize(resultTable$mddl, 0, 1)^2 +
                          normalize(resultTable$mda, 0, 1)^2)
  resultTable$testId = testId
  resultTable$tech = techName
  resultTable$paramIdx = validTestIdx
  resultTable$param = (apply(format(params), 1, paste, collapse=","))[validTestIdx]

  return(resultTable)
}
loadSeriesFile <- function(seriesFolder){
  #select all series file in the data folder
  seriesFile = list.files(path = seriesFolder, full.names = TRUE)
  seriesList = list()
  for(i in 1:length(seriesFile))
    seriesList[[i]] = get(load(seriesFile[i]))
  return(seriesList)
}
exec <- function(F, series, param){
  result = c()
  tryCatch({
      result = F(series, param)
  }, warning = function(w){
    #write(paste('WARN:', w, '\n'), stderr())
  }, error = function(e){
    #write(paste('ERRO:', e, '\n'), stderr())
  })
  if(is.null(result)) result = rep(0, length(series))
  return(result)
}

foreachParam <- function(F, series, params){
  if(is.null(params)) return(exec(series, F, NULL))
  return(t(apply(params, 1, function(x) exec(F, series, x))))
}

gridSearch <- function(F, params, seriesList, modelFolder, techName, cores = 1){

  doMC::registerDoMC(cores)

  resultTable =	foreach::foreach(i=1:length(seriesList), .combine='rbind') %dopar% {
 	  st = Sys.time()
    cat(paste('ts:',i,'- begin\n'))

    #load series object from RData file
    seriesObj = seriesList[[i]]
    params$m = unlist(seriesObj$det.embDim)
    params$d = unlist(seriesObj$det.sepDim)

    #Grid Search for better params
    resultSeries = foreachParam(F, seriesObj$series, params)

    #evaluate results with know deterministic component
    obs      = seriesObj$det.series
    rTable   = evaluateResult(obs, resultSeries, params, techName, i)
	  bestIdx  = which.min(rTable$dist)

    #save model result into model folder
    model = list( model.name = techName, F = F,
                  best.param = params[rTable$paramIdx[bestIdx],],
                  eval = rTable[bestIdx,] )
    save(model, file=paste(modelFolder, '/', techName, '_',
                           formatC(i, width=2, format='d', flag='0') ,'.RData',sep=''))

    cat(paste('ts:',i,'- done - ', Sys.time()-st,'\n'))
    rTable
  }
  return(resultTable)
}
