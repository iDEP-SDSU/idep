source('businessLogic/LogicManager.R')
LogicManager <- Logic.Manager$new()

source('controllers/Ctl.ReactiveVars.R')
source('controllers/Ctl.LoadData.R')
source('controllers/Ctl.PreProcess.R')
source('controllers/Ctl.Heatmap.R')
source('controllers/Ctl.Report.R')
source('controllers/Ctl.Kmeans.R')
source('controllers/Ctl.PCA.R')


LoadDataCtrl <- Ctl.LoadData$new()
PreProcessCtrl <- Ctl.PreProcess$new()
HeatmapCtrl <- Ctl.Heatmap$new()
ReportCtrl <- Ctl.Report$new()
ReactVarsCtrl <- Ctl.ReactiveVars$new()
KmeansCtrl <- Ctl.Kmeans$new()
PCACtrl <- Ctl.PCA$new()



