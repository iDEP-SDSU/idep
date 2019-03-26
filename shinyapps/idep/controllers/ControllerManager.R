source('businessLogic/LogicManager.R')
LogicManager <- Logic.Manager$new()

source('controllers/Ctl.LoadData.R')
source('controllers/Ctl.PreProcess.R')
source('controllers/Ctl.Heatmap.R')


LoadDataCtrl <- Ctl.LoadData$new()
PreProcessCtrl <- Ctl.PreProcess$new()
HeatmapCtrl <- Ctl.Heatmap$new()