source('localization/zh.R')
source('localization/en.R')

if(CONFIG_SERVER_LANGUAGE == 'ZH'){
	TxtLibrary <- Localization.ZH 
}else{
	TxtLibrary <- Localization.EN
}

source('views/View.LoadData.R')
source('views/View.PreProcess.R')
source('views/View.Heatmap.R')
source('views/View.Kmeans.R')
source('views/View.Report.R')






