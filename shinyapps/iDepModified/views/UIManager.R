source('localization/zh.R')
source('localization/en.R')

if(CONFIG_SERVER_LANGUAGE == 'ZH'){
	TxtLibrary <- Localization.ZH 
}else{
	TxtLibrary <- Localization.EN
}

source('views/View.LoadData.R')






