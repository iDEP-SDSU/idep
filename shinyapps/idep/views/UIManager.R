#` View Manager
#` This file contains a reference to any workable ui.
#` 
#` 
#` 
#` 
#` 
#` 
#` 
#` 
#` 

source('localization/zh.R')
source('localization/en.R')

if(CONFIG_SERVER_LANGUAGE == 'zh'){
	TxtMsg <- Localization.ZH$new() 
}else {
	TxtMsg <- Localization.EN$new()
}

source('views/View.LoadData.R')





View.Check.UI <- function(){
    
}

View.Get.UI <- function(){

}

View.Get.UI.List <- function(){

}










