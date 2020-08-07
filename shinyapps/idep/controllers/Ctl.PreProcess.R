library('R6')
library(shiny)
library(shinyBS)
library(plotly)
library(processx)


Ctl.PreProcess <- R6Class("Ctl.PreProcess")

 
Ctl.PreProcess$set("public", "RawSampleInfoPreprocess",
	function(rawDesign, geneNames){
		return(
			LogicManager$PreProcessing$RawSampleInfoPreprocess(rawDesign, geneNames)
		)
	}
)


# Preprocessing tab, first plot
Ctl.PreProcess$set("public", "GetTotalReadCountsPlot",
	function(rawCount){
		# check data 
		if(is.null(rawCount)){
			return(NULL)
		}

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(rawCount)) )

		# fetch the plot
		return(LogicManager$Display$GetReadcountBarPlot(rawCount, groups))
	}
)


# Preprocessing tab, second plot
Ctl.PreProcess$set("public", "GetTransformedDataBoxPlot",
	function(transformedData){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(transformedData)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		transformedData <- as.data.frame(transformedData)

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(transformedData)) )

		# fetch the plot
		return(LogicManager$Display$GetTransformedDataBoxPlot(transformedData, groups))
	}
)

# Preprocessing tab, third plot
Ctl.PreProcess$set("public", "GetTransformedDataDensityPlot",
	function(transformedData){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(transformedData)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		transformedData <- as.data.frame(transformedData)

		# calculate group
		groups = as.factor( LogicManager$PreProcessing$DetectGroups(colnames(transformedData)) )

		# fetch the plot
		return(LogicManager$Display$GetTransformedDataDensityPlot(transformedData, groups))
	}
)

# Preprocessing tab, fourth plot
Ctl.PreProcess$set("public", "GetTransformedDataScatterPlot",
	function(transformedData){
		#
		# this function should be called after LogicManager$PreProcessing$RawDataPreprocess() function has been called. 
		#

		if(is.null(transformedData)){
			# if we cannot get result, then exit with null
			return(NULL)
		} 

		# Get transfered data
		transformedData <- as.data.frame(transformedData)

		# fetch the plot
		return(LogicManager$Display$GetTransformedDataScatterPlot(transformedData))
	}
)

Ctl.PreProcess$set("public", "GetGuessSpeciesResult",
	function(convertedIDResult){
		if(is.null(convertedIDResult)){
			return(as.data.frame("ID not recognized."))
		} 
		else { 
			return(convertedIDResult$speciesMatched)
		}
	}
)


Ctl.PreProcess$set("public", "ConvertedIDResult",
	function(GeneIDs, selectedOrg){
		if(is.null(GeneIDs)){
			return(null)
		}

		withProgress(
			message = "Convert Gene IDs", 
			detail="Converting gene IDs", 
			{
				convertedId <- LogicManager$PreProcessing$GetConvertID(GeneIDs, selectedOrg)
				incProgress(1, detail = paste("Done"))
			}
		)
		return(convertedId)
	}
)

Ctl.PreProcess$set("public", "GetAllGeneInfomation",
	function(Reactive_ConvertedIDResult, input){
		withProgress(message="Looking up gene annotation", {
			if(is.null(Reactive_ConvertedIDResult)){
				return(NULL)
			}

			ensemblIDs <- Reactive_ConvertedIDResult$ensemblIDs
			species <- Reactive_ConvertedIDResult$species
			selectOrg <- input$selectOrg

			if(is.null(ensemblIDs) || is.null(species) || is.null(selectOrg)){
				return(NULL)
			}

			return( LogicManager$PreProcessing$GetGenesInfomationByEnsemblIDs(ensemblIDs, species, selectOrg) )
		})
	}
)


Ctl.PreProcess$set("public", "InitChoiceSelectOrgUI",
	function(){
		return(LogicManager$DB$GetSpeciesChoice())
	}
)

# refer to convertedData() in 0.81 code
# Convert Transformed data based on Convert ID result. 
# If no conversion applied on ID, then use transformed data directly
Ctl.PreProcess$set("public", "GetConvertedTransformedData",
	function(input, Reactive_ConvertedIDResult, Reactive_PreProcessResult){
		withProgress(message="Converting data ... ", {
			# if no preprocess result, return null
			if(is.null(Reactive_PreProcessResult)){
				return(NULL)
			}

			# if no converted id result
			# or 'no id conversion' is selected
			# then use transformed data directly
			if(is.null(Reactive_ConvertedIDResult)){
				return(Reactive_PreProcessResult$dat)
			}

			if(input$isNoIDConversion){
				return(Reactive_PreProcessResult$dat)
			}

			transformedData <- Reactive_PreProcessResult$dat
			conversionTable <- Reactive_ConvertedIDResult$conversionTable
			incProgress(1, "Done.")
		})
		return(LogicManager$PreProcessing$ApplyConvertIDToGivenData(transformedData, conversionTable))
	}
)

# refer to convertedCounts() in 0.81 code
# Convert raw readcount data based on Convert ID result. 
# If no conversion applied on ID, then use raw readcount data directly
Ctl.PreProcess$set("public", "GetConvertedRawReadcountData",
	function(input, Reactive_ConvertedIDResult, Reactive_PreProcessResult){
		withProgress(message="Converting data ... ", {
			# if no preprocess result, return null
			if(is.null(Reactive_PreProcessResult)){
				return(NULL)
			}

			# if no converted id result
			# or 'no id conversion' is selected
			# then use raw read count data directly
			if(is.null(Reactive_ConvertedIDResult)){
				return(Reactive_PreProcessResult$rawCount)
			}

			if(input$isNoIDConversion){
				return(Reactive_PreProcessResult$rawCount)
			}

			rawReadCount <- Reactive_PreProcessResult$rawCount
			conversionTable <- Reactive_ConvertedIDResult$conversionTable
			incProgress(1, "Done.")
		})
		return(LogicManager$PreProcessing$ApplyConvertIDToGivenData(rawReadCount, conversionTable))
	}
)

# refer to ConvertedPvals() in 0.81 code
# Convert pvals data based on Convert ID result. 
# If no conversion applied on ID, then use raw pvals directly
Ctl.PreProcess$set("public", "GetConvertedPvals",
	function(input, Reactive_ConvertedIDResult, Reactive_PreProcessResult){
		withProgress(message="Converting data ... ", {
			# if no preprocess result, return null
			if(is.null(Reactive_PreProcessResult)){
				return(NULL)
			}

			# if no converted id result
			# or 'no id conversion' is selected
			# then use pvals data directly
			if(is.null(Reactive_ConvertedIDResult)){
				return(Reactive_PreProcessResult$pvals)
			}

			if(input$isNoIDConversion){
				return(Reactive_PreProcessResult$pvals)
			}

			pvals <- Reactive_PreProcessResult$pvals
			conversionTable <- Reactive_ConvertedIDResult$conversionTable
			incProgress(1, "Done.")
		})
		return(LogicManager$PreProcessing$ApplyConvertIDToGivenData(pvals, conversionTable))
	}
)

# GetSingleGenePlot
Ctl.PreProcess$set("public", "GetSingleGenePlot",
	function(input, Reactive_ConvertedTransformedData, Reactive_AllGeneInfo){

		mdf <- LogicManager$PreProcessing$GenerateDataForSingleGenePlot(
			Reactive_ConvertedTransformedData, Reactive_AllGeneInfo, input$selectOrg, 
			input$txt_PreProcess_SearchedGeneID)

		if(input$is_PreProcess_ShowIndividualSamples == 1){
			ymax <- max(mdf$value)
			p <- LogicManager$Display$GetBarPlotSingleGeneOfIndividualSamples(mdf, ymax)
		}else{
			summarizedData <- LogicManager$PreProcessing$GenerateDataForAllSamplesSingleGenePlot(mdf)
			isUseSD <- ifelse(input$is_PreProcess_useSD == 1, TRUE, FALSE)
			p <- LogicManager$Display$GetBarPlotSingleGeneOfAllSamples(summarizedData, isUseSD)
		}

		return(p)
	}
)


Ctl.PreProcess$set("public", "SaveSingleGenesPlotEpsInTempFile",
	function(input, Reactive_ConvertedTransformedData, Reactive_AllGeneInfo, file){
		cairo_ps(file, width = 8, height = 6, points = 8 )
		self$GetSingleGenePlot(input, Reactive_ConvertedTransformedData, Reactive_AllGeneInfo)
		dev.off()
	}
)


Ctl.PreProcess$set("public", "GetDataTableOfConvetedTransformedData",
	function(input, Reactive_AllGeneInfo, Reactive_ConvertedTransformedData){
		if(input$selectOrg == "NEW"| ncol(Reactive_AllGeneInfo) == 1){
			return( round(Reactive_ConvertedTransformedData, 2) )
		}else{
			tb <- merge( 
				Reactive_AllGeneInfo[,c('ensembl_gene_id','symbol')], 
				round(Reactive_ConvertedTransformedData,2),
				by.x="ensembl_gene_id", 
				by.y ="row.names", 
				all.y=T
			)

			return(tb)
		}
	}
)

Ctl.PreProcess$set("public", "SaveConvetedTransformedDataInTempFile",
	function(input, file, Reactive_AllGeneInfo, Reactive_ConvertedTransformedData, Reactive_ConvertedIDResult){
		
		withProgress(message = "Download Processed Data",
		{
			incProgress(1/5, "Prepare data ... ")

			if(input$selectOrg == "NEW" | ncol(Reactive_AllGeneInfo) == 1 ){
				# sometimes users upload unknow species but not choosing "NEW".
				result <- Reactive_ConvertedTransformedData
			}else{
				result <- LogicManager$PreProcessing$FormatProcessedTransformedDataForDownload(Reactive_AllGeneInfo, Reactive_ConvertedTransformedData, Reactive_ConvertedIDResult)
			}
			
			incProgress(4/5, "Prepare download file ... ")

			write.csv(result, file, row.names=FALSE)

			incProgress(1, "Done")
		})
	}
)

Ctl.PreProcess$set("public", "SaveConvetedReadCountDataInTempFile",
	function(input, file, Reactive_AllGeneInfo, Reactive_ConvertedRawReadcountData, Reactive_ConvertedIDResult){

		withProgress(message = "Download Processed Count Data",
		{
			incProgress(1/5, "Prepare data ... ")

			if(input$selectOrg == "NEW" | ncol(Reactive_AllGeneInfo) == 1){
				result <- Reactive_ConvertedRawReadcountData
			}else{
				result <- LogicManager$PreProcessing$FormatProcessedRawReadcountDataForDownload(Reactive_AllGeneInfo, Reactive_ConvertedRawReadcountData, Reactive_ConvertedIDResult)
			}

			incProgress(4/5, "Prepare download file ... ")

			write.csv(result, file, row.names=FALSE)

			incProgress(1, "Done")
		})
	}
)



### This function need:
### 1. orca support: https://github.com/plotly/orca
### 2. zip command exist

Ctl.PreProcess$set("public", "SaveAllPlotsInTempFile",
	function(fn, p1, p2, p3, p4){
		if( is.null(p1) || is.null(p2) || is.null(p3) || is.null(p4) ){
			return(NULL)
		}

		# orca cannot export file/to/path directly
		# 1. use orca generate file
		# 2. move to the given directry
		withProgress(message = "Download High Resolution Plots",
		{
			incProgress(0/8, "Render plot 1 of 4 ... ")
			orca(p1, format = "svg", file = '1.svg')
			incProgress(2/8, "Render plot 2 of 4 ... ")
			orca(p2, format = "svg", file = '2.svg')
			incProgress(4/8, "Render plot 3 of 4 ... ")
			orca(p3, format = "svg", file = '3.svg')
			incProgress(6/8, "Render plot 4 of 4 ... ")
			orca(p4, format = "svg", file = '4.svg')

			incProgress(7/8, "Prepare plot zip files ... ")
			zip('tmp.zip', c('1.svg','2.svg','3.svg','4.svg'))
			file.copy('tmp.zip', fn)
			
			file.remove(c('1.svg','2.svg','3.svg','4.svg','tmp.zip'))
			incProgress(8/8, "Done")
		})
	}
)

