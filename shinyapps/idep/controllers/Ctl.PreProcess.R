library('R6')
library(shiny)
library(shinyBS)
library(plotly)
source('businessLogic/LogicManager.R')

Ctl.PreProcess <- R6Class("Ctl.PreProcess")
LogicManager <- Logic.Manager$new()

Ctl.PreProcess$set("public", "PreProcessResult",
	function(input, session, storeVariableList){
		withProgress(message="Reading and pre-processing", {

			# Pre checking
			if(is.null(storeVariableList$RawData)){
				#if no raw data, then return null
				return(NULL)
			}



			if(input$dataFileFormat == 1){
				incProgress(1/3, "Pre-processing counts data")
				minCount = input$numMinCounts
				minSample = input$numNMinSamplesInCountCase
				logStart = input$numCountsLogStart
				incProgress(1/2,"transforming raw counts")
			}else{
				# actrually should be input.dataFileFormat == 2
				# However, when input.dataFileFormat == 3, minCound and minSample won't use.
				# So we simply it into 'if else' statement.
				incProgress(1/3,"Pre-processing data")
				minCount = input$numMinFPKM
				minSample = input$numNMinSampleInFPKMCase
				logStart = input$numFPKMLogStart
			}

			#if preprocess is not trigger, do it.
			PreProcessResult <- 
				LogicManager$PreProcessing$RawDataPreprocess(
					storeVariableList$RawData, 
					input$selectMissingValueImputationMethod,
					input$dataFileFormat,
					minCount,
					minSample,
					as.numeric(input$selectCountsTransform),
					logStart,
					input$isApplyLogTransFPKM,
					input$isNoFDR			## This parm is got in load data tab.
				)
			incProgress(1, "Done.")
			return(PreProcessResult)
		})
	}
)

 
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

Ctl.PreProcess$set("public", "InitSelectOrgUI",
	function(){
		# Create a list for Select Input options
		orgInfo <- LogicManager$DB$OrgInfo

		speciesChoice <- setNames(as.list( orgInfo$id ), orgInfo$name2 )
		# add a defult element to list    # new element name       value
		speciesChoice <- append( setNames( "NEW","**NEW SPECIES**"), speciesChoice  )
		speciesChoice <- append( setNames( "BestMatch","Best matching species"), speciesChoice  )

		# move one element to the 2nd place
		move2 <- function(i) c(speciesChoice[1:2],speciesChoice[i],speciesChoice[-c(1,2,i)])
		i= which( names(speciesChoice) == "Glycine max"); speciesChoice <- move2(i)
		i= which( names(speciesChoice) =="Zea mays"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) =="Arabidopsis thaliana"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Saccharomyces cerevisiae"); speciesChoice <- move2(i)
		i= which(names(speciesChoice)  == "Caenorhabditis elegans"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) =="Zebrafish" ); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Cow" ); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Rat" ); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Mouse"); speciesChoice <- move2(i)
		i= which(names(speciesChoice) == "Human"); speciesChoice <- move2(i)

		return(speciesChoice)
	}
)