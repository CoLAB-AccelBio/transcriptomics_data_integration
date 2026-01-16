#############################################
# Load libraries
#############################################

library(plotly)
library(dplyr)


#############################################
# Functions
#############################################

# Function to convert data frame into sunburst compatible data
as.sunburstDF <- function(DF, value_column = NULL, add_root = FALSE){
  require(data.table)
  
  colNamesDF <- names(DF)
  
  if(is.data.table(DF)){
    DT <- copy(DF)
  } else {
    DT <- data.table(DF, stringsAsFactors = FALSE)
  }
  
  if(add_root){
    DT[, root := "Total"]  
  }
  
  colNamesDT <- names(DT)
  hierarchy_columns <- setdiff(colNamesDT, value_column)
  DT[, (hierarchy_columns) := lapply(.SD, as.factor), .SDcols = hierarchy_columns]
  
  if(is.null(value_column) && add_root){
    setcolorder(DT, c("root", colNamesDF))
  } else if(!is.null(value_column) && !add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c(setdiff(colNamesDF, value_column), "values"))
  } else if(!is.null(value_column) && add_root) {
    setnames(DT, value_column, "values", skip_absent=TRUE)
    setcolorder(DT, c("root", setdiff(colNamesDF, value_column), "values"))
  }
  
  hierarchyList <- list()
  
  for(i in seq_along(hierarchy_columns)){
    current_columns <- colNamesDT[1:i]
    if(is.null(value_column)){
      currentDT <- unique(DT[, ..current_columns][, values := .N, by = current_columns], by = current_columns)
    } else {
      currentDT <- DT[, lapply(.SD, sum, na.rm = TRUE), by=current_columns, .SDcols = "values"]
    }
    setnames(currentDT, length(current_columns), "labels")
    hierarchyList[[i]] <- currentDT
  }
  
  hierarchyDT <- rbindlist(hierarchyList, use.names = TRUE, fill = TRUE)
  
  parent_columns <- setdiff(names(hierarchyDT), c("labels", "values", value_column))
  hierarchyDT[, parents := apply(.SD, 1, function(x){fifelse(all(is.na(x)), yes = NA_character_, no = paste(x[!is.na(x)], sep = ":", collapse = " - "))}), .SDcols = parent_columns]
  hierarchyDT[, ids := apply(.SD, 1, function(x){paste(x[!is.na(x)], collapse = " - ")}), .SDcols = c("parents", "labels")]
  hierarchyDT[, c(parent_columns) := NULL]
  return(hierarchyDT)
}

##### A wrapper to saveWidget which compensates for arguable BUG in saveWidget which requires `file` to be in current working directory (see post https://github.com/ramnathv/htmlwidgets/issues/299 )
saveWidgetFix <- function ( widget, file, ...) {
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  htmlwidgets::saveWidget(widget,file=file,...)
}


#############################################
# Load data
#############################################

# Read data from the Excel file
file_path <- "/Users/jacek_marzec/Library/CloudStorage/OneDrive-SharedLibraries-AssociaçãoAccelbio/AccelBio Central Hub - Projects_gastric_cancer/data/GSE84433/target_GSE84433.txt"
data <- read.table( file_path, sep="\t", as.is=TRUE, header=TRUE)
colnames(data) <- make.names(colnames(data))

# For demonstration purposes, creating a sample data frame
# 4 categories
data <- data.frame(
  Category = data$M_risk,
  Subcategory1 = data$T_stage,
  Subcategory2 = data$N_stage,
  Subcategory3 = data$Survival
)



# Summarize data to count occurrences
# 1 category
summarized_data_1_cat <- data %>%
  count(Category, name = "Value")

# 2 categories
summarized_data_2_cat <- data %>%
  count(Category, Subcategory1, name = "Value")

# 3 categories
summarized_data_3_cat <- data %>%
  count(Category, Subcategory1, Subcategory2, name = "Value")

# 4 categories
summarized_data_4_cat <- data %>%
  count(Category, Subcategory1, Subcategory2, Subcategory3, name = "Value")

# Convert data frame into sunburst compatible data
sunburstDF_1_cat <- as.sunburstDF(summarized_data_1_cat, value_column = "Value", add_root = TRUE)
sunburstDF_2_cat <- as.sunburstDF(summarized_data_2_cat, value_column = "Value", add_root = TRUE)
sunburstDF_3_cat <- as.sunburstDF(summarized_data_3_cat, value_column = "Value", add_root = TRUE)
sunburstDF_4_cat <- as.sunburstDF(summarized_data_4_cat, value_column = "Value", add_root = TRUE)


#############################################
# Create plot
#############################################

# Create the sunburst plot
sunburst_plot_1_cat <- plot_ly(data = sunburstDF_1_cat, ids = ~ids, labels= ~labels, parents = ~parents, values= ~values, type='sunburst', branchvalues = 'total')
sunburst_plot_2_cat <- plot_ly(data = sunburstDF_2_cat, ids = ~ids, labels= ~labels, parents = ~parents, values= ~values, type='sunburst', branchvalues = 'total')
sunburst_plot_3_cat <- plot_ly(data = sunburstDF_3_cat, ids = ~ids, labels= ~labels, parents = ~parents, values= ~values, type='sunburst', branchvalues = 'total')
sunburst_plot_4_cat <- plot_ly(data = sunburstDF_4_cat, ids = ~ids, labels= ~labels, parents = ~parents, values= ~values, type='sunburst', branchvalues = 'total')


##### Save interactive plot as html file
saveWidgetFix(sunburst_plot_1_cat, file = "/Users/jacek_marzec/Library/CloudStorage/OneDrive-SharedLibraries-AssociaçãoAccelbio/AccelBio Central Hub - Projects_gastric_cancer/data/GSE84433/target_GSE84433_sunburst_plot_1_cat.html")
saveWidgetFix(sunburst_plot_2_cat, file = "/Users/jacek_marzec/Library/CloudStorage/OneDrive-SharedLibraries-AssociaçãoAccelbio/AccelBio Central Hub - Projects_gastric_cancer/data/GSE84433/target_GSE84433_sunburst_plot_2_cat.html")
saveWidgetFix(sunburst_plot_3_cat, file = "/Users/jacek_marzec/Library/CloudStorage/OneDrive-SharedLibraries-AssociaçãoAccelbio/AccelBio Central Hub - Projects_gastric_cancer/data/GSE84433/target_GSE84433_sunburst_plot_3_cat.html")
saveWidgetFix(sunburst_plot_4_cat, file = "/Users/jacek_marzec/Library/CloudStorage/OneDrive-SharedLibraries-AssociaçãoAccelbio/AccelBio Central Hub - Projects_gastric_cancer/data/GSE84433/target_GSE84433_sunburst_plot_4_cat.html")
