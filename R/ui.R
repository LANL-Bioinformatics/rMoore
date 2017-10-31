# Load the ggplot2 package which provides
# the 'mpg' dataset.

fluidPage(
  ## hack to  align action button  with textInput box  # use 25px if textInput label != ""
  tags$head(   tags$style(type='text/css', 'button#update{ margin-top: 20px;}') ),
  titlePanel( 'Search biomedical literature in Europe PMC'),
  #helpText("Multiple accession ids must be separated by ';'"),
  fluidRow( column(10, textInput("acc_num",   value= "" , width= "100%", label= "", placeholder="accession number(s); multiple accession numbers must be separated by space" )) ,
            # color like submitButton  
            column(2, actionButton("update", "Search",  style="color: #fff; background-color: #337ab7; border-color: #2e6da4" ))
            #column(2, numericInput("obs", "Number of articles", 100 ))
            ),
  br(),
  
  tabsetPanel(type = "tabs", 
              tabPanel("Articles",    DT::dataTableOutput("table")),
              tabPanel("Plot",     dygraphs::dygraphOutput('plot') )         
  )  
)
