library(shiny)
library(DT)
library(dygraphs)
library(xts)

getData<-function(acc_num, limit) {
  library(readr)
  library(dplyr)
  library(europepmc)
  
  url1  <- "ftp://ftp.ncbi.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
  
  x <- read_delim(url1, delim = "\t", na = "-", quote = "")
  
  #replace  #Organism/Name in column 1
  names(x)[1] <- "Organism"
  
  #21.6K with genome publication
  
  table(is.na( x$`Pubmed ID`))
  
  #The replicon field isn’t great, but you can use grep to find matches
  #  SEARCH for accession in Replicon
  
  # y <- filter( x, grepl( "NC_003143",  Replicons) )
  
  y <- filter( x, grepl( acc_num,  Replicons) )
  
  
  ## Run queries using Ropensci’s package
  
  res <- vector("list", 4)
  
  ## full names Or abbreviations?
  names(res) <- c("Genome paper", "Cites genome paper", "Mentions genome accessions",  "Matches keywords")
  names(res) <- c("*", "G", "A", "K")
  
  
  #1. Get Genome publication
  
  q1 <- paste0( "EXT_ID:", y$`Pubmed ID` )
 
  
  res[[1]] <- epmc_search( q1, limit = limit)
  
  
  #  100 is default limit in epmc_search  (change if needed)
  res$citedByCount
  
  
  #2. Cites Genome publication
  
  q2 <- paste0( "cites:", y$`Pubmed ID`, "_MED")
  
  res[[2]] <- epmc_search( q2, limit = limit)
  
  
  #3. Mentions genome accession in Full text...
  
  # create a function to parse replicon field and add WGS or bioproject acc (rarely cited).
  ## add SRC:MED to return publications with Pubmed ID
  
  y$Replicons
 
  format_acc <- function(y){
    ids <- gsub(".*chromosome:([^.]+).*/([^.]+).*", "\\1 OR \\2", y$Replicons)
    if(!is.na( y$WGS))   ids <- paste0(ids, " OR ", substr(y$WGS, 1,5), "*")
    paste0( "(", ids, ") AND SRC:MED")
  }
  
  q3 <- format_acc(y)
 
  res[[3]] <- epmc_search( q3, limit = limit)
  
  
  # 4   Keywords
  
  # I would maybe drop this, since it often returns too many unrelated articles.  
  #In this example, the organism name is missing the strain PR1, so you get 15 results
  #instead of 3 results (although most organisms have strain names)
  
  q4 <-    paste( y$Organism, "genome AND SRC:MED")        
  #[1] "Algoriphagus machipongonensis genome AND SRC:MED"
  res[[4]] <- epmc_search(q4, limit = limit)
  
  # COMBINE
  y <- bind_rows(res, .id = "Evidence")
  z <- group_by_(y,  .dots=setdiff(names(y), "Evidence")) %>% summarize( Evidence = paste(Evidence, collapse=", "))
  
  ## again, abbreviations might be easier...
  table(z$Evidence)
  authors_etal(z$authorString)
  bib_format(z)
  return (z)
  
}
shinyServer(function(input, output) {
  x <- eventReactive( input$update, {
    withProgress({     
      setProgress(message = "Searching Europe PMC...")
      # CHECK query
      req(input$acc_num)
      query <- input$acc_num
      #n <- input$obs
      n=100
      if(n > 2500) n <- 2500
      x <- new.env(hash=T, parent=emptyenv())
      ids<-strsplit(query, "\\s+")[[1]]
      for(i in ids) {
        #cat(file=stderr(), i,"test\n")
        y = getData(i,n)
        assign(i,y,x)
      }
      #x = getData(query,n)
      x
    })
  })
  
  # Filter data based on selections
  output$table <-DT::renderDataTable({
    df <-x()
    i <- 0
    for(y in ls(df)) {
      #cat(file=stderr(), y,"got\n")
      if(i == 0) {
        df1 <- df[[y]]
      } else {
        df1 <- bind_rows(df1,df[[y]])
      }
      i <- i+1
    }
    DT::datatable(DT_format(df1), escape = c(1,5) )
  })
  
  output$plot <- renderDygraph({
    df <-x()
    i <- 0
    for(y in ls(df)) {
      #cat(file=stderr(), y,"got\n")
      if(i == 0) {
        df1 <- df[[y]]
        df1["acc"] <- NA
        df1$acc <- y
      } else {
        df2 <- df[[y]]
        df2["acc"] <- NA
        df2$acc <- y
        df1 <- bind_rows(df1,df2)
      }
      i <- i+1
    }
    z <- table(df1$acc,df1$pubYear)
    y <- as.data.frame.matrix(z)
    #if plotting total citations rather than noisy citations each year
    y <- t(apply(y, 1, cumsum))
    z <- data.frame( year=as.numeric( colnames(y)), t(y), check.names =FALSE)
    
    dygraph(  z ) %>%
      dyAxis("x" , label ="Year", drawGrid=FALSE)  %>%
      dyAxis("y" , label ="Total citations", drawGrid=FALSE)  %>%
      dyOptions(logscale= log )  %>%
      dyHighlight(highlightSeriesOpts = list(strokeWidth = 3))  %>%
      dyLegend(width=350) %>%
      dyCSS( "dygraph.css")
   
  })
  
})


#' Format a Javascript DataTable
#'
#' Format a data.frame to display using datatable in the DT library
#'
#' @param x Europe PMC search results
#' @param authors Number of authors to display before adding et al.
#' @param issue Include issue number with journal citation
#' @param links Add html links to PubMed ID, Journal, and Cited by counts, default TRUE
#'
#' @return a data.frame with pmid, authors, year, title, journal and cited by counts
#' @note Requires the \code{DT} package for displaying tables
#' @author Chris Stubben
#'
#' @examples
#' data(yp)
#' DT_format(yp[6:8,])
#' \dontrun{
#' x1 <- DT_format(yp)
#' library(DT)
#' datatable(x1, escape = c(1,5) , caption= "Publications with Yersinia pestis virulence in title")
#' }
#' @export
DT_format <- function(x, authors=3, issue=FALSE, links=TRUE ){
  n1 <- grep("^author", names(x) )  #authorString or authors
  authors <- authors_etal( x[[n1]], authors=authors)
  
  n1 <- grep("pubYear|year", names(x) )  #pubYear or year
  year <- x[[n1]]
  n1 <- grep("^title", names(x) )
  title <- x[[n1]]
  
  #combine journal volume pages
  journal <- journal_cite(x, issue=issue)
  
  n1 <- grep("cited", names(x) )
  citedBy <- x[[n1]]
  pmid <- x$pmid
  
  # hyperlinks
  if(links){
    citedBy <- ifelse(citedBy == 0, 0, paste('<a href="http://europepmc.org/search?query=cites%3A',
                                             x$pmid, '_MED" target="_blank">', citedBy,  '</a>', sep=""))
    pmid <- paste0('<a href="http://europepmc.org/abstract/MED/', pmid, '" target="_blank">', pmid,  '</a>')
    # some dois missing
    journal <-  ifelse(is.na(x$doi), journal, paste0('<a href="http://dx.doi.org/', x$doi, '" target="_blank">', journal,  '</a>') )
    
  }
  x <- data.frame(pmid, authors, year, title, journal, citedBy)
  x
}

#' Format a bibliography
#'
#' Format Europe PMC search results into a bibliography
#'
#' @param x epmc_search results
#' @param authors Number of authors to display before adding et al.
#' @param issue Include issue number with citation
#' @param links Add Markdown links to Journal, Cited by counts and PubMed ID, default FALSE
#' @param cited Include Cited By:<count>
#' @param pmid Include PubMed:<ID>
#' @param number Number references in bibliography
#'
#' @return a vector
#' @note Currently, the references are formatted using author, year, title and journal
#' @author Chris Stubben
#' @seealso \code{\link{strwrap}} for formatting
#' @examples
#' data(yp)
#' bib_format(yp[1:5,])
#' cat(strwrap(bib_format(yp[1:5,], number=TRUE, pmid=TRUE), exdent=3), sep="\n")
#' @export

bib_format <- function(x, authors=3, issue=TRUE, links=FALSE, cited=FALSE, pmid=FALSE, number=FALSE){
  n1 <- grep("^author", names(x) )  #authorString or authors
  authors <- authors_etal( x[[n1]], authors=authors)
  #combine journal volume issue pages
  journal <- journal_cite(x, issue=issue)
  
  
  n1 <- grep("pubYear|year", names(x) )  #pubYear or year
  year <- x[[n1]]
  n1 <- grep("^title", names(x) )
  title <- x[[n1]]
  
  n1 <- grep("cited", names(x) )
  citedBy <- x[[n1]]
  
  # Markdown links
  if(links){
    citedBy <- ifelse(citedBy == 0, 0,
                      paste('[', citedBy, '](http://europepmc.org/search?query=cites%3A', x$pmid, '_MED)', sep=""))
    x$pmid <- paste0('[', x$pmid, '](http://europepmc.org/abstract/MED/', x$pmid, ')' )
    journal <-  ifelse(is.na(x$doi), journal, paste0('[', journal, '](http://dx.doi.org/', x$doi, ')') )
  }
  ## ADD additional formats as option?
  y <-  paste0(authors, " ", year, ". ", title, " ", journal)
  
  #  year in parentheses or journal in Markdown italics
  # y <-  paste0(authors, " (", year, "). ", title, " ", journal, ".")
  #  y <-  paste0(authors, " ", year, ". ",   title, " *", journal, "*.")
  
  if(cited) y <- paste0(y, " Cited By:", citedBy)
  if(pmid)  y <- paste0(y, " PubMed:", x$pmid)
  if(number){
    n1 <- sprintf(paste0("%", nchar(length(y)), "s"), 1:length(y) )
    y <- paste(n1, y, sep=". ")
  }
  y
}

#' Format author strings
#'
#' Format delimited lists of authors and replace others with et al.
#'
#' Replaces 2 or more authors with et al.
#'
#' @param x A vector of author names
#' @param authors Number of authors to display before adding et al.
#' @param split Character used for splitting author strings, default comma
#'
#' @return a vector
#'
#' @author Chris Stubben
#'
#' @examples
#' authors_etal("Kawasaki S, Mizuguchi K,Sato M, Kono T, Shimizu H.")
#' authors_etal("Kawasaki S, Mizuguchi K,Sato M, Kono T, Shimizu H.", 2)
#' @export

authors_etal <- function(x, authors=3, split=", *"){
  y <- strsplit(x, split)
  sapply(y, function(x){
    if(length(x) > (authors + 1))  x <- c(x[1:authors], "et al.")
    paste(x, collapse=", ")
  })
}

#' Format journal citations
#'
#' Format journal title, volume, issue and pages into single string
#'
#' Formats citations with or without the issue number like MBio 5(1):13  or MBio 5:13
#'
#' @param x a data.frame with output from \code{epmc_search}
#' @param n a vector with 4 columns positions for journal, volume, issue and page,
#'         will match \code{epmc_search} output by default
#' @param issue include issue with citation, default TRUE
#'
#' @return a vector
#'
#' @author Chris Stubben
#'
#' @examples
#' data(yp)
#' journal_cite(yp[1:5,])
#' journal_cite(yp[1:5,], issue=FALSE)
#' @export

journal_cite <- function( x, n, issue=TRUE){
  if(missing(n))
  {
    # search_core output
    n <- pmatch(c("nlm_ta", "volume", "issue", "page"), names(x))
    if(any(is.na(n))){  # epmc_search from europepmc
      n <- pmatch(c("journalTitle", "journalVolume", "issue", "page"), names(x))
    }
    if( any(is.na(n))) stop("Cannot match journal, volume, issue, page columns")
  }
  
  if(issue){
    y <- apply(x[,n], 1, function(y) paste(y[1], " ", y[2], "(", y[3], "):", y[4], sep ="") )
    y <- gsub("(NA)", "", y, fixed=TRUE)  # missing issue
  }else{
    y <- apply(x[,n], 1, function(y) paste(y[1], " ", y[2], ":", y[4], sep ="") )
  }
  y <- gsub(":NA$", "", y )  # missing pages
  y <- gsub(" NA", " ", y )  # missing journal?
  y
}

#' Publication year time series
#'
#' Create a time-series object from publication years
#'
#' @param x a vector of years, or optionally epmc_search results
#' @param start year of first observation, default minimum year in series
#' @param end year of last observation, default today's year
#' @param sum return cumulative sum
#'
#' @return a time-series object
#' @seealso \code{\link{ts}}
#' @author Chris Stubben
#'
#' @examples
#' data(yp)
#' y <- year_ts(yp)
#' plot(y, xlab="Year published", ylab="Articles per year", las=1,
#'     main="Publications with Yersinia pestis virulence in title")
#' @export

year_ts <- function(x, start, end, sum=FALSE){
  if(is.data.frame(x)){
    n1 <- grep("pubYear|year", names(x) )  #pubYear or year
    x <- x[[n1]]
  }
  y <- as.numeric(x)
  if(length(y)==0) stop("No pub years found")
  if(missing(start))  start <- min(y)
  if(missing(end))  end <- format(Sys.Date(), "%Y")
  y <- factor(y, start:end)
  y <- table(y)
  if(sum) y <- cumsum(y)
  
  z <- stats::ts(as.vector(y), start= start )
  z
}