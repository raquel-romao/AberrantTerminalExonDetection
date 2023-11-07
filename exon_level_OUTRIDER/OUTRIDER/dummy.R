#'---
#' title: Final html
#' author: Ines Scheller, raquelromao
#' wb:
#'  input:
#'   - summaryHTML: '`sm expand("{htmlOutputPath}/{tissue}_Summary.html", htmlOutputPath=config["htmlOutputPath"], tissue=config["tissues"])`'
#' output:
#'  html_document
#'---
#+ input


#+ echo=FALSE, results="asis"
summaryHTMLFiles <- snakemake@input$summaryHTML
summaryTissues <- gsub("_Summary.html$", "", basename(summaryHTMLFiles))

cat("<h1>Summary for each tissue</h1><p>")
devNull <- sapply(summaryTissues, function(name){
  cat(paste0(
    "<a href='Summary/", name, "_Summary.html'>",
    name, " summary</a>",
    "</br>"))
})
cat("</p>")
