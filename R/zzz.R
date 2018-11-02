#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_attr
#' @importFrom utils unzip download.file
#' @importFrom tools file_path_sans_ext
.onLoad <- function(libname, pkgname) {

  td <- file.path(tempdir(),'torsten')
  dir.create(td,showWarnings = FALSE)
  TH <- find.package('torstenHeaders')

  getURLS <- xml2::read_html('https://github.com/metrumresearchgroup/TorstenHeaders/tree/master/inst')%>%
    rvest::html_nodes(xpath = '//*[contains(concat( " ", @class, " " ), concat( " ", "css-truncate-target", " " ))]//span//a')%>%
    rvest::html_attr(name = 'href')

  getURLS <- sprintf('https://github.com%s.zip',gsub('tree','archive',getURLS))

  names(getURLS) <- sapply(strsplit(getURLS,'/'),'[',5)

  zipnames <- sapply(names(getURLS),function(x) paste(x,tools::file_path_sans_ext(basename(getURLS[x])),sep='-'))

  if(length(list.files(file.path(TH,'stan')))==0){

    utils::download.file(url = getURLS['stan'],destfile = file.path(td,'/stan.zip'),quiet = TRUE)
    utils::unzip(file.path(td,'stan.zip'),exdir = td)

    file.copy(file.path(td,zipnames['stan'],'src/stan'),
              TH,
              overwrite = TRUE,
              recursive = TRUE)
  }

  if(length(list.files(file.path(TH,'math')))==0){

    download.file(url = getURLS['math'],destfile = file.path(td,'/math.zip'),quiet = TRUE)
    utils::unzip(file.path(td,'math.zip'),exdir = td)

    suppressWarnings({
      dir.create(file.path(TH,'math'))
      file.copy(file.path(td,zipnames['math'],'stan'),
                file.path(TH,'math'),
                overwrite = TRUE,
                recursive = TRUE)
    })

  }

  unlink(file.path(td,zipnames['math']),recursive = TRUE,force=TRUE)
  unlink(file.path(td,zipnames['stan']),recursive = TRUE,force=TRUE)
}

#' @importFrom xml2 read_html
#' @importFrom rvest html_nodes html_attr
#' @importFrom utils unzip download.file
#' @importFrom tools file_path_sans_ext
.onAttach <- function(libname, pkgname) {

  td <- file.path(tempdir(),'torsten')
  dir.create(td,showWarnings = FALSE)
  TH <- find.package('torstenHeaders')

  getURLS <- xml2::read_html('https://github.com/metrumresearchgroup/TorstenHeaders/tree/master/inst')%>%
    rvest::html_nodes(xpath = '//*[contains(concat( " ", @class, " " ), concat( " ", "css-truncate-target", " " ))]//span//a')%>%
    rvest::html_attr(name = 'href')

  getURLS <- sprintf('https://github.com%s.zip',gsub('tree','archive',getURLS))

  names(getURLS) <- sapply(strsplit(getURLS,'/'),'[',5)

  zipnames <- sapply(names(getURLS),function(x) paste(x,tools::file_path_sans_ext(basename(getURLS[x])),sep='-'))

  if(length(list.files(file.path(TH,'stan')))==0){

    utils::download.file(url = getURLS['stan'],destfile = file.path(td,'/stan.zip'),quiet = TRUE)
    utils::unzip(file.path(td,'stan.zip'),exdir = td)

    file.copy(file.path(td,zipnames['stan'],'src/stan'),
              TH,
              overwrite = TRUE,
              recursive = TRUE)
  }

  if(length(list.files(file.path(TH,'math')))==0){

    download.file(url = getURLS['math'],destfile = file.path(td,'/math.zip'),quiet = TRUE)
    utils::unzip(file.path(td,'math.zip'),exdir = td)

    suppressWarnings({
      file.copy(file.path(td,zipnames['math'],'stan'),
                file.path(TH,'math'),
                overwrite = TRUE,
                recursive = TRUE)
    })


  }

  unlink(file.path(td,zipnames['math']),recursive = TRUE,force=TRUE)
  unlink(file.path(td,zipnames['stan']),recursive = TRUE,force=TRUE)
}
