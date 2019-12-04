pkgVersionCRAN = function(pkg, cran_url='http://cran.r-project.org/web/packages/')
{

  # Create URL
  cran_pkg_loc = paste0(cran_url,pkg)

  # Try to establish a connection
  suppressWarnings( conn <- try( url(cran_pkg_loc) , silent=TRUE ) )

  # If connection, try to parse values, otherwise return NULL
  if ( all( class(conn) != "try-error") ) {
    suppressWarnings( cran_pkg_page <- try( readLines(conn) , silent=TRUE ) )
    close(conn)
  } else {
    return(NULL)
  }

  # Extract version info
  version_line = cran_pkg_page[grep("Version:",cran_pkg_page)+1]
  as.package_version(gsub("<(td|\\/td)>","",version_line))
}
