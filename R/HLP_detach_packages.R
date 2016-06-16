
##### Detach packages
# It is possible to have multiple versions of a package loaded at once (for example, if you have a development 
# version and a stable version in different libraries). To detach guarantee that all copies are 
# detached, use this function.
detach_package <- function(pkg, character.only = FALSE)
{
  for(p in pkg) {
    if(!character.only)
    {
      p <- deparse(substitute(p))
    }
    search_item <- paste("package", p, sep = ":")
    while(search_item %in% search())
    {
      detach(search_item, unload = TRUE, character.only = TRUE)
    }
  } # end for-loop
} # end function definition

