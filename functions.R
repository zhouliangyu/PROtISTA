# trim the last isoform character
trim_iso <- function(id)
{
	splitted <- unlist(strsplit(id, "\\."))
	mixed <- splitted[2]
	result <- ""
	for (i in 1:nchar(mixed))
	{
		currChar <-substr(mixed, i, i)
		if (grepl("[0-9]", currChar))
			result <- paste(result, currChar, sep="", collapse="")
	}
	result <- paste(splitted[1], ".", result, sep="", collapse="")
	return(result)
}

