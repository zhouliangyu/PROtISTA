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

# calculate STQ numbers and positions based on AA seq
stq_num_pos <- function(seq)
{
	stqs <- unlist(gregexpr("[ST]Q", seq))
	if (length(stqs) == 1 && stqs == -1)
		return(c("0", "", ""))
	else
	{
		num <- as.character(length(stqs))
		pos <- paste(stqs, collapse=",")
		aas <- character(length(stqs))
		aasCounter <- 1
		for (i in stqs)
		{
			aas[aasCounter] <- substr(seq, i, i)
			aasCounter <- aasCounter+1
		}
		aas <- paste(aas, collapse=",")
		return(c(num, pos, aas))
	}
	
}

# temporary code
df$stq_num <- 0
df$stq_pos <- ""
df$stq_aa <- ""
for (i in 1:nrow(df))
{
	cat(i, "\r")
	vec <- stq_num_pos(df$pro_seq[i])
	df$stq_num[i] <- as.numeric(vec[1])
	df$stq_pos[i] <- vec[2]
	df$stq_aa[i] <- vec[3]
}
cat("\n")

