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

# General function to return postions of a given pair of AA from a AA seq
# e.g. calculate STQ numbers and positions based on AA seq
aa_num_pos <- function(seq, pattern="[ST]Q")
{
	stqs <- unlist(gregexpr(pattern, seq))
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

# calculate every pair of AA occupacy in any given proteome sequence vector
pair_mat <- function(pro_seq_vec)
{
	proAAList <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
				   "Q", "R", "S", "T", "V", "W", "Y")
	aaPairMat <- data.frame(first_aa=character(400), second_aa=character(400),
					pair_sum=integer(400), stringsAsFactors=FALSE)
	aa_seq_vec <- pro_seq_vec
	currRow <- 1
	for (i in proAAList)
	{
		for(j in proAAList)
		{
			aaPairMat$first_aa[currRow] <- i
			aaPairMat$second_aa[currRow] <- j
			currSum <- 0
			currPair <- paste(c(i,j), collapse="")
			for (k in 1:length(aa_seq_vec))
			{
				cat("Pair:", currPair, "(", currRow, ")", " Seq:", k, "\r")
				result <- unlist(gregexpr(currPair, aa_seq_vec[k]))
				resultLen <- length(result)
				currSum <- currSum + ifelse((resultLen == 1 && result == -1), 0, resultLen)
			}
			aaPairMat$pair_sum[currRow] <- currSum
			currRow <- currRow + 1
		}
	}
	cat("\n")
	return(aaPairMat)
}


# ===== temporary code ======
