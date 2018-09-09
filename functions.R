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

# calculate relative coordinates
# gap is 181 for budding yeast, 190 for Ce and 155 for mouse
# coordinate format: x1,y1;x2,y2;...xn,yn;
# for the ease of python, scale_to means 0~9
calc_coordinates <- function(seq, pattern="[ST]Q", gap=181, scale_to=5)
{
	locs <- unlist(gregexpr(pattern, seq))
	locsLen <- length(locs)
	if (locsLen == 1 && locs == -1)
		return(as.character(NA))
	else
	{
		result <- ""
		seqLen <- nchar(seq)
		locs <- c(-gap, locs, gap+seqLen)
		scaleFactor <- gap / scale_to
		for (i in 2:(locsLen+1))
		{
			currX <- locs[i] - locs[i-1]
			currY <- locs[i+1] - locs[i]
			currX <- ifelse(currX >= gap, scale_to-1, floor(currX/scaleFactor))
			currY <- ifelse(currY >= gap, scale_to-1, floor(currY/scaleFactor))
			result <- paste(c(result,currX,",",currY,";"), sep="", collapse="")
		}
		return(result)
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

# trim the last asterisk symbol from all protein sequences
trim_last_asterisk <- function(seq)
{
	seqLen <- nchar(seq)
	lastChar <- substr(seq, seqLen, seqLen)
	return (ifelse(lastChar == "*", substr(seq, 1, seqLen-1), seq))
}

# parse coordinates to (normalized) vector
parse_coordinate_to_vec <- function(coordinateString, matSize=5)
{
	splitted <- unlist(strsplit(coordinateString, ";"))
	result <- matrix(rep(0, matSize*matSize), ncol=matSize)
	for (i in splitted)
	{
		secondSplitted <- unlist(strsplit(i, ","))
		result[as.integer(secondSplitted[1])+1, as.integer(secondSplitted[2])+1] =
			result[as.integer(secondSplitted[1])+1, as.integer(secondSplitted[2])+1] + 1
	}
	result <- result / sum(result)
	dim(result) <- NULL
	return(result)
}

# ===== temporary code ======
for (i in 1:nrow(df))
{
	if (i == 1)
		temp <- parse_coordinate_to_vec(df$coordinates[i], matSize=5)
	else
		temp <- rbind(temp, parse_coordinate_to_vec(df$coordinates[i],matSize=5))
}



