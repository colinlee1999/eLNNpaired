# returns FDR, FNDR, FNR, FPR
rough_analyze <- function(m)
{	
	npositive = 0
	npositive_recognize = 0
	npositive_not_recognize = 0
	false_positive = 0

	nnegative = 0
	nnegative_recognize = 0
	nnegative_not_recognize = 0
	false_negative = 0

	G = nrow(m)

	for (i in 1:G)
	{
		if (m[i,1]==1)
		{
			npositive = npositive +1
		}
		if (m[i,2]==1)
		{
			npositive_recognize = npositive_recognize + 1
		}
		if (m[i,1]==1 && m[i,2]!=1)
		{
			npositive_not_recognize = npositive_not_recognize + 1
		}
		if (m[i,1]==2)
		{
			nnegative = nnegative + 1
		}
		if (m[i,2]==2)
		{
			nnegative_recognize = nnegative_recognize + 1
		}
		if (m[i,1] == 2 && m[i,2]!=2)
		{
			nnegative_not_recognize = nnegative_not_recognize + 1
		}
		if (m[i,1]!=1 && m[i,2] == 1)
		{
			false_positive = false_positive + 1
		}
		if (m[i,1]!=2 && m[i,2] == 2)
		{
			false_negative = false_negative + 1
		}
	}

	return (c(G,npositive,npositive_recognize,npositive_not_recognize,false_positive,nnegative,nnegative_recognize,nnegative_not_recognize,false_negative))
}