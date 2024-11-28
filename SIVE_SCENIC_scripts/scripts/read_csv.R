remove_space = function(title) {
	return(unlist(strsplit(title, "\\ "))[1])
}

get_intersection = function(csv_1, csv_2) {
	csv_1 = read.csv(csv_1)
	csv_2 = read.csv(csv_2)
	s_1 = unlist(lapply(as.list(csv_1$X), remove_space))
	s_2 = unlist(lapply(as.list(csv_2$X), remove_space))
	return(intersect(s_1, s_2))
}


get_rows = function(csv_1, csv_2, filename) {
	get_int = get_intersection(csv_1, csv_2)
	csv_1 = read.csv(csv_1)
	csv_2 = read.csv(csv_2)
	s_1 = unlist(lapply(as.list(csv_1$X), remove_space))
        s_2 = unlist(lapply(as.list(csv_2$X), remove_space))
	rownames(csv_1) = s_1
	rownames(csv_2) = s_2
	csv_1 = csv_1[get_int, ]
	csv_2 = csv_2[get_int, ]
	

	len = length(colnames(csv_1))
	len2  = length(colnames(csv_2))
	stopifnot(len == len2)
	csv_1 = as.matrix(csv_1[, c(2:len)])
	csv_2 = as.matrix(csv_2[, c(2:len)])	
	write.csv(abs(csv_1 - csv_2), filename)
	return(abs(csv_1 - csv_2))
}
