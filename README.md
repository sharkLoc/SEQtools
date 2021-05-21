# SEQtools
filt and stats PE fastq files

##usage
Usage: SEQtools [options]

	-i	:<char>   fq1 file
	-I	:<char>   fq2 file
	-l	:<int>    low quality threshold (default: [10])
	-m	:<int>    filter reads with low average quality less than [30]
	-q	:<float>  low quality rate (default: [0.50])
	-n	:<float>  N rate threshold (default: [0.05])
	-o	:<char>   clean fq1 file name
	-O	:<char>   clean fq2 file name
	-h	:<char>   show this help
