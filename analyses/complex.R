##Rscript /science/willerslev/users-shared/science-snm-willerslev-fvr124/lf_scripts/DepthFromComp_victor_hugh_tsk2.R /willerslev/scratch/lf/misc/preseq_results/020368P-HM0130-020368E1UL1D-U-8IZKN-TGAGCC-TAGCGA-191011-191011.preseq.txt 0.34 40 26842761  4812384


##these are a constants
genomesize<-3095677412
duplvlcutoff <- 0.75

Args<-commandArgs(T)
filename<-Args[1]
DoC<-as.numeric(Args[2])
MeanLength<-as.numeric(Args[3])
totalseqreads<-as.numeric(Args[4])
seqreads<-as.numeric(Args[5])

#print(file.size(gsub("preseq.txt","log",filename)))
#print(file.size(filename))

a<-read.table(filename, h=T)
if(any(is.na(a[,2]))){
    cat(totalseqreads, seqreads, "NaN", "NaN", "\n")
        }else if(!all(diff(a[,2])>=0)){
	    cat(totalseqreads, seqreads, "NaN", "NaN", "\n")
	    }else{

    ## y is the number of reads we require given the 'length of genome' and length of each read
        y<-DoC*genomesize/MeanLength
      ##duplication level
        duplevel <- 1-a[,2]/a[,1]

    ##below is the index that tells us the closest number of reads that matches our coveragte
        index <- which.min(abs(a[,2]-y))
    if (is.na(duplevel[index])){
      cat(totalseqreads, seqreads, "NaN", "NaN", "\n")

    }else{
    if(duplevel[index]>duplvlcutoff){
            cat(totalseqreads, seqreads, -log(0),duplvlcutoff,"\n")
	        }else{
		        ##if matching duplication lvl of reqnr of distinct reads > cutoff, then choose reads given by cutoff

        ##x is now the number of sequenced reads which will be closes to our y
	        x<-a[index,1]

        ##not used anymore
	        ##ReqReads<-totalseqreads*x/seqreads
		        RemReads<-(totalseqreads*x/seqreads)-totalseqreads
			        cat(totalseqreads,seqreads,RemReads,duplevel[index],"\n")
				    }
				    }
          }
