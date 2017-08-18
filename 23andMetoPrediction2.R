#read VCF files, matching GWAS data
#source("http://bioconductor.org/biocLite.R") 
#biocLite("VariantAnnotation") #install the package

#check package 

        if (!require("VariantAnnotation",character.only = TRUE))
        {
                source("http://bioconductor.org/biocLite.R") 
                biocLite("VariantAnnotation") 
                
                if(!require("VariantAnnotation",character.only = TRUE)) stop("Package not found")
        }
library("VariantAnnotation") #load the package


#Generate CSF files based on 23andMe data 

PatientID=as.character(commandArgs(trailingOnly = T)[[1]][1])
InputAddress=as.character(commandArgs(trailingOnly = T)[[2]][1])
OutputAddress=as.character(commandArgs(trailingOnly = T)[[3]][1])
OutputAddress=paste0(OutputAddress,"/Patient_",PatientID)
ReferenceAddress="./reference/"

system(paste("python ./23andme/23andmeVersionBuild.py", InputAddress))

print('version detected')

system(paste("python ./23andme/23andmeVersionBuild_convert_liftover2.py", InputAddress,OutputAddress,ReferenceAddress))

print('lifted over and coverted to CVF')


#read vcf
VCF=paste0(OutputAddress,'/',list.files(OutputAddress)[grepl(pattern = '.vcf',list.files(OutputAddress))])
Target=readVcf(VCF,"hg19")
print('read VCF')
#geno(Target)

#geno(header(Target))

#header information 
#header(Target)

#samples(header(Target))

#geno(header(Target))

#Genomic position 
#head(rowRanges(Target),3)
#ref(Target)[1:5]
#alt(Target)[1:5]
#qual(Target)[1:5]


#Genotype 

#geno(Target)


#sapply(geno(Target),class)
#geno(header(Target))["DS",]

#DS <-geno(Target)$DS
#dim(DS)
GT <-geno(Target)$GT

#GT[1:3,]
#length(GT)
GTvariants=unlist(lapply(strsplit(GT,'[/]'),FUN = function(x){sum(as.integer(x))}))
#length(GTvariants)
#head(GTvariants)
names(GTvariants)=rownames(GT)

#Get the GWAS informaiton 
GWAScatalog=read.csv('./GWASCatalogLinkedToPNMe.csv')
#head(GWAScatalog)
#names(GWAScatalog)
GWASsubInfor=GWAScatalog[,c('PatientLikeMeName',"SNPS")]
#head(GWASsubInfor)

#matching , 3 situation, mathcing >0, matching =0, non-matching 
DiseaseReport=GWASsubInfor
DiseaseReport$Risk='Unknown'
DiseaseReport=DiseaseReport[!is.na(DiseaseReport$PatientLikeMeName),]
#head(DiseaseReport)
#nrow(DiseaseReport)
VecScore=GTvariants[as.character(DiseaseReport$SNPS)]
DiseaseReport$Risk[which(VecScore>0)]='Positive'
DiseaseReport$Risk[which(VecScore==0)]='Negative'

#sum(DiseaseReport$Risk=='Positive')
#sum(DiseaseReport$Risk=='Negative')
#sum(DiseaseReport$Risk=='Unknown')

#Aggregate the results
DiseaseReport=split(DiseaseReport,f = DiseaseReport$PatientLikeMeName)
#length(DiseaseReport)


GetType=function(x){
        AllTypes=unique(x$Risk)
        x$PatientLikeMeName=as.character(x$PatientLikeMeName)
       if("Positive"%in%AllTypes){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Positive'))} 
       if((!"Positive"%in%AllTypes)&("Negative"%in%AllTypes)){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Negative'))}
        if((AllTypes=='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Unknown'))}
}

DiseaseReport=lapply(DiseaseReport,FUN = GetType)
DiseaseReport=Reduce(DiseaseReport,f = rbind)
rownames(DiseaseReport)=NULL
DiseaseReport=as.data.frame(DiseaseReport)
#head(DiseaseReport)
#sum(DiseaseReport$Risk=='Unknown')

#tail(DiseaseReport)
write.csv(DiseaseReport,file = paste0(OutputAddress,'/GeneticPhenotypeReport_','Patient_',PatientID,'.csv'))
print('finished output ')

