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

# 
# #Generate CSF files based on 23andMe data 
# 
PatientID=as.character(commandArgs(trailingOnly = T)[[1]][1])

InputAddress=as.character(commandArgs(trailingOnly = T)[[2]][1])
InputAddress=paste0(InputAddress,"Patient_",PatientID)

OutputAddress=as.character(commandArgs(trailingOnly = T)[[3]][1])

Threshold=as.numeric(commandArgs(trailingOnly = T)[[4]][1])
#OutputAddress=paste0(InputAddress,"/Patient_",PatientID)
ReferenceAddress="./reference/"

print(paste("Patient ID:", PatientID))
print(paste("InputAddress:", InputAddress))
print(paste("OutputAddress:", OutputAddress))
print(paste("Threshold:",Threshold))
# 


#read vcf
VCF=paste0(InputAddress,'/',list.files(InputAddress)[grepl(pattern = '.vcf',list.files(InputAddress))])
Target=readVcf(VCF,"hg19")
print('read VCF')
#geno(Target)



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


GetType=function(x,Threhold){
        AllTypes=x$Risk
        x$PatientLikeMeName=as.character(x$PatientLikeMeName)
       if(all(AllTypes=='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Unknown'))}
        if(!all(AllTypes=='Unknown')){
       if(sum(AllTypes=="Positive")>=Threhold*sum(AllTypes!='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Positive'))} 
       if(sum(AllTypes=="Positive")<Threhold*sum(AllTypes!='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Negative'))}
                }
        
}

DiseaseReport=lapply(DiseaseReport,FUN = function(x){GetType(x,Threshold)})
DiseaseReport=Reduce(DiseaseReport,f = rbind)
rownames(DiseaseReport)=NULL
DiseaseReport=as.data.frame(DiseaseReport)
#head(DiseaseReport)
#sum(DiseaseReport$Risk=='Unknown')
table(DiseaseReport$Risk)

#tail(DiseaseReport)
write.csv(DiseaseReport,file = paste0(OutputAddress,'/GeneticPhenotypeReport_','Patient_',PatientID,'.csv'))
print('finished output ')

