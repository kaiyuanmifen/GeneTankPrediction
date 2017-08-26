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
InputAddress=paste0(InputAddress,"/Patient_",PatientID)

OutputAddress=as.character(commandArgs(trailingOnly = T)[[3]][1])
OutputAddress=paste0(OutputAddress,"/Patient_",PatientID)
# PatientID='010'
# InputAddress='/Users/dliu/BoxSync/GeneTank/IntegratedPipeline/VCF_files/Patient_010'
# OutputAddress='/Users/dliu/BoxSync/GeneTank/IntegratedPipeline/Prediction'

#Threshold=as.numeric(commandArgs(trailingOnly = T)[[4]][1])
Threshold=c(1,0.5)
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
head(GT)
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

GetScore=function(x){
        AllTypes=x$Risk
        x$PatientLikeMeName=as.character(x$PatientLikeMeName)
        if(all(AllTypes=='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Score='Unknown'))}
        if(!all(AllTypes=='Unknown')){
                return(c(PhenotypeName=x$PatientLikeMeName[1],Score=sum(AllTypes=="Positive")/sum(AllTypes!='Unknown')))
               
        }
        
}


# GetType=function(x,Threhold){
#         AllTypes=x$Risk
#         x$PatientLikeMeName=as.character(x$PatientLikeMeName)
#        if(all(AllTypes=='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Unknown'))}
#         if(!all(AllTypes=='Unknown')){
#        if(sum(AllTypes=="Positive")>=Threhold[1]*sum(AllTypes!='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='High Risk'))} 
#        if((sum(AllTypes=="Positive")<Threhold[1]*sum(AllTypes!='Unknown'))&(sum(AllTypes=="Positive")>=Threhold[2]*sum(AllTypes!='Unknown'))){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Moderate Risk'))}
#        if(sum(AllTypes=="Positive")<Threhold[2]*sum(AllTypes!='Unknown')){return(c(PhenotypeName=x$PatientLikeMeName[1],Risk='Low Risk'))}
#                 }
#         
# }

DiseaseReport=lapply(DiseaseReport,FUN = function(x){GetScore(x)})
head(DiseaseReport)

#DiseaseReport=lapply(DiseaseReport,FUN = function(x){GetType(x,Threshold)})
DiseaseReport=Reduce(DiseaseReport,f = rbind)
rownames(DiseaseReport)=NULL
DiseaseReport=as.data.frame(DiseaseReport)
DiseaseReport$Score=as.numeric(as.character(DiseaseReport$Score))
#head(DiseaseReport)
#sum(DiseaseReport$Risk=='Unknown')
#DiseaseReport$Risk='Unknown'
#DiseaseReport$Risk[DiseaseReport$Score>=Threshold[1]]='High Risk'

#DiseaseReport$Risk[DiseaseReport$Score<Threshold[1]&DiseaseReport$Score>=Threshold[2]]='Moderate Risk'
#DiseaseReport$Risk[DiseaseReport$Score<Threshold[2]]='Low Risk'
#head(DiseaseReport)
#DiseaseReport[is.na(DiseaseReport$Score),]
#table(DiseaseReport$Risk)
#tail(DiseaseReport)

#ClinVar
##Read ClinVar disease data 
X=load('ClinVAR.Rdata')
ClinVar=get(X)
ClinVar=ClinVar[!is.na(ClinVar$ID),]
#head(ClinVar)
ClinVarReport=ClinVar
head(ClinVarReport$ID)
head(GTvariants)
VecScore=rep(NA,nrow(ClinVarReport))
Vec=match(names(GTvariants),ClinVarReport$ID)
VecScore[Vec[!is.na(Vec)]]=GTvariants[which(!is.na(Vec))]
head(VecScore)
ClinVarReport$Risk='Unknown'
ClinVarReport$Risk[which(VecScore>0)]='Positive'
ClinVarReport$Risk[which(VecScore==0)]='Negative'

head(ClinVarReport)
ClinVarReport=split(ClinVarReport,f = as.factor(as.character(ClinVarReport$PatientLikeMeName)))
length(ClinVarReport)
ClinVarReport=lapply(ClinVarReport,FUN = GetScore)
ClinVarReport=Reduce(ClinVarReport,f = rbind)
ClinVarReport=as.data.frame(ClinVarReport)
ClinVarReport$Score=as.numeric(as.character(ClinVarReport$Score))
tail(ClinVarReport)

rownames(ClinVarReport)=NULL


#Combine twp reports 
PatientLikeMe=read.csv('PatientsLikeMe.csv')
head(PatientLikeMe)
FinalReport=data.frame(Phenotype=PatientLikeMe$Phenotype_Name,GWAS_Score=NA,ClinVarScore=NA)
FinalReport$Phenotype=tolower(FinalReport$Phenotype)
Vec=match(as.character(FinalReport$Phenotype),as.character(DiseaseReport$PhenotypeName))
FinalReport$GWAS_Score[Vec[!is.na(Vec)]]=DiseaseReport$Score[!is.na(Vec)]

Vec=match(as.character(FinalReport$Phenotype),ClinVarReport$PhenotypeName)
FinalReport$ClinVarScore[Vec[!is.na(Vec)]]=ClinVarReport$Score[!is.na(Vec)]

FinalReport$MaxScore=apply(as.matrix(FinalReport[,c(2,3)]),MARGIN = 1,function(x){max(x,na.rm = T)})
FinalReport$MaxScore[is.infinite(FinalReport$MaxScore)]=NA
table(FinalReport$MaxScore)
FinalReport$Risk="Unknown"
FinalReport$Risk[FinalReport$MaxScore>=Threshold[1]]='High Risk'
FinalReport$Risk[FinalReport$MaxScore<Threshold[1]&FinalReport$MaxScore>=Threshold[2]]='Moderate Risk'
FinalReport$Risk[FinalReport$MaxScore<Threshold[2]]='Low Risk'
table(FinalReport$Risk)
#head(FinalReport)
#tail(DiseaseReport)
dir.create(OutputAddress, showWarnings = FALSE)
write.csv(FinalReport,file = paste0(OutputAddress,'/GeneticPhenotypeReport_','Patient_',PatientID,'.csv'))
print('finished output ')

