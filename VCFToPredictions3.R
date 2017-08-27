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


        
        
#Get patient information         
PatientID=as.character(commandArgs(trailingOnly = T)[[1]][1])

InputAddress=as.character(commandArgs(trailingOnly = T)[[2]][1])
InputAddress=paste0(InputAddress,"/Patient_",PatientID)

OutputAddress=as.character(commandArgs(trailingOnly = T)[[3]][1])
OutputAddress=paste0(OutputAddress,"/Patient_",PatientID)
 # PatientID='010'
 # InputAddress='/Users/dliu/BoxSync/GeneTank/IntegratedPipeline/VCF_files/Patient_010'
 # OutputAddress='/Users/dliu/BoxSync/GeneTank/IntegratedPipeline/Prediction'

#Threshold=as.numeric(commandArgs(trailingOnly = T)[[4]][1])
Threshold=c(0,-2)
#OutputAddress=paste0(InputAddress,"/Patient_",PatientID)
ReferenceAddress="./reference/"

print(paste("Patient ID:", PatientID))
print(paste("InputAddress:", InputAddress))
print(paste("OutputAddress:", OutputAddress))
print(paste("Threshold:",Threshold))




#read vcf
VCF=paste0(InputAddress,'/',list.files(InputAddress)[grepl(pattern = '.vcf',list.files(InputAddress))])
Target=readVcf(VCF,"hg19")
print('read VCF')
#geno(Target)



#Tidy up Genotype 
GT <-geno(Target)$GT
tail(GT)
class(GT)
GT=as.data.frame(GT)
GT$ref=as.vector(ref(Target))
GT$alt=as.vector(unlist(alt(Target)))
##There some SNP with only one allele , assuming the other one is zero

GT$Dosage=unlist(lapply(strsplit(as.character(GT$GENOTYPE),'[/]'),FUN = function(x)(sum(as.integer(x)))))
GetGeno=function(i){
        c(rep(GT$ref[i],2-GT$Dosage[i]),rep(GT$alt[i],GT$Dosage[i]))
}


#re-organist the genotype for searching SNPedia
GT$Genotype2=unlist(lapply(1:nrow(GT),FUN=function(x){paste(GetGeno(x),collapse = ';')}))
GT$Genotype2=paste0('(',GT$Genotype2,")")
GT$Genotype2=paste0(rownames(GT),GT$Genotype2)
tail(GT)

#Search Snpedia 

if (!require("SNPediaR",character.only = TRUE))
{
        source("http://bioconductor.org/biocLite.R") 
        biocLite("SNPediaR") 
        
        if(!require("SNPediaR",character.only = TRUE)) stop("Package not found")
}
library("SNPediaR") #load the package


#only a very small portion of SNPs are avaible in SNPedia 26th Aug 2017 
X=load('SNPavailableInSnpedia.Rdata')
AvilSNPs=get(X)
GT=GT[rownames(GT)%in%tolower(AvilSNPs),]

print(paste(sum(rownames(GT)%in%tolower(AvilSNPs)),"SNPs available in SNPedia"))

Vec=floor(seq(from=1,to=nrow(GT),length.out = 10))

res<-list()
for (i in (1:(length(Vec)-1))){
        print(paste("working on SNP",Vec[i],'to',Vec[i+1]))
        res=c(res,getPages(GT$Genotype2)[Vec[i]:Vec[i+1]])
        }

##remove those without matches 
res=res[sapply(res,length)>0]
length(res)
GetSnpTags=function(x){extractTags(x,tags = c("rsid", "allele1","allele2",'magnitude',"repute",'Gene',"summary"))}

SnpediaResult=t(sapply(res,GetSnpTags))
SnpediaResult=as.data.frame(SnpediaResult)

head(SnpediaResult)
SnpediaResult<-SnpediaResult[!is.na(SnpediaResult$summary),]

print('Data collection fron Snpedia done')


#Bad =-1, good =+1, NA =0
SnpediaResult$reputeScore=0
SnpediaResult$reputeScore[SnpediaResult$repute=='Good']=1
SnpediaResult$reputeScore[SnpediaResult$repute=='Bad']=-1
SnpediaResult$reputeScore[is.na(SnpediaResult$repute)]=0
SnpediaResult$HealthScore=SnpediaResult$reputeScore*as.numeric(as.character(SnpediaResult$magnitude))
#SnpediaResult$HealthScore[is.na(SnpediaResult$HealthScore)]=0



#Match with patient like me disease 
PatientLikeMe=read.csv('PatientsLikeMe.csv')
print('loading PatientLikeMe data')



head(PatientLikeMe)
nrow(PatientLikeMe)
FinalReport=PatientLikeMe
FinalReport$SnpediaScore=0
FinalReport$Risk='Unknown'
for (i in 1:nrow(PatientLikeMe)){
Vec=grep(PatientLikeMe$Phenotype_Name[i],SnpediaResult$summary,ignore.case = T)
HealthScore=sum(SnpediaResult[Vec,'HealthScore'])
FinalReport$SnpediaScore[i]=HealthScore
if(is.na(HealthScore)){FinalReport$Risk[i]='Unknown'}
if(!is.na(HealthScore)){
if(HealthScore>=Threshold[1]){FinalReport$Risk[i]='Low Risk'}
if(HealthScore<Threshold[1]&HealthScore>=Threshold[2]){FinalReport$Risk[i]='Moderate Risk'}
if(HealthScore<Threshold[2]){FinalReport$Risk[i]='High Risk'}
        }
}
table(FinalReport$Risk)
table(FinalReport$SnpediaScore)
head(FinalReport)


#Category of phenotypes 
FinalReport$Category=sample(c(1:5),size=nrow(FinalReport),replace = T)

#head(FinalReport)
#tail(DiseaseReport)
dir.create(OutputAddress, showWarnings = FALSE)
write.csv(FinalReport,file = paste0(OutputAddress,'/GeneticPhenotypeReport_','Patient_',PatientID,'.csv'))
print('finished output ')

