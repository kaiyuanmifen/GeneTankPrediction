
if (!require('Matrix',character.only = TRUE))
{
        install.packages('Matrix')
        if(!require('Matrix',character.only = TRUE)) stop("Package not found")
}
library('Matrix') #load the package

PatientID=as.character(commandArgs(trailingOnly = T)[[1]][1])
Address=as.character(commandArgs(trailingOnly = T)[[2]][1])

OutputAddress=as.character(commandArgs(trailingOnly = T)[[3]][1])

OutputAddress=paste0(OutputAddress,"Patient_",PatientID)

print(paste('working on address',Address))
print(paste('working on patient',PatientID))
#give prediction 
X=load('./DiseaseNetZhao2014.Rdata')
DNet=get(X)
AllPhenotype=unique(c(as.character(DNet$Disease1),as.character(DNet$Disease2)))

#Address='./Paitent_7.csv'
#PatientNumber=7""
PatientData=read.csv(Address)
tail(PatientData)
#head(PatientData)

Target=PatientData[PatientData$Phenotype%in%AllPhenotype,]
Target=Target[Target$Value==1,]

PatientVector=rep(0,length(AllPhenotype))
names(PatientVector)=AllPhenotype
PatientVector[names(PatientVector)%in%Target$Phenotype]=1
#head(Target)
#similary shated 

DNet$i=match(DNet$Disease1,AllPhenotype)
DNet$j=match(DNet$Disease2,AllPhenotype)
DNetmat=sparseMatrix(i = DNet$i,j =DNet$j,dims = c(length(AllPhenotype),length(AllPhenotype)),
                     dimnames =list(AllPhenotype,AllPhenotype) ,x = DNet$symptom.similarity.score)


#Similarity based prediciton

Report=Target
Report=Report[,-1]
rownames(Report)=NULL
Report$type='Reported'
Report$Value=as.numeric(Report$Value)
Report$Value='NA'




dim(DNetmat)
length(PatientVector)
Predicted=as.vector(DNetmat%*%PatientVector)
names(Predicted)=AllPhenotype
Predicted=Predicted[Predicted>0]
Predicted=Predicted[!names(Predicted)%in%Report$Phenotype]
Predicted=(Predicted-min(Predicted))/(max(Predicted)-min(Predicted))#normalized 
Predicted=signif(Predicted,digits = 3)
#length(Predicted)
#hist(Predicted)
# Predicted=unique(c(as.character(DNet$Disease2[DNet$Disease1%in%Target$Phenotype]),
#                    as.character(DNet$Disease1[DNet$Disease2%in%Target$Phenotype])))
# Predicted=Predicted[!Predicted%in%Target$Phenotype]


#head(Report)
Report=rbind(Report,data.frame(Phenotype=names(Predicted),Value=as.vector(Predicted),type='Predicted'))

dir.create(OutputAddress, showWarnings = FALSE)
#Report$Value='Postiive'

write.csv(Report,file = paste0(OutputAddress,'/PhenoPre_report_patient_',PatientID,'.csv'))

print('done')
      