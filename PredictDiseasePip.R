
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
#head(PatientData)

Target=PatientData[PatientData$Phenotype%in%AllPhenotype,]
Target=Target[Target$Value==1,]
#head(Target)
Predicted=unique(c(as.character(DNet$Disease2[DNet$Disease1%in%Target$Phenotype]),
                   as.character(DNet$Disease1[DNet$Disease2%in%Target$Phenotype])))
Predicted=Predicted[!Predicted%in%Target$Phenotype]

Report=Target
Report=Report[,-1]
rownames(Report)=NULL
Report$type='Reported'

#head(Report)
Report=rbind(Report,data.frame(Phenotype=Predicted,Value=1,type='Predicted'))
dir.create(OutputAddress, showWarnings = FALSE)


write.csv(Report,file = paste0(OutputAddress,'/PhenoPre_report_patient_',PatientID,'.csv'))

print('done')
      