#version 3: this version increase the threshold to mark disease as "Positive"


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



