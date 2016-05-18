################################################################################################################################
# Automatize Regression Script provided by AZ&MK                                                                               #                         
# - read in data from excel workbooks                                                                                          #    
#	- format data for regression (wide to long)                                                                                  #
#	- logistic regression (drm function from drc package)                                                                        #
#  	- profile likelihood CI                                                                                                    #
#	- plots                                                                                                                      #
#       - additional results output: EC10,20,90		                                                                             #
#	- save output in excel workbooks	                                                                                           #
#                                                                                                                              #
# UPDATE 15.12.2015 >> analysis summarized per chemical (labs as replicates)                                                   #
# UPDATE 08.12.2015 >> profile likelihood EC50 (the concentration corresponding to the min of goodnessfit)                     #
# UPDATE 08.12.2015 >> add to plots: profile likelihood EC50 (the concentration corresponding to the min of goodnessfit)       #
# UPDATE 13.11.2015 >> allow missing values in the "conc" column                                                               #
#                                                                                                                              #
#                                                                                                                              #
# Diana Coman Schmid                                                                                                           #
# diana.comanschmid@eawag.ch                                                                                                   #
################################################################################################################################

# clear the workspace
rm(list = ls())

# load libraries
library("XLConnect")
library("reshape2")
library("drc")
library("lattice")
library("sfsmisc")


# define the working directory
# when copy file path from Windows Explorer please replace \ with //

dat.path <- "D://Melanie//input//"# your directory with the input files only (xlsx format)
res.path <- "D://Melanie//output//"# your directory where the output results will be stored


dat.path <- "/media/dianacs/data1/Melanie/input/"# your directory with the input files only (xlsx format)
res.path <- "/media/dianacs/data1/Melanie/output/"# your directory where the output results will be stored


# dat.path <- "/media/dianacs/data1/Melanie/input/"
# res.path <- "/media/dianacs/data1/Melanie/output/"
# collect the file names available in the workind directory
list.files(dat.path)
file.list = list.files(dat.path, pattern = "*.xlsx")
names(file.list) <- file.list


# read in the excel workbooks containing sheets for each chemical/lab
dat.in <- list()
for (f in names(file.list)){
  wb <- loadWorkbook(file.path(dat.path,f))
  ws = readWorksheet(wb, sheet = getSheets(wb))
  dat.in[[f]] <- ws
}


# remove spaces from excel workbook and worksheet names  
names(dat.in) <- gsub(" ","_",names(dat.in))
for (n in names(dat.in)){
  names(dat.in[[n]]) <- gsub(" ","_",names(dat.in[[n]]))
}

names(dat.in)
# format input (wide to long format)
# add replication information
# fit model
# automatically save results to excel workbooks-worksheets with names corresponding to the input files
# reports: 
# tables: coefficients (EC50) and model estimate,SD, CI,slope; EC 10,20,90
# the profile likelihood CI estimates are appended
# figures: model fit, residuals vs. fitted and normal QQ-plot
# the profile likelihood plots are saved separately to PDFs


dat.out <- list()
for (b in names(dat.in)){
  dat.out.t <- list()
  res.wb <- loadWorkbook(file.path(res.path,paste("Res_",b,sep="")),create=TRUE)
  for (s in names(dat.in[[b]])){
    cat(paste("Processing Workbook >>", b,"  Worksheet >>",s, "\n", sep=''))
    df <- dat.in[[b]][[s]][,-2]# omit the "stdev.conc" column
    ### Update 15.12.2015
    ### ! SUMMARY Analysis (per lab)>> uses different input file
    df.long <- na.omit(melt(df,id="conc"))# switch from wide to long format and omit "NA" entries
    rep <- unique((substring(as.character(df.long$variable), 1, 4)))# get the number of bioReps per Lab are available
    df.long$bioRep <- as.factor(paste("bioRep",substring(as.character(df.long$variable), 1, 4),sep=""))# add column with biological replication info
    colnames(df.long) <-c("conc","var","pViab","bioRep")
    ### 
    
    # fit the model
    tryCatch({
      mod1 <- drm(pViab~conc, bioRep, data=df.long, 
                  weight = c(1/df.long$pViab),
                  fct=LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "EC50")), 
                  lowerl=c(NA,(0-0.01),(100-0.01),NA), upperl=c(NA,(0+0.01),(100+0.01),NA),
                  na.action=na.exclude)
    },error=function(e){cat("Error :", conditionMessage(e),"\n")})
    dat.out.t[[s]] <- summary(mod1)
    
    # save coefficient estimates to worksheets-workbooks
    cat(paste("Saving model estimates to Workbook >>", b,"  Worksheet >>",s, "\n", sep=''))
    createSheet(res.wb, name = paste(s,"_EC50",sep=""))
    createSheet(res.wb, name = paste(s,"_EC10_20_90",sep=""))
    EC50.e <- as.data.frame(mod1$coefficients)
    EC50.e.rn <- rownames(EC50.e)
    EC50.e <- cbind(EC50.e.rn,EC50.e)
    colnames(EC50.e) <- c("Coeff.","Model Estim.")
    ECs.e <- ED(mod1,c(10,20,90),interval="delta")
    ECs.e.rn <- rownames(ECs.e)
    ECs.e <- cbind(ECs.e.rn,ECs.e)
    colnames(ECs.e) <- c("Coeff.","Model Estim.","Std.Error","Lower","Upper")
    
    writeWorksheet(res.wb,EC50.e,paste(s,"_EC50",sep=""))# save coefficients of each model (EC50)
    writeWorksheet(res.wb,ECs.e,paste(s,"_EC10_20_90",sep=""))# save coefficients of each model (EC10,20,90)
    
    
    setSheetColor(res.wb,paste(s,"_EC50",sep=""),XLC$COLOR.LIGHT_TURQUOISE)
    setSheetColor(res.wb,paste(s,"_EC10_20_90",sep=""),XLC$COLOR.GREY_25_PERCENT)
    
    
    # save figures to worksheets-workbooks
    cat(paste("Saving figures to Workbook >>", b,"  Worksheet >>",s, "\n", sep=''))
    createName(res.wb, name = paste(s,"_EC50",sep=""), formula = paste(s,"_EC50","!$B$40",sep=""))
    png(file = paste(s,".png",sep=""),width=1000,height=700)
    devAskNewPage(ask = FALSE)
    
    par(mfrow=c(2,2),mar=c(5.1, 4.1, 4.1, 2.1))
    
    plot(mod1,type="all",col=rainbow(length(rep)),
         ylim=c(-10,120), xlim=c(0.1,100),
         main=s,cex.main=0.80,
         ylab="% cell viability",xlab="concentration [mg/L]",
         legendPos = c(120, 120),cex.legend=0.75)
    
    plot(resid(mod1)~predict(mod1),pch=18,col="gray",main="Residuals vs. Fitted",cex.main=0.80,
         ylab="residuals",xlab="fitted",)# plot residuals 
    abline(h=0,col="orange")
    
    qqnorm(resid(mod1),pch=18,col="gray",main="Normal QQ-plot",cex.main=0.80)# QQ plot 
    qqline(resid(mod1),col="orange")
    dev.off()
    
    
    addImage(res.wb,paste(s,".png",sep=""),name=paste(s,"_EC50",sep=""),originalSize=TRUE)
    
    # calculate 95 % CI with profile likelihood 
    
    pdf(file.path(res.path,paste(s,"_ProfileL_CI.pdf",sep="")),width=11, height=8.5,pointsize=12, paper='special')
    
    par(mfrow=c(3,2))
    rep.n <- unique(df.long$bioRep)
    names(rep.n) <- unique(df.long$bioRep)

    
    for (br in names(rep.n)){
      cat(paste("Calculating 95% Profile Likelihood CI for Workbook >>", b,"  Worksheet >>",s,"Replicate",br, "\n", sep=''))
      conc <- seq(0.1, 100, 0.1) 
      goodnessfit <- rep(NA, length(conc)) #you create a vector with the same size (NA) as steps in "conc" (otherwise potentially memory problem)
      for(i in 1:length(conc)) {
        fixEC50 <- conc[i]
        tryCatch({
          mod1.pl <- drm(pViab~conc, data=df.long[which(df.long$bioRep == br),], 
                         weight = c(1/df.long[which(df.long$bioRep == br),"pViab"]),
                         fct=LL.4(names=c("Slope", "Lower Limit", "Upper Limit", "EC50")), 
                         lowerl=c(NA,(0-0.01),(100-0.01),(fixEC50-0.01)), 
                         upperl=c(NA,(0+0.01),(100+0.01),(fixEC50+.01)))
          goodnessfit[i] <- sqrt(sum(residuals(mod1.pl)^2))
        },error=function(e){cat("Error :", conditionMessage(e),"\n")})
      }
      
      index1 <- which(goodnessfit<(min(goodnessfit)+qchisq(p=0.95, df=3))) #which of the dots of my goodnessfit are below the line = 95 % CI 
      conc[min(index1)]#get upper and lower 95% CI => put outcome to the respective data in the excel file (see line 107)
      conc[max(index1)]
      
      pl_ci <- as.data.frame(rbind(cbind(paste("EC50_PL95%CI_Lower Limit:",br,sep=""),as.character(conc[min(index1)])),
                                   cbind(paste("EC50_PL95%CI_Upper Limit:",br,sep=""),as.character(conc[max(index1)]))))
      
      ## UPDATE 08.12.2015 >> profile likelihood EC50 (the concentration corresponding to the min of goodnessfit)
      pl_ci <- rbind(pl_ci,cbind(paste("EC50_PL:",br,sep=""),as.character(conc[which(goodnessfit == min(goodnessfit))])))
      
      appendWorksheet(res.wb,pl_ci, paste(s,"_EC50",sep=""),header=FALSE)
      
      plot(conc, goodnessfit, main= paste("Profile Likelihood CI",s,":",br,sep=" "),ylim=c(0.1,1000)) 
      abline(h=min(goodnessfit), col="red")
      abline(h=min(goodnessfit)+(qchisq(p=0.95, df=3)))
 
      ## UPDATE 08.12.2015 >> add to plots: profile likelihood EC50 (the concentration corresponding to the min of goodnessfit)
      abline(v=conc[which(goodnessfit == min(goodnessfit))],col="cyan")
      text(10+conc[which(goodnessfit == min(goodnessfit))],0,
           paste("EC50_PL=",conc[which(goodnessfit == min(goodnessfit))],sep=""),col="cyan",cex=0.75)
      
    }
    dev.off()
    
    
  }
  dat.out[[b]] <- dat.out.t 
  
  saveWorkbook(res.wb)
}

