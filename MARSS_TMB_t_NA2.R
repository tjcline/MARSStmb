####This is a fairly general MARSS model to estimated using TMB
#Model and Code by Tim Cline 
#4/12/2018
#I am consistently adding functionality to the model and wrapper function to handle a wider variety of cases so check back for updated code

library(TMB)
library(adnuts)
library(snowfall)
library(Matrix)
StoreOldWD<-getwd()
setwd('~/Documents/RStuff/MARSS_TMB')
#
MARSS_TMB <- "// A General Multivariate Autoregressive State Space Model
#include <TMB.hpp>

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(obs); /*  timeSteps x stateDim*/
  DATA_MATRIX(Z);
  DATA_MATRIX(ProCovars);
  DATA_MATRIX(ObsCovars);
  DATA_INTEGER(EVCVpro);
  DATA_INTEGER(EVCVobs);
  PARAMETER_VECTOR(logsdObs);
  PARAMETER_VECTOR(cholCorrObs);
  PARAMETER_MATRIX(A);
  PARAMETER_MATRIX(D);
  PARAMETER_VECTOR(logsdPro);
  PARAMETER_VECTOR(cholCorrPro);
  PARAMETER_MATRIX(B);
  PARAMETER_MATRIX(U);
  PARAMETER_MATRIX(C);
  PARAMETER_MATRIX(states); /* State */

  int timeSteps=obs.row(0).size();
  int nTimeSeries=obs.col(0).size();
  int nStates=Z.row(0).size();

  vector<Type> sdPro=exp(logsdPro);
  vector<Type> sdObs=exp(logsdObs);

  using namespace density;
  //Here is the setup of the Process Variance values (warning its nasty)
  matrix<Type> proCorr(nStates,nStates);
  proCorr.fill(cholCorrPro(0));
  proCorr.diagonal().fill(Type(1.0));
  
  matrix<Type> vec1(nStates,1);
  vec1.fill(1.0);
  matrix<Type> FullCorrMatPro(nStates,nStates);
  FullCorrMatPro = vec1.asDiagonal();
  

  UNSTRUCTURED_CORR_t<Type> corMatGenObs(cholCorrObs);// This is the full CormatObs
  matrix<Type> FullCorrMatObs=corMatGenObs.cov();
  if(EVCVobs==1){
   	 matrix<Type> obsCorr(nTimeSeries,nTimeSeries);
 	 for(int i=0;i<nTimeSeries;i++){
 	 for(int j=0;j<nTimeSeries;j++){
 	 	if(i==j){
 			obsCorr(i,j)=Type(1.0);
 		}else{
 			obsCorr(i,j)=cholCorrObs(0); // just use one as they are all the same
 		}
 	 }}
	 MVNORM_t<Type> corMatGenObs(obsCorr);
	 matrix<Type> FullCorrMatObs=corMatGenObs.cov();
  }

  /* Define likelihood */
  Type ans=0;
  matrix<Type> predStates = B*states + C*ProCovars;
  for(int i=1;i<timeSteps;i++){
	vector<Type> predDiffer = states.col(i) - (predStates.col(i-1)+U);
    if(nStates>1){
		if(EVCVpro==1){
			MVNORM_t<Type> corMatGenPro(proCorr);
			FullCorrMatPro=corMatGenPro.cov();
			ans+= VECSCALE(MVNORM(proCorr),sdPro)(predDiffer); // Process likelihood
		}else{
	  	    UNSTRUCTURED_CORR_t<Type> corMatGenPro(cholCorrPro);// This is the full CormatPro
	  	    FullCorrMatPro=corMatGenPro.cov();
			ans+= VECSCALE(corMatGenPro,sdPro)(predDiffer); // Process likelihood
		}
	}else{
		ans -= dnorm(predDiffer(0),Type(0.0),sdPro(0),1);
	}
  }

  //Observation likelihood with multivariate normal error.
  matrix<Type> obsPred = Z*states + D*ObsCovars;
  for(int timeIndex=0;timeIndex<timeSteps;timeIndex++){
     int nonNAcount = 0; //start at zero NA values
 	 vector<int> GoodVals(nTimeSeries);
	 for(int NAIndex=0;NAIndex<nTimeSeries;NAIndex++){//loop over all time series for this time step
	    if(!isNA(obs.col(timeIndex)(NAIndex))){//if value is not NA
			GoodVals(nonNAcount) = NAIndex; //add position to good values (good values only stored in beginning of vector)
			nonNAcount++; //increment the values of
		}
	 }
	 if(nonNAcount<nTimeSeries){//if NA values present do this subsetting of the observation portion of the parameter list.
		if(nonNAcount>0){
			matrix<Type> subCorr(nonNAcount,nonNAcount);
			vector<Type> subSds(nonNAcount);
			vector<Type> subData(nonNAcount);
	 		vector<Type> subPred(nonNAcount);
			vector<Type> subA(nonNAcount);

			for(int subIndex=0;subIndex<nonNAcount;subIndex++){
				subData(subIndex) = obs.col(timeIndex)(GoodVals(subIndex));
				subPred(subIndex) = obsPred.col(timeIndex)(GoodVals(subIndex));
				subA(subIndex) = A.col(0)(GoodVals(subIndex));
				subSds(subIndex) = sdObs(GoodVals(subIndex));
				for(int k=0; k<nonNAcount;k++){
					subCorr(subIndex,k) = FullCorrMatObs(GoodVals(subIndex),GoodVals(k));
				}
			}
  			vector<Type> subDiffer(nonNAcount);
			subDiffer = subData-(subPred+subA);
			ans += VECSCALE(MVNORM(subCorr),subSds)(subDiffer);
		}
 	}else{ //if no NA values present, compute a normal observation likelihood
   		vector<Type> obsDiffer = obs.col(timeIndex)-(obsPred.col(timeIndex)+A);
 		ans += VECSCALE(corMatGenObs,sdObs)(obsDiffer);
 	}//end of data likelihood for this time step
 }//end of loop over time steps

  //REPORTING

  matrix<Type> R(nTimeSeries,nTimeSeries);
  matrix<Type> doSD(nTimeSeries,1);
  doSD = sdObs;
  R = doSD.asDiagonal() * FullCorrMatObs * doSD.asDiagonal();

  matrix<Type> Q(nStates,nStates);
  matrix<Type> dpSD(nStates,1);
  dpSD = sdPro;
  Q = dpSD.asDiagonal() * FullCorrMatPro * dpSD.asDiagonal();

  REPORT(B);
  REPORT(U);
  REPORT(C);
  REPORT(Q);
  REPORT(A);
  REPORT(D);
  REPORT(R);
  REPORT(states);

  ADREPORT(B);
  ADREPORT(U);
  ADREPORT(C);
  ADREPORT(Q);
  ADREPORT(A);
  ADREPORT(D);
  ADREPORT(R);
  ADREPORT(states);

  return ans;
}"
#A nasty bit of code that basically checks if the the model above is EXACTLY the same as the already compiled version. Saves a small amount of time for rerunning code.
if(file.exists('MARSS_TMB.cpp')){
  t1<-read.table('MARSS_TMB.cpp',sep="\t",stringsAsFactors=F,comment.char='', quote = "")
  write(MARSS_TMB, file = "MARSS_TMB_temp.cpp")
  t2<-read.table('MARSS_TMB_temp.cpp',sep="\t",stringsAsFactors=F,comment.char='', quote = "")
  file.remove('MARSS_TMB_temp.cpp')
  #t1<-data.frame(t1[-which(t1==''),])
  #t2<-data.frame(strsplit(dfa_model,fixed=T,'\n'))
  #t2<-data.frame(t2[-which(t2==''),])
  if(identical(t1,t2)){
    dyn.load(dynlib("MARSS_TMB"))
  }else{
    write(MARSS_TMB, file = "MARSS_TMB.cpp")
    compile("MARSS_TMB.cpp")
    dyn.load(dynlib("MARSS_TMB"))
  }
  }else{
    write(MARSS_TMB, file = "MARSS_TMB.cpp")
    compile("MARSS_TMB.cpp")
    dyn.load(dynlib("MARSS_TMB"))
}
dyn.load(dynlib("MARSS_TMB"))
#Wrapper function
MARSStmb<-function(obsIn,model=NULL,map=NULL,control=list(EstSE=TRUE)){
	#obsIn<-ZL_OA3
	#model=list(Z=Z.1,B='diagonal and unequal',R='diagonal and equal',Q='diagonal and unequal',A='zero',U='unequal')	
	#map<-NULL
	#Check parameter Inputs, fill as necessary, and convert text to matrix
	parList<-c('Z','D','d','A','R','B','C','c','U','Q')
	defaults<-list(Z='unequal',D='zero','d'=matrix(0,nrow=1,ncol=ncol(obsIn)),A='zero',R='diagonal and equal',
		B='rw',C='zero','c'=matrix(0,nrow=1,ncol=ncol(obsIn)),U='zero',Q='diagonal and unequal')
	if(is.null(model)){
		model<-defaults
	}else{
		#Fill in missing parameters with defaults
		for(p in 1:length(parList)){
			if(!(parList[p] %in% names(model))){
				model[[parList[p]]]<-defaults[[parList[p]]]
			}
		}
	}
	
	#initialize a blank list for map if not entered
	if(is.null(map)) map <- list()
	#if(is.null(control)) control <-list(EstSE=FALSE)
	
	#Zmatrix stuff
	if(identical(model$Z,'unequal')){ model$Z<-diag(nrow(obsIn)) #default 1 state per series
		}else if(identical(model$Z,'equal')){ model$Z<-matrix(1,nrow(obsIn),1)}
	
	#Dmatrix stuff (observation covariates)
	if(identical(model$D,'zero')){
		model$D<-matrix(0,nrow(obsIn),nrow(model$d))
		map$D<-as.factor(matrix(NA,nrow(obsIn),nrow(model$d)))
		}else if(identical(model$D,'equal')){
			model$D<-matrix(0.5,nrow(obsIn),nrow(model$d))
			map$D<-as.factor(matrix(sort(rep(seq(1,nrow(model$d)),nrow(obsIn))),nrow(obsIn),nrow(model$d)))	
		}else if(identical(model$D,'unequal')){
			model$D<-matrix(0.5,nrow(obsIn),nrow(model$d))
			map$D<-as.factor(matrix(seq(1,nrow(obsIn)*nrow(model$d)),nrow(obsIn),nrow(model$d)))	
		}else if(identical(model$D,'individual')){
			model$D<-matrix(0.5,nrow(obsIn),nrow(model$d))
			map$D<-as.factor(matrix(seq(1,nrow(obsIn)*nrow(model$d)),nrow(obsIn),nrow(model$d)))
		}else{
			if(!('D' %in% names(map))){
				print("Need to specify 'D' in map");break
			}
	}
		
	#Cmatrix stuff	
	if(identical(model$C,'zero')){
		model$C<-matrix(0,ncol(model$Z),nrow(model$c))
		map$C<-as.factor(matrix(NA,ncol(model$Z),nrow(model$c)))
		}else if(identical(model$C,'equal')){
			model$C<-matrix(0.5,ncol(model$Z),nrow(model$c))
			map$C<-as.factor(matrix(sort(rep(seq(1,nrow(model$c)),ncol(model$Z))),ncol(model$Z),nrow(model$c)))	
		}else if(identical(model$C,'unequal')){
			model$C<-matrix(0.5,ncol(model$Z),nrow(model$c))
			map$C<-as.factor(matrix(seq(1,ncol(model$Z)*nrow(model$c)),ncol(model$Z),nrow(model$c)))	
		}else if(identical(model$C,'individual')){
			model$C<-matrix(0.5,ncol(model$Z),nrow(model$c))
			map$C<-as.factor(matrix(seq(1,ncol(model$Z)*nrow(model$c)),ncol(model$Z),nrow(model$c)))
		}else{
			if(!('C' %in% names(map))){
				print("Need to specify 'C' in map");break
			}
	}
	
	if(identical(model$A,'unequal')){
		model$A<-matrix(0,nrow(obsIn),1)
		map$A<-as.factor(matrix(seq(1,nrow(obsIn)),nrow(obsIn),1))
		}else if(identical(model$A,'equal')){
			model$A<-matrix(0,nrow(obsIn),1)
			map$A<-as.factor(matrix(1,nrow(obsIn),1))
		}else if(identical(model$A,'zero')){
			model$A<-matrix(0,nrow(obsIn),1)
			map$A<-as.factor(matrix(NA,nrow(obsIn),1))
		}else{
			if(!('A' %in% names(map))){
				print("Need to specify 'A' in map");break
			}
	}
		
	if(identical(model$U,'unequal')){
		model$U<-matrix(0,ncol(model$Z),1)
		map$U<-as.factor(matrix(seq(1,ncol(model$Z)),ncol(model$Z),1))
		}else if(identical(model$U,'equal')){
			model$U<-matrix(0,ncol(model$Z),1)
			map$U<-as.factor(matrix(1,ncol(model$Z),1))
		}else if(identical(model$U,'zero')){
			model$U<-matrix(0,ncol(model$Z),1)
			map$U<-as.factor(matrix(NA,ncol(model$Z),1))
		}else{
			if(!('U' %in% names(map))){
				print("Need to specify 'U' in map");break
			}
	}
		
	if(identical(model$B,'rw')){
		model$B<-diag(ncol(model$Z))
		map$B<-as.factor(matrix(NA,ncol(model$Z),ncol(model$Z)))
		}else if(identical(model$B,'diagonal and unequal')){
			model$B<-diag(0.5,ncol(model$Z))
			map$B<-matrix(NA,ncol(model$Z),ncol(model$Z))
			diag(map$B)<-seq(1,ncol(model$Z))
			map$B<-as.factor(map$B)
		}else if(identical(model$B,'unconstrained')){
			model$B<-diag(ncol(model$Z))
			map$B<-as.factor(matrix(seq(1,ncol(model$Z)^2),ncol(model$Z),ncol(model$Z)))
		}else{
			if(!('B' %in% names(map))){
				print("Need to specify 'B' in map");break
			}
	}
		
	if(identical(model$R,'diagonal and equal')){
		model$cholCorrObs<-rep(0,nrow(obsIn)*(nrow(obsIn)-1)/2)
		model$logsdObs<-log(rep(0.1,nrow(obsIn)))
		map$logsdObs<-as.factor(rep(1,nrow(obsIn)))
		map$cholCorrObs <-as.factor(rep(NA,nrow(obsIn)*(nrow(obsIn)-1)/2))
		}else if(identical(model$R,'diagonal and unequal')){
			model$cholCorrObs<-rep(0,nrow(obsIn)*(nrow(obsIn)-1)/2)
			model$logsdObs<-log(rep(0.5,nrow(obsIn)))
			map$logsdObs<-as.factor(seq(1,nrow(obsIn)))
			map$cholCorrObs <-as.factor(rep(NA,nrow(obsIn)*(nrow(obsIn)-1)/2))
		}else if(identical(model$R,'unconstrained')){
			model$cholCorrObs<-rep(0,nrow(obsIn)*(nrow(obsIn)-1)/2)
			model$logsdObs<-log(rep(0.5,nrow(obsIn)))
			map$logsdObs<-as.factor(seq(1,nrow(obsIn)))
			map$cholCorrObs <-as.factor(seq(1,nrow(obsIn)*(nrow(obsIn)-1)/2))
		}else if(identical(model$R,'zero')){
			model$cholCorrObs<-rep(0,nrow(obsIn)*(nrow(obsIn)-1)/2)
			model$logsdObs<-log(rep(1e-6,nrow(obsIn)))
			map$logsdObs<-as.factor(rep(NA,nrow(obsIn)))
			map$cholCorrObs <-as.factor(rep(NA,nrow(obsIn)*(nrow(obsIn)-1)/2))
		}else if(identical(model$R,'equalvarcov')){
			model$cholCorrObs<-rep(0,nrow(obsIn)*(nrow(obsIn)-1)/2)
			model$logsdObs<-log(rep(0.1,nrow(obsIn)))
			map$logsdObs<-as.factor(rep(1,nrow(obsIn)))
			map$cholCorrObs <-as.factor(rep(1,nrow(obsIn)*(nrow(obsIn)-1)/2))
		}else{
			if(!('R' %in% names(map))){
				print("Need to specify 'R' in map");break
			}
	}
		
	if(identical(model$Q,'diagonal and equal')){
		if(ncol(model$Z)>1){
			model$cholCorrPro<-rep(0,ncol(model$Z)*(ncol(model$Z)-1)/2)
			map$cholCorrPro <-as.factor(rep(NA,ncol(model$Z)*(ncol(model$Z)-1)/2))
		}else{
			model$cholCorrPro<-matrix(1,1,1)
			map$cholCorrPro<-as.factor(matrix(NA,1,1))
		}
		model$logsdPro<-log(rep(0.5,ncol(model$Z)))
		map$logsdPro<-as.factor(rep(1,ncol(model$Z)))
	}else if(identical(model$Q,'diagonal and unequal')){
			if(ncol(model$Z)>1){
				model$cholCorrPro<-rep(0,ncol(model$Z)*(ncol(model$Z)-1)/2)
				map$cholCorrPro <-as.factor(rep(NA,ncol(model$Z)*(ncol(model$Z)-1)/2))
			}else{
				model$cholCorrPro<-matrix(1,1,1)
				map$cholCorrPro<-as.factor(matrix(NA,1,1))
			}
			model$logsdPro<-log(rep(0.5,ncol(model$Z)))
			map$logsdPro<-as.factor(seq(1,ncol(model$Z)))
		}else if(identical(model$Q,'unconstrained')){
			if(ncol(model$Z)>1){
				model$cholCorrPro<-rep(0,ncol(model$Z)*(ncol(model$Z)-1)/2)
				map$cholCorrPro <-as.factor(seq(1,ncol(model$Z)*(ncol(model$Z)-1)/2))
			}else{
				model$cholCorrPro<-matrix(1,1,1)
				map$cholCorrPro<-as.factor(matrix(NA,1,1))
			}
			model$logsdPro<-log(rep(0.5,ncol(model$Z)))
			map$logsdPro<-as.factor(seq(1,ncol(model$Z)))
		}else if(identical(model$Q,'zero')){
			if(ncol(model$Z)>1){
				model$cholCorrPro<-rep(0,ncol(model$Z)*(ncol(model$Z)-1)/2)
				map$cholCorrPro <-as.factor(rep(NA,ncol(model$Z)*(ncol(model$Z)-1)/2))
			}else{
				model$cholCorrPro<-matrix(1,1,1)
				map$cholCorrPro<-as.factor(matrix(NA,1,1))
			}
			model$logsdPro<-log(rep(0.5,ncol(model$Z)))
			map$logsdPro<-as.factor(rep(NA,ncol(model$Z)))
		}else if(identical(model$Q,'equalvarcov')){
			if(ncol(model$Z)>1){
				model$cholCorrPro<-rep(0,ncol(model$Z)*(ncol(model$Z)-1)/2)
				map$cholCorrPro <-as.factor(rep(1,ncol(model$Z)*(ncol(model$Z)-1)/2))
			}else{
				model$cholCorrPro<-matrix(1,1,1)
				map$cholCorrPro<-as.factor(matrix(NA,1,1))
			}
			model$logsdPro<-log(rep(0.5,ncol(model$Z)))
			map$logsdPro<-as.factor(rep(1,ncol(model$Z)))
		}else{
			if(!('Q' %in% names(map))){
				print("Need to specify 'Q' in map");break
			}
	}
	
	data <- list(obs=obsIn,
		Z=model$Z,
		ProCovars=model$c,
		ObsCovars=model$d)
	data$EVCVobs<-ifelse(identical(model$R,'equalvarcov'),1,0)
	data$EVCVpro<-ifelse(identical(model$Q,'equalvarcov'),1,0)
	#Creates the input parameter list
	parameters <- model[c('logsdObs','cholCorrObs','A','D','logsdPro','cholCorrPro','B','U','C')]
	map<-map[c('logsdObs','cholCorrObs','A','D','logsdPro','cholCorrPro','B','U','C')]
	parameters$states <- matrix(0,ncol=ncol(obsIn),nrow=ncol(model$Z))
	
	
	#Creates the model object and runs the optimization
	obj1 <- MakeADFun(data,parameters,random="states",DLL="MARSS_TMB",silent=T,map=map)
	obj.mcmc <- MakeADFun(data,parameters,DLL="MARSS_TMB",silent=T,map=map)
	if(sum(!(is.na(map$B)))>0){
		lowerB<-matrix(-Inf,ncol(model$B),ncol(model$B));diag(lowerB)<--1
		lowerB<-as.vector(lowerB);lowerB<-lowerB[!is.na(map$B)]
		#names(lowerB)<-rep('B',length(lowerB))
		upperB<-matrix(Inf,ncol(model$B),ncol(model$B));diag(upperB)<-1
		upperB<-as.vector(upperB);upperB<-upperB[!is.na(map$B)]
		#names(upperB)<-rep('B',length(upperB))
		
		
		uniMap<-lapply(map,FUN=function(x) unique(na.omit(x)))
		lowerBounds<-c(rep(-Inf,length(as.numeric(unlist(uniMap[!(names(uniMap) %in% c('B','U','C'))])))),lowerB,rep(-Inf,length(as.numeric(unlist(uniMap[names(uniMap) %in% c('U','C')])))))
		upperBounds<-c(rep(Inf,length(as.numeric(unlist(uniMap[!(names(uniMap) %in% c('B','U','C'))])))),upperB,rep(Inf,length(as.numeric(unlist(uniMap[names(uniMap) %in% c('U','C')])))))
		
		opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr,lower=lowerBounds,upper=upperBounds,control=list(iter.max=2,eval.max=2))
		obj1$control=list(trace=1,REPORT=1,reltol=1e-12,maxit=4000)
		obj1$fn()
		obj1$gr()
		obj1$method="L-BFGS-B"
		obj1$par=opt1$par
		obj1$lower=lowerBounds
		obj1$upper=upperBounds
		system.time(opt1 <- do.call("optim",obj1))
		
		#profi<-tmbprofile(obj1,2,ytol=4,ystep=0.01)
		#plot(profi$value~profi$logsdPro)
		#mcmc1<-sample_tmb(obj1,iter=1000,lower=lowerBounds,upper=upperBounds,init=NULL,parallel=T,cores=3,path=getwd(),laplace=TRUE)
		#mcmc1<-sample_tmb(obj1,iter=1000,lower=lowerBounds,upper=upperBounds,init=NULL,laplace=T)
		
		#e1<-extract_samples(mcmc1)
		#launch_shinytmb(mcmc1)
	}else{
			opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr,control=list(iter.max=2000,eval.max=2000))
			#newtonOption(obj1,smartsearch=TRUE)
			obj1$control=list(trace=1,REPORT=1,reltol=1e-12,maxit=2000)
			obj1$fn()
			obj1$gr()
			obj1$method='BFGS'
			obj1$par=opt1$par
			obj1$lower
			system.time(opt1 <- do.call("optim",obj1))
	}
	
	#opt1 <- do.call("optim",obj1) #opt1 contained the optimization information like objective_function_value. Convergence criteria.
	#lowerPar<-c(rep(-1000,ncol(obs)),rep(-100,ncol(obs)*(ncol(obs)-1)/2),rep(-2,NumStates^2),rep(-10,nrow(obs)),rep(-1,nrow(obs)))
	#upperPar<-c(rep(1000,ncol(obs)),rep(100,ncol(obs)*(ncol(obs)-1)/2),rep(2,NumStates^2),rep(10,nrow(obs)),rep(1,nrow(obs)))
	#opt1 <- nlminb(obj1$par,obj1$fn,obj1$gr,lower=lowerPar,upper=upperPar)
	#opt1 <- nlminb(opt1$par,obj1$fn,obj1$gr)
	pl1 <- obj1$env$parList()#This contains all of your parameter estimates RAW as they come out of the optimizer
	mlePars<-obj1$report()
	if(control$EstSE){
		sdr<-sdreport(obj1)
		SES<-list()
		SEnames<-unique(names(sdr$value))
		for(i in 1:length(SEnames)){
			SEindex<-which(names(sdr$value) == SEnames[i])
			if(SEnames[i] %in% c('B','Q')){
				SES[[SEnames[i]]]<- matrix(sdr$sd[SEindex],nrow(model$U),nrow(model$U))
			}else if(SEnames[i] %in% c('U','C')){
				SES[[SEnames[i]]]<- matrix(sdr$sd[SEindex],nrow(model$U))
			}else if(SEnames[i] %in% c('A','D')){
				SES[[SEnames[i]]]<- matrix(sdr$sd[SEindex],nrow(model$A))
			}else if(SEnames[i] %in% c('R')){
				SES[[SEnames[i]]]<- matrix(sdr$sd[SEindex],nrow(model$A),nrow(model$A))
			}else if(SEnames[i] %in% c('states')){
				SES[[SEnames[i]]]<- matrix(sdr$sd[SEindex],nrow(model$U))
			}
		}
	}
	
	#if(EstSE){pl1$R <- matrix(sdr$value[which(names(sdr$value)=='FullCovMat')],nrow=length(logsdObs),ncol=length(logsdObs))
		#}else{pl1$R<- diag(exp(pl1$logsdObs)) %*% obj1$report()$FullCorrMat  %*% diag(exp(pl1$logsdObs))}
	
	#Fits for each time series
	FitSeries<- data$Z %*% mlePars$states + mlePars$D %*% data$ObsCovars + matrix(mlePars$A,nrow(data$Z),ncol(mlePars$states)) #matrix(0,nrow=1,ncol=ncol(obsIn))
	
	#Compute AIC.
	k<-length(opt1$par)
	n<-length(na.omit(as.vector(data$obs)))
	if(sum(!(is.na(map$B)))>0){
		AIC<-2*k + 2*opt1$value;AIC
		AICc <- AIC + (2*k^2+2*k)/(n-k-1)
		BIC <- log(n)*k +2*opt1$value
	}else{
		AIC<-2*k + 2*opt1$value;AIC
		AICc <- AIC + (2*k^2+2*k)/(n-k-1)
		BIC <- log(n)*k +2*opt1$value
	}
	
	
	#print(AIC)
	# if(EstSE){return(list(Optimization = opt1, Estimates = pl1, Fits = FitSeries,AIC=AIC,StdErr=SES,ParCorrs=sdr$cov.fixed))
# 	}else{return(list(Optimization = opt1, Estimates = mlePars, Fits = FitSeries,AIC=AIC))}
	mlePars$Optimization<- opt1
	mlePars$Fits <- FitSeries
	mlePars$AIC <- AIC
	mlePars$AICc <- AICc
	mlePars$BIC <- BIC
	mlePars$SES <- SES
	if(mlePars$Optimization$convergence != 0) print('Warning: Potential optimization error')
	return(mlePars)	
}


setwd(StoreOldWD)


