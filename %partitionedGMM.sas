/*****************************************************************/
/**	SAS Code: %partitionedGMM								    **/
/** Programmer: Kyle Irimata									**/
/** Description: Performs partitioned GMM regression			**/
/*****************************************************************/

%macro partitionedGMM(ds=., file=, timeVar=, outVar=, predVarTD=, idVar=, alpha=0.05, predVarTI=., distr=bin, optim=NLPCG, MC=LWY);

*Check if there is a library;
%if &ds. ^=. %then %do;
	LIBNAME DS &ds.;  
	DATA mydata; 
		SET DS.&file.; 
	RUN;
%end;
%else %do;
	DATA mydata;
		SET &file.;
	RUN;
%end;

TITLE "Partitioned GMM";
PROC SORT DATA=mydata OUT=mydatasorted; 
BY &timeVar.; RUN;

*Check if there is a time independent variable;
%if &predVarTI. ^=.  %then %do;
	%let predVar = &predVarTD &predVarTI;
%end;
%else %do;
	%let predVar = &predVarTD;
%end;

*Use either the LWY or LS moment check;
%if &MC=LWY %then %do; *LWY approach;

*Obtain residuals;
%if &distr.=normal %then %do;
	proc reg data=mydatasorted NOPRINT;
	BY &timeVar.;
	MODEL &outVar. = &predVar.;
	OUTPUT OUT=outpool3 PREDICTED=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = 1/mu; 
	rsdraw = &outVar.-mu; 
%end;
%else %if &distr.=bin %then %do;
	PROC logistic DATA=mydatasorted NOPRINT; 
	BY &timeVar.;
	MODEL &outVar. (event='1') = &predVar. / aggregate scale=none;
	OUTPUT OUT=outpool3 P=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = mu*(1-mu);
	rsdraw = &outVar.-mu; 
%end;
%else %do;
	%put ERROR must be normal or binomial distributions; 
	%return;
%end;

PROC SORT DATA=outpool3 OUT=outpool3 ;
  BY &idVar. &timeVar.; RUN;
quit;

PROC IML;
use outpool3;                                                           
%if &predVarTI. ^=. %then %do;
	read all VARIABLES {&predVarTD. &predVarTI.} into Zmat; 
	read all var {&predVarTI.} into timeInd;
	read all var {&predVarTD.} into timeDep;
%end;
%else %do;
	read all VARIABLES {&predVarTD.} into Zmat;
%end;
read all var {wt} into wt;
read all var {rsdraw} into rsd;
read all var {&idVar.} INTO ID;
read all var {&timeVar.} INTO time;
close outpool3;

N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
Np = ncol(Zmat); *Number of predictors, not including intercept;

%if &predVarTI. ^=. %then %do;
	NpTI = ncol(timeInd);*The last NPTI columns of Zmat will be time independent;
	NpTD = ncol(timeDep); *The first NPTD columns of Zmat will be time dependent;
%end;

*Find the correlation between X (a) and Y (rsd) across the T time points;
start rho(a,rsd) global(N,T);
abm = j(N,2*T,.);
abm[,1:T] = shape(rsd,N);	* N x T - First T columns are the values of rsd sorted into T columns (by time);
abm[,T+1:2*T] = shape(a,N); * Remaining T columns are the values of a sorted into T columns;
corr = corr(abm);  
rho = corr[1:T,T+1:2*T];    * T x T - Only take the correlations between the X and Y in the T time points;
return(rho);
finish rho;


*Standard deviation for each correlation?;
start stddev(a,rsd) global(N,T);
bm = shape(rsd,N);    		 * N x T;
bdev = bm-j(N,1,1)*bm[:,];   * bdev N x T,   bm[:,] Col Mean is a row vector   1 x T - each row of bm minus column means;
bdev2 = bdev#bdev;      
am = shape(a,N);   
adev = am-j(N,1,1)*am[:,];  *N x T;
adev2 = adev#adev;      
stddev = sqrt( (1/N)*t(bdev2)*adev2 );   * T x T;
return(stddev);
finish stddev;

* corrected standardization;
start stdzn(x) global(N,T);
xrows = shape(x,N);   *N x T - by shape default columns are T=nrow(x)/N;
y = xrows - xrows[:,];  *N x T - Each value minus the column mean (mean for that time point);
vcv = (1/(N-1))*t(y)*y; * T x T;
v = diag(vcv); * T x T diagonal elements of vcv;
sinv = sqrt(inv(v));
x2 = y*sinv;   *N x T;
x2  = shape(x2,N*T,1); *N*T x 1 vector of standardized values;
return(x2);
finish stdzn;

pvec = j(Np*T*T,1,.);sevec = j(Np*T*T,1,.);  * pvec   (T*T) x Np;
r4out = j(T,T*Np,.); se4out = j(T,T*Np,.); z4out = j(T,T*Np,.); p4out = j(T,T*Np,.); 

y = rsd;
y_std = stdzn(y);

DO i=1 TO Np;
x = wt#Zmat[,i]; 			 * (N*T) x 1;
x_std = stdzn(x);

*Find p-values for the correlation of X_i and Y;
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;

*Fill the corresponding T columns of each matrix with the calculated values;
r4out[,T*(i-1)+1:T*i] = r;
se4out[,T*(i-1)+1:T*i] = se;
z4out[,T*(i-1)+1:T*i] = z;
p4out[,T*(i-1)+1:T*i] = p;

DO j = 1 TO T;
p[j,j] = 1;  * Not going to test diagonal elements, set to 1;
END;
*Takes the values in p4out and se4out, across row, then down column and creates vectors;
pvec[T*T*(i-1)+1:T*T*i,1] = shape(p,T*T,1);    *(T*T) x 1;
sevec[T*T*(i-1)+1:T*T*i,1] = shape(se,T*T,1);   *(T*T) x 1;
END;


TypeVec2 = (pvec >= &alpha.*j(Np*T*T,1,1) );
Type2 = shape(TypeVec2, Np, T*T);

*Individual test approach;
x = wt; 			 * (N*T) x 1;
x_std = stdzn(x);
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;
T_wt = p>&alpha.;
T_wt = shape(T_wt,1); 
TypeMtx3 = j(Np+1,T*T,.);
TypeMtx3[1,] = T_wt;
TypeMtx3[2:Np+1,] = Type2;
Type2[1,] = shape(I(T),1);

*Adjust the last NpTI rows of TypeMtx3 to account for time independent covariates;
%if &predVarTI. ^=. %then %do;
	TypeMtx3= TypeMtx3[1:(NPTD+1),]; *Subset TypeMtx3 to include only time dependent;
	TypeMtxIND = repeat(shape(I(T),1),NPTI); *Create rows of typeMtx for the time independent;
	Np = NpTD;
%end;

/* define helper functions ROW and COL */
start row(x);  /* return matrix m such that m[i,j] = i */
   return( repeat( T(1:x), 1, x ));
finish;
start col(x);  /* return matrix m such that m[i,j] = j */
   return( repeat(1:x, x) );
finish;

*Helper matrices;
rt = row(T);
ct = col(T);

*Indices of upper diagonal for removal;
upperIdx = loc(ct>rt);

*Remove backwards conditions;
typeMtxNB = TypeMtx3;
typeMtxNB[,upperIdx] = 0;

*Create the combination TypeMtx, of maximum size to be pared down later, excluding the intercept;
TypeMtxCombo = j((Np*T),T*T,0); 

TypeMtxCombo2 = j((Np*T),T*T,0);
*Loop across the different types, starting with current;
DO i = 0 TO (T-1);
	*Identify the indices for the appropriate setting;
	Idx = loc(rt-ct = i);
	IdxDif = setdif(1:T*T, Idx);
	temp = typeMtxNB[2:(Np+1),];
	temp[,IdxDif] = 0;
	TypeMtxCombo[(i*Np+1):((i+1)*Np),] = temp;

	*Shift the values to accomodate shifted Y-X relationship;
	if(i>0) then DO; 
		dummyZero = j(Np, i, 0) || temp;
		TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = dummyZero[,1:(T*T)];
	END;
	
	ELSE TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = temp;

END; 

*Identify predictors or relationships with no valid moment conditions;
neq = TypeMtxCombo2[,+];
*Print a note about omitted moment conditions;
if min(neq)=0 then do;
	zeroPred = loc(neq=0);
	Note = "There are no valid moment conditions for " +char(ncol(zeroPred))+" covariate relationship(s).";
	print Note[label="Moment Condition Notes"];
	print "These covariate relationships will be omitted in the analysis.";
end;
else print "All covariate relationships will be evaluated."[label="Moment Condition Notes"];

keepPred = loc(neq>0);

TypeMtxCombo3 = TypeMtxCombo2[keepPred,];


*Append the intercept conditions to the matrix;
*Append the time independent conditions to the matrix, if available;
%if &predVarTI. ^=. %then %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxIND // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxIND // TypeMtxCombo3; *Type I intercept;
	
%end;
%else %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxCombo3; *Type I intercept;
%end;

*If there are time independent variables, they are represented after the intercept, before the time dependent;
CREATE TypeMtxCombo3 from TypeMtxCombo3;
APPEND from TypeMtxCombo3;
CLOSE TypeMtxCombo3;
print "Each row of TypeMtx is the shifted type vector for each of the predictors, by individual test";
print TypeMtxCombo3;

*Create a printable form of the type matrix;
Type2[,upperIdx] = 0;
Type4out = j(T,T*(Np),.);
DO i=1 to Np;
Type4out[,T*(i-1)+1:T*i] = shape(Type2[i,],T,T);
END;

CREATE Type4out from Type4out;
APPEND from Type4out;
CLOSE Type4out;

*Create the modified data set;
USE mydatasorted;
read all var {&predVarTD.} into xnew;
%if &predVarTI. ^=. %then %do;
	read all VAR {&idVar. &outVar. &timeVar. &predVarTI.} INTO otherVars;
%end;
%else %do;
	read all VAR {&idVar. &outVar. &timeVar.} INTO otherVars;
%end;
CLOSE mydatasorted;

*Creating the new  variables;
X2 = j(N*T,Np*(T-1),0);
DO i=1 TO T-1;
	X2[(i*N+1):N*T, ((i-1)*Np+1):i*Np] = Xnew[1:(T-i)*N,];
END;

X2 = Xnew || X2;

*Remove the predictors with no valid moments;
X3 = X2[,keepPred];

*Add in the ID, outcome and time (and possible time independent variables;
X3 = otherVars || X3;

*Create the adjusted variable names;
predNames = t(repeat(t({&predVarTD.}),T));
lagName = j(1,Np*T,0);
DO i=1 to (T-1);
	lagName[1,(i*Np+1):(i+1)*Np]=j(1,Np,i);
END;
varNames = catx("_", predNames, char(lagName));

*Retain only the variable names for those with valid moment conditions;
varNames = varNames[,keepPred];
%if &predVarTI. ^=. %then %do;
	varNames = {&predVarTI.} ||  varNames;
%end;

*Add in the names for the ID, outcome and time;
varNames2 = {&idVar.} || {&outVar.} || {&timeVar.} ||  varNames;

*Sort the data;
call sort(X3 , {1 3});

CREATE Mydata3 from X3[c=varNames2];
APPEND from X3;
close Mydata3;

*Module as demonstrated at http://blogs.sas.com/content/iml/2016/01/18/create-macro-list-values.html;
start CreateMacro(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   /* convert numeric to character */
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       /* delimit and concatenate */
   s = substr(s, 1, nleng(s)-nleng(delimiter)); /* remove delimiter at end */
   call symputx(macroName, s);      /* create macro variable */
finish;

call CreateMacro(varNames, "adjPredVar");

QUIT;

%end;

%else %if &MC=LS %then %do; *Lai and Small Approach; 

*Create the lagged data sets;
proc iml;
use mydata;                                                           
%if &predVarTI. ^=. %then %do;
	read all VARIABLES {&predVarTD. &predVarTI.} into Zmat; 
	read all var {&predVarTI.} into timeInd;
	read all var {&predVarTD.} into timeDep;
%end;
%else %do;
	read all VARIABLES {&predVarTD.} into Zmat;
%end;
read all var {&idVar.} INTO ID;
read all var {&timeVar.} INTO time;
close mydata;

USE mydatasorted;
read all var {&predVarTD.} into xnew;
%if &predVarTI. ^=. %then %do;
	read all VAR {&idVar. &outVar. &timeVar. &predVarTI.} INTO otherVars;
%end;
%else %do;
	read all VAR {&idvar. &outvar. &timevar.} INTO otherVars;
%end;
CLOSE mydatasorted;

N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
Np = ncol(Zmat); *Number of predictors, not including intercept;


%if &predVarTI. ^=. %then %do;
	NpTI = ncol(timeInd);*The last NPTI columns of Zmat will be time independent;
	NpTD = ncol(timeDep); *The first NPTD columns of Zmat will be time dependent;
%end;
%else %do;
	NpTD = Np;
%end;

*Module as demonstrated at http://blogs.sas.com/content/iml/2016/01/18/create-macro-list-values.html;
*Creates modified macro variables that contain the shifted covariate names;
start CreateMacro(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   /* convert numeric to character */
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       /* delimit and concatenate */
   s = substr(s, 1, nleng(s)-nleng(delimiter)); /* remove delimiter at end */
   call symputx(macroName, s);      /* create macro variable */
finish;


*Creating the modified data;
X2 = j(N*T,NpTD*(T-1),0);
DO i=1 TO T-1;
	X2[(i*N+1):N*T, ((i-1)*NpTD+1):i*NpTD] = Xnew[1:(T-i)*N,];
END;

X2 = Xnew || X2;


*Add in the ID, outcome and time (and possible time independent variables;
X3 = otherVars || X2;

*Create the adjusted variable names;
predNames = t(repeat(t({&predVarTD.}),T));
lagName = j(1,NpTD*T,0);
DO i=1 to (T-1);
	lagName[1,(i*NpTD+1):(i+1)*NpTD]=j(1,NpTD,i);
END;
varNames = catx("_", predNames, char(lagName));

*Add in the time-independent variable names;
%if &predVarTI. ^=. %then %do;
	varNames = {&predVarTI.} ||  varNames;
%end;


*Add in the names for the ID, outcome and time;
varNames2 = {&idVar.} || {&outVar.} || {&timeVar.} ||  varNames;

*Sort the data;
call sort(X3 , {1 3});

*Create an adjusted SAS data set;
CREATE mydata3 from X3[c=varNames2];
APPEND from X3;
close mydata3;

*Save the adjusted covariate variable names as a macro variable;
call CreateMacro(varNames, "adjPredVar");

*Create the type matrix with no backwards and all others assumed valid;
/* define helper functions ROW and COL */
start row(x);  /* return matrix m such that m[i,j] = i */
   return( repeat( T(1:x), 1, x ));
finish;
start col(x);  /* return matrix m such that m[i,j] = j */
   return( repeat(1:x, x) );
finish;

*Helper matrices;
rt = row(T);
ct = col(T);

*Indices of upper diagonal for removal;
upperIdx = loc(ct>rt);
TypeMtx = j(Np,T*T,1);
TypeMtx[,upperIdx] = 0;
*Make the intercept time independent;
TypeMtx[1,] =  shape(I(T),1);

*Adjust the NpTI rows of TypeMtx to account for time independent covariates;
%if &predVarTI. ^=. %then %do;
	TypeMtx[2:NpTI+1,] = repeat(shape(I(T),1),NPTI); *Create rows of typeMtx for the time independent;
%end;
*Create the lagged TypeMtxCombo3;



*Create the lagged covariate matrix;
%if &predVarTI. ^=. %then %do; *TypeMatrix for data with time-independent;
	TypeMtxCombo3 = repeat(shape(I(T),1),Np+1); *Base typeMtxCombo3, which includes time-independent to be built on;
	modMoment = I(T);
	DO i=1 to (T-1); *Modify the lagged moments;
		modMoment[,i] = 0;
		TypeMtxCombo3 = TypeMtxCombo3 // repeat(shape(modMoment,1),NPTD);
	END;

%end;
%else %do; *TypeMatrix for data with time-dependent only;
	TypeMtxCombo3 = repeat(shape(I(T),1),Np+1); *Base typeMtxCombo3 to be built on;
	modMoment = I(T);
	DO i=1 to (T-1); *Modify the lagged moments;
		modMoment[,i] = 0;
		TypeMtxCombo3 = TypeMtxCombo3 // repeat(shape(modMoment,1),Np);
	END;
%end;

*Save the typeMtx as a dataset;
CREATE TypeMtxCombo3 from TypeMtxCombo3;
APPEND from TypeMtxCombo3;
CLOSE TypeMtxCombo3;
print TypeMtxCombo3;


quit;

%end;

%else %do;
	%put ERROR Must select either LWY or LS moment selection; 
	%return;
%end;



*Obtain initial GEE estimates for optimization;
ods select none;
proc genmod data=mydata3 descend;
	class &idVar. &timeVar.;
	model &outVar. = &adjPredVar. / dist=&distr.;
	repeated subject=&idVar. / within=&timeVar. corr=ind corrw;
	OUTPUT OUT=GEEout XBETA=xb RESRAW = rraw;
	ods output GEEEmpPEst=betaGEE;
RUN;
ods output close;
ods select all; 


*Obtain GMM Estimates;
PROC IML;
*Do not import los_2 since there are no moment conditions; 
USE mydata3;
READ all VAR {&adjPredVar.} INTO Zmat;
READ all VAR {&outVar.} INTO yvec;
READ all VAR {&idVar.} INTO ID;
READ all VAR {&timeVar.} INTO time;
CLOSE mydata3;

N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
Np = ncol(Zmat); *Number of predictors, not including intercept;
Pn=Np+1;                      * number of covariates/parameters TO estimate, including intercept;

*Create the X matrix;
int = j(N*T,1,1);
Xmat =j(N*T,Pn,.); Xmat[,1]=int; Xmat[,2:Pn]=Zmat; *Xmat is matrix of intercept in column 1 (vector of 1) and Zmat;

*Use GEE estimates as starting values;
USE betaGEE;
READ all VAR {Estimate} INTO beta0;
CLOSE betaGEE;
beta0 = t(beta0);

*Use the identified covariate types from V40 code (multiple test);
USE TypeMtxCombo3;
READ all INTO TypeMtx;
CLOSE TypeMtxCombo3;

*Number of valid moment conditions;
neq = TypeMtx[,+];

nloc = j(1,Pn+1,0);          * nloc containing the starting/END positions of reg eqs brought by each covariate  ;
DO p =1 TO Pn;
  nloc[p+1] = sum(neq[1:p]); 
END;

nreg = sum(neq); * total number of reg_eq, 26;

Wn = I(nreg);                * I(26) initial weight matrix ;

*Module to calculate the objective function;
START TSGMM(beta) global(Pn,T,N,Xmat,yvec,nreg,TypeMtx,nloc,Wn);      
Gn = j(nreg,1,0);             * TO initialize Gn whenever objfun is called, b/c Gn depends on beta ;   
eq = j(nreg,N,0);              * 26x1625 reg_eqs for a type I covariate, for all persons  ;
S = j(nreg,nreg,0); 

DO i = 1 TO N;               * DO loops for each person;
  x = Xmat[(i-1)*T+1:i*T,];  * T x Pn, READ data for each person across all times;
  y = yvec[(i-1)*T+1:i*T];   * T outcomes for person i;

  if "&distr." = "bin" then 
  	mu = exp(x*t(beta)) / ( 1+exp(x*t(beta)) );
  else if "&distr." = "normal" then
 	mu = x*t(beta);

  Rsd = y - mu;                * T x 1;
  DO p = 1 TO Pn; * Do loops for each predictor;
  	if "&distr." = "bin" then
    	D = x[,p]#mu#(1- mu);
	else if "&distr." = "normal" then
		D= x[,p]#(mu##(-1));
    Eqmtx = Rsd*t(D);	*Gives x#mu#(1-mu)#(y-mu); *Do the time points not matter?;
    eq[nloc[p]+1:nloc[p+1],i] = Eqmtx[loc(TypeMtx[p,]^=0)];    *Save calculations to eq;
  END; *End loop across predictors;
  S = S + eq[,i]*t(eq[,i]);
END;     * TO CLOSE "DO i = 1 TO N";
Wn = ginv(S/N);		* Update the weight matrix;
Gn = eq[,:];               * row mean => col vector  neq x 1;
f = t(Gn)*Wn*Gn;           * the objective fn TO be minimized; 
RETURN(f);
FINISH TSGMM;


tc = {2000 2000}; optn = {0 1};              * optn[1]=0 specifies a minimization problem;
*Optimization takes: response code, row vector of optimal pts, module, starting value, optimization details, threshold conditions;
  *Run with the GEE estimates as starting parameters;
if "&optim."="NLPNRA" then
	CALL NLPNRA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson;
else if "&optim."="NLPNRR" then
	CALL NLPNRR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson ridge;
else if "&optim."="NLPQN" then
	CALL NLPQN(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization dual quasi-newton;
else if "&optim."="NLPCG" then
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
else if "&optim."="NLPQUA" then
	CALL NLPQUA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization quadratic;
else if "&optim."="NLPTR" then
	CALL NLPTR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization trust region;
else if "&optim."="NLPNMS" then
	CALL NLPNMS(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization nelder-mead simplex;
else do;
	print "Must use valid optimization, will run with NLPCG";
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
end;
beta = xres;


* Module to calculate the asymptotic variance;
DG = j(nreg,Pn,.);    * nreg x Pn;
DO k = 1 TO Pn; 
  DGi = j(nreg,N,0);          
  DO i = 1 TO N;               * DO loops for each for each person;
    x = Xmat[(i-1)*T+1:i*T,];  * T x Pn, READ data for each person across all times;
    y = yvec[(i-1)*T+1:i*T];  
	if "&distr." = "bin" then 
    	mu = exp(x*t(beta)) / ( 1+exp(x*t(beta)) );
	else if "&distr." = "normal" then
		mu = x*t(beta);
    Rsd = y - mu;                * T x 1;
	if "&distr." = "bin" then do;
    	Dk =  x[,k]#mu#(1- mu);
		Dkz =  x[,k]#(1- 2*mu);
		end;
	else if "&distr." = "normal" then do;
		Dk =  x[,k]#(mu##(-1));
		Dkz =  x[,k]#(-mu##(-2));
		end;
    DO p = 1 TO Pn;
	  if "&distr." = "bin" then 
      	Dp = x[,p]#mu#(1- mu);
	  else if "&distr." = "normal" then
		Dp = x[,p]#(mu##(-1));
      Dkzp = Dkz#Dp;                  * T x 1;
	  DGmtx = Rsd*t(Dkzp)-Dk*t(Dp);
      DGi[nloc[p]+1:nloc[p+1],i] = DGmtx[loc(TypeMtx[p,]^=0)];   
    END;
  END;
  DG[,k]= DGi[,:];    * row mean => col vector  neq x 1;
END;   



AsymVar = (1/N)*ginv(t(DG)*Wn*DG);   * Pn x Pn  note Wn is the inverse ;
AVvec = vecdiag(AsymVar); *Take only the diagonal elements for variance;
StdDev = sqrt(AVvec);

*Calculate the p-Values;
zvalue = t(beta)/StdDev;
pvalue = 2*(1-cdf('normal',abs(zvalue)));

*Create a matrix of results;
Outmtx = j(Pn,4,.);
ColLabel={"Estimate" "StdDev" "Zvalue" "Pvalue"};
Varnames=t({"Intercept"} || {&adjPredVar.});
OutTitle = "Analysis of Partial GMM Estimates";

Outmtx[,1]=t(beta);
Outmtx[,2]=StdDev;
Outmtx[,3]=zvalue;
Outmtx[,4]=pvalue;
PRINT Outmtx[colname=colLabel rowname=Varnames label=OutTitle ];

*Calculate residuals with GMM estimates;
resvec = yvec - Xmat*t(beta);

*Save the results as a dataset;
CREATE outdat1 FROM Outmtx;
APPEND FROM Outmtx;
CLOSE outdat1;

QUIT;

%mend partitionedGMM; *End macro code;
