/*****************************************************************/
/**	SAS Code: Medicare Data Analysis						    **/
/** Programmer: Kyle Irimata									**/
/** Description: Example code for using the %partitionedGMM		**/
/**    macro as applied to the Medicare readmission study		**/
/** NOTE: The %partitionedGMM macro must be run to load			**/
/**    the macro into the current sessions BEFORE running the	**/
/**    code given below											**/
/*****************************************************************/

*Read in the data;
*Replace 'C:\' with the local folder in which MedicareData is saved;
libname ds 'C:\';
data Medicare;
	set ds.MedicareData;
run;

*The outcome is labeled as biRadmit;
*The time-dependent covariates are labeled as NDX, NPR, LOS and DX101;

*Run the partitioned GMM model by calling the macro;
*This code will utilize the Lalonde, Wilson and Yin moment check;
%partitionedGMM(
	file=Medicare,
	timeVar=time, 
	outVar=biRadmit,
	predVarTD=NDX NPR LOS DX101,
	idVar=PNUM_R,
	alpha=0.05,
	optim=NLPCG,
	mc=LWY)
;

*Run the partitioned GMM model by calling the macro;
*This code will utilize the Lai and Small moment check;
%partitionedGMM(
	file=Medicare,
	timeVar=time, 
	outVar=biRadmit,
	predVarTD=NDX NPR LOS DX101,
	idVar=PNUM_R,
	alpha=0.05,
	optim=NLPCG,
	mc=LS)
;
