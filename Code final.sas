options nosource;
options nonotes;

/*-------------------------------------------------------DGP Indépendance----------------------------------------------*/

proc iml;
call randseed(550);

/* Metrics Initialization */
p_fit = 0;
over = 0;
l_over = 0;
g_under = 0;
b_under = 0;
fail = 0;
sev = 0;

/* Data Generating Process */
N = 100;
mu = j(50,1,0);
Cov = I(50);
iter = 1000;
do i=1 to iter;
X = randnormal(N,mu,cov);
eps = normal(j(N,1,1));
Y = 1.5*X[,1]+1.3*X[,2]-1.7*X[,3]+1.6*X[,4]-1.4*X[,5]+1.2*X[,6]+eps*0.01;

/* Saving our data in a table*/
R = Y||X;
cnames="Y"//("X1":"X50")`;
create test from R[colname=cnames];
append from R;
close test;

/* Signal-to-noise ratio */
beta = inv(X`*X)*X`*Y;
Yhat = X*beta;
epshat = Y-Yhat;
varYhat = ((Yhat-Yhat[:])`*(Yhat-Yhat[:]))/(N-6);
varepshat = (epshat`*epshat)/(N-6);
SNR = varYhat/varepshat;

/* Glmselect procedure */
submit;
proc glmselect data=test outdesign=mod noprint;
	model Y = X1-X50
	/selection=lasso(stop=sbc choose=press lscoeffs);
run;

/* Metrics calculation by using results from the glmselect procedure */
data mod; set mod; drop y;
proc transpose data=mod out=smodel;
run;
endsubmit;
use smodel;
read all var {_name_} into obj; 
close;
ref = {"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};
resu = xsect(ref,obj);
cardi_sm = nrow(obj);
cardi_res = nrow(t(resu));
cardi_ref = nrow(ref);
diff = cardi_sm - cardi_ref;
if cardi_res = cardi_sm & cardi_ref = cardi_sm then do;
p_fit = p_fit+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res = cardi_sm then do;
g_under = g_under+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res < cardi_sm then do;
b_under = b_under+(1/iter)*100; end;
if cardi_ref < cardi_sm & cardi_res = 7 then do;
over = over+(1/iter)*100; end;
if cardi_res = 7 & diff > 4 then do;
l_over = l_over+(1/iter)*100; end;
if cardi_ref <= cardi_sm & cardi_res < cardi_ref then do;
fail = fail+(1/iter)*100; end;
if cardi_sm <= 3 then do;
sev = sev+(1/iter)*100; end;
end;

/* Creation of table in order to display our results, using the name of the "choose" criterion */
metri_name = {'Overfit','Large_Overfit','Bad_Underfit','Good_Underfit','Severity','Perfect_fit','Fail'};
total_metri = over+b_under+g_under+p_fit+fail;
Frequency = over//l_over//b_under//g_under//sev//p_fit//fail;
total = j(nrow(Frequency),1,total_metri);
create Resu_press var {"metri_name", "Frequency", "total"};
append; 
close Resu_press;
quit;

/*------------------------------------------DGP Multicolinéarité-------------------------------------------------------*/

/* Interne */

proc iml;
call randseed(550);

/* Metrics Initialization */
p_fit = 0;
over = 0;
l_over = 0;
g_under = 0;
b_under = 0;
fail = 0;
sev = 0;

/* Data Generating Process */
N = 100;
cov1 = toeplitz({1 0.9 0.8 0.7 0.6 0.5}); /*specify the type of colinearity, here only with 6 features*/
mu1 = j(6,1,0);
mu2 = j(44,1,0);
cov2 = I(44);
iter = 1000;
do i=1 to iter;
X1 = randnormal(N,mu1,cov1);
X2 = randnormal(N,mu2,cov2);
eps = normal(j(N,1,1));
Y = 1.5*X1[,1]+1.3*X1[,2]-1.7*X1[,3]+1.6*X1[,4]-1.4*X1[,5]+1.2*X1[,6]+eps*0.01;

/* Saving our data in a table */
RRR = Y||X1||X2;
cnnames = "Y"//("X1":"X6")`//("X7":"X50")`;
create testa from RRR[colname=cnnames];
append from RRR;
close testa;

/* Signal-to-noise ratio */
beta = inv(X1`*X1)*X1`*Y;
Yhat = X1*beta;
epshat = Y-Yhat;
varYhat = ((Yhat-Yhat[:])`*(Yhat-Yhat[:]))/(N-6);
varepshat = (epshat`*epshat)/(N-6);
SNR = varYhat/varepshat;

/* Glmselect procedure */
submit;
proc glmselect data=testa outdesign=mod noprint;
	model Y = X1-X50
	/selection=stepwise(stop=cv choose=press);
run;

/* Metrics calculation by using results from the glmselect procedure */
data mod; set mod; drop y;
proc transpose data=mod out=smodel;
run;
endsubmit;
use smodel;
read all var {_name_} into obj; 
close;
ref = {"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};
resu = xsect(ref,obj);
cardi_sm = nrow(obj);
cardi_res = nrow(t(resu));
cardi_ref = nrow(ref);
diff = cardi_sm - cardi_ref;
if cardi_res = cardi_sm & cardi_ref = cardi_sm then do;
p_fit = p_fit+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res = cardi_sm then do;
g_under = g_under+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res < cardi_sm then do;
b_under = b_under+(1/iter)*100; end;
if cardi_ref < cardi_sm & cardi_res = 7 then do;
over = over+(1/iter)*100; end;
if cardi_res = 7 & diff > 4 then do;
l_over = l_over+(1/iter)*100; end;
if cardi_ref <= cardi_sm & cardi_res < cardi_ref then do;
fail = fail+(1/iter)*100; end;
if cardi_sm <= 3 then do;
sev = sev+(1/iter)*100; end;
end;

/* Creation of table in order to display our results, using the name of the "choose" criterion */
metri_name = {'Overfit','Large_Overfit','Bad_Underfit','Good_Underfit','Severity','Perfect_fit','Fail'};
total_metri = over+b_under+g_under+p_fit+fail;
Frequency = over//l_over//b_under//g_under//sev//p_fit//fail;
total = j(nrow(Frequency),1,total_metri);
create Resu_press var {"metri_name", "Frequency", "total"};
append; 
close Resu_press;
quit;

/* Externe */

proc iml;
call randseed(550);

/* Metrics Initialization */
p_fit = 0;
over = 0;
l_over = 0;
g_under = 0;
b_under = 0;
fail = 0;
sev = 0;

/* Data Generating Process */
N = 100;
cov1 = toeplitz({1 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.54 0.53 0.51 0.50 0.45 }); 
mu1 = j(15,1,0);
mu2 = j(35,1,0);
cov2 = I(35);
iter = 1000;
do i=1 to iter;
X1 = randnormal(N,mu1,cov1);
X2 = randnormal(N,mu2,cov2);
eps = normal(j(N,1,1));
Y = 1.5*X1[,1]+1.3*X1[,2]-1.7*X1[,3]+1.6*X1[,4]-1.4*X1[,5]+1.2*X1[,6]+eps*0.01;

/* Saving our data in a table */
RRR = Y||X1||X2;
cnnames = "Y"//("X1":"X6")`//("X7":"X50")`;
create testa from RRR[colname=cnnames];
append from RRR;
close testa;

/* Signal-to-noise ratio */
beta = inv(X1`*X1)*X1`*Y;
Yhat = X1*beta;
epshat = Y-Yhat;
varYhat = ((Yhat-Yhat[:])`*(Yhat-Yhat[:]))/(N-6);
varepshat = (epshat`*epshat)/(N-6);
SNR = varYhat/varepshat;

/* Glmselect procedure */
submit;
proc glmselect data=testa outdesign=mod noprint;
	model Y = X1-X50
	/selection=elsticnet(stop=cv choose=cv);
run;

/* Metrics calculation by using results from the glmselect procedure */
data mod; set mod; drop y;
proc transpose data=mod out=smodel;
run;
endsubmit;
use smodel;
read all var {_name_} into obj; 
close;
ref = {"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};
resu = xsect(ref,obj);
cardi_sm = nrow(obj);
cardi_res = nrow(t(resu));
cardi_ref = nrow(ref);
diff = cardi_sm - cardi_ref;
if cardi_res = cardi_sm & cardi_ref = cardi_sm then do;
p_fit = p_fit+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res = cardi_sm then do;
g_under = g_under+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res < cardi_sm then do;
b_under = b_under+(1/iter)*100; end;
if cardi_ref < cardi_sm & cardi_res = 7 then do;
over = over+(1/iter)*100; end;
if cardi_res = 7 & diff > 4 then do;
l_over = l_over+(1/iter)*100; end;
if cardi_ref <= cardi_sm & cardi_res < cardi_ref then do;
fail = fail+(1/iter)*100; end;
if cardi_sm <= 3 then do;
sev = sev+(1/iter)*100; end;
end;

/* Creation of table in order to display our results, using the name of the "choose" criterion */
metri_name = {'Overfit','Large_Overfit','Bad_Underfit','Good_Underfit','Severity','Perfect_fit','Fail'};
total_metri = over+b_under+g_under+p_fit+fail;
Frequency = over//l_over//b_under//g_under//sev//p_fit//fail;
total = j(nrow(Frequency),1,total_metri);
create Resu_cv var {"metri_name", "Frequency", "total"};
append; 
close Resu_cv;
quit;


/*-------------------------------------------------------DGP Outlier---------------------------------------------------*/

proc iml;
call randseed(550);
/* Metrics Initialization */
p_fit = 0;
over = 0;
l_over = 0;
g_under = 0;
b_under = 0;
fail = 0;
sev = 0;

/* Data Generating Process */
/* Loop to generate outliers in our data */
N = 100;
mu = j(50,1,0);
Cov = I(50);
iter = 1000;
do m=1 to iter;
do rep=1 to 1;
	X = RandNormal(N, mu, cov);
	do i=1 to 45; /*this line specify how much of our features are subject to outliers*/
    	do j=1 to nrow(X);
			u = uniform(0);
			if u<0.9 then X[j,i]=normal(0);
			else X[j,i]=normal(0)+5;
		end;
	end;
end;
eps = normal(j(N,1,1));
Y = 1.5*X[,1]+1.3*X[,2]-1.7*X[,3]+1.6*X[,4]-1.4*X[,5]+1.2*X[,6]+eps*0.01;

/* Saving our data in a table */
R = Y||X;
cnames = "Y" // ("X1":"X50")`;
create test from R[colname=cnames];
append from R;
close test;

/* Signal-to-noise ratio */
beta = inv(X`*X)*X`*Y;
Yhat = X*beta;
epshat = Y-Yhat;
varYhat = ((Yhat-Yhat[:])`*(Yhat-Yhat[:]))/(N-6);
varepshat = (epshat`*epshat)/(N-6);
SNR = varYhat/varepshat;

/* Glmselect procedure */
submit;
proc glmselect data=test outdesign=mod noprint;
	model Y = X1-X50
	/selection=stepwise(stop=cv choose=press);
run;

/* Metrics calculation by using results from the glmselect procedure */
data mod; set mod; drop y;
proc transpose data=mod out=smodel;
run;
endsubmit;
use smodel;
read all var {_name_} into obj; 
close;
ref = {"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};
resu = xsect(ref,obj);
cardi_sm = nrow(obj);
cardi_res = nrow(t(resu));
cardi_ref = nrow(ref);
diff = cardi_sm - cardi_ref;
if cardi_res = cardi_sm & cardi_ref = cardi_sm then do;
p_fit = p_fit+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res = cardi_sm then do;
g_under = g_under+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res < cardi_sm then do;
b_under = b_under+(1/iter)*100; end;
if cardi_ref < cardi_sm & cardi_res = 7 then do;
over = over+(1/iter)*100; end;
if cardi_res = 7 & diff > 4 then do;
l_over = l_over+(1/iter)*100; end;
if cardi_ref <= cardi_sm & cardi_res < cardi_ref then do;
fail = fail+(1/iter)*100; end;
if cardi_sm <= 3 then do;
sev = sev+(1/iter)*100; end;
end;

/* Creation of table in order to display our results, using the name of the "choose" criterion */
metri_name = {'Overfit','Large_Overfit','Bad_Underfit','Good_Underfit','Severity','Perfect_fit','Fail'};
total_metri = over+b_under+g_under+p_fit+fail;
Frequency = over//l_over//b_under//g_under//sev//p_fit//fail;
total = j(nrow(Frequency),1,total_metri);
create Resu_press var {"metri_name", "Frequency", "total"};
append; 
close Resu_press;
quit;

/*------------------------------------------DGP Outlier et Multicolinéarité--------------------------------------------*/

proc iml;
call randseed(550);

p_fit = 0;
over = 0;
l_over = 0;
g_under = 0;
b_under = 0;
fail = 0;
sev = 0;

/* Data Generating Process */
/* Loop to generate outliers in our data */
N = 100;
mu = j(50,1,0);
Cov = I(50);
iter = 1000;
do m=1 to iter;
do rep=1 to 1;
X = RandNormal(N, mu, cov);
	do i=1 to 40; /*this line specify how much of our features are subject to outliers*/
    	do j=1 to nrow(X);
			u = uniform(0);
			if u<0.9 then X[j,i]=normal(0);
			else X[j,i]=normal(0)+5;
		end;
	end;
end;

/* Saving our data in a table */
R=X;
cnames=("X1":"X50")`;
create testi from R[colname=cnames];
append from R;
close testi;

use testi; 
read all var _NUM_ into K; 
close testi;

/* Creating 2 matrix, lili with independence from X7 to X50 and
   lala with colinearity between the first 6 features*/
L=K[, 7:50];
cccnames=("X7":"X50")`;
create lili from L[colname=cccnames];
append from L;
close lili;

P=K[, 1:6];
ccnames=("X1":"X6")`;
create lala from P[colname=ccnames];
append from P;
close lala;


/*Iman Connover */
load module=(Imanconovertransform);/* load the function definition */
 
/* Step 1: Read in the data */
varNames = {'X1', 'X2', 'X3', 'X4', 'X5', 'X6'};
use lala; read all var varNames into XX; close;
 
/* Step 2: specify target rank correlation */
C = toeplitz({1 0.95 0.94 0.93 0.92 0.90});
W = ImanConoverTransform(XX, C);

/* write new data to a SAS data set */
newvarNames = {'X1', 'X2', 'X3', 'X4', 'X5', 'X6'};
create SimCorr from W[colname=newvarNames];  append from W;  close;

/*Merging the lili and Simcorr in order to have one dataset with all features, 
  using matrix X44 and X55*/
use Simcorr; 
read all var _NUM_ into X44; 
close Simcorr;

use lili; 
read all var _NUM_ into X55; 
close lili;

XF = X44||X55;
eps = normal(j(N,1,1));
Y = 1.5*XF[,1]+1.3*XF[,2]-1.7*XF[,3]+1.6*XF[,4]-1.4*XF[,5]+1.2*XF[,6]+eps*0.01;

/* Saving our data in a table */
R = Y||XF;
cnames = "Y" // ("X1":"X50")`;
create testa from R[colname=cnames];
append from R;
close testa;

/* Signal-to-noise ratio*/
beta = inv(XF`*XF)*XF`*Y;
Yhat = XF*beta;
epshat = Y-Yhat;
varYhat = ((Yhat-Yhat[:])`*(Yhat-Yhat[:]))/(N-6);
varepshat = (epshat`*epshat)/(N-6);
SNR = varYhat/varepshat;

/* Glmselect Procedure */
submit;
proc glmselect data=testa outdesign=mod noprint;
	model Y = X1-X50
	/selection=stepwise(stop=aicc choose=aic);
run;

/* Metrics calculation by using results from the glmselect procedure */
data mod; set mod; drop y;
proc transpose data=mod out=smodel;
run;
endsubmit;
use smodel;
read all var {_name_} into obj; 
close;
ref = {"Intercept", "X1", "X2", "X3", "X4", "X5", "X6"};
resu = xsect(ref,obj);
cardi_sm = nrow(obj);
cardi_res = nrow(t(resu));
cardi_ref = nrow(ref);
diff = cardi_sm - cardi_ref;
if cardi_res = cardi_sm & cardi_ref = cardi_sm then do;
p_fit = p_fit+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res = cardi_sm then do;
g_under = g_under+(1/iter)*100; end;
if cardi_ref > cardi_sm & cardi_res < cardi_sm then do;
b_under = b_under+(1/iter)*100; end;
if cardi_ref < cardi_sm & cardi_res = 7 then do;
over = over+(1/iter)*100; end;
if cardi_res = 7 & diff > 4 then do;
l_over = l_over+(1/iter)*100; end;
if cardi_ref <= cardi_sm & cardi_res < cardi_ref then do;
fail = fail+(1/iter)*100; end;
if cardi_sm <= 3 then do;
sev = sev+(1/iter)*100; end;
end;

/* Creation of table in order to display our results, using the name of the "choose" criterion */
metri_name = {'Overfit','Large_Overfit','Bad_Underfit','Good_Underfit','Severity','Perfect_fit','Fail'};
total_metri = over+b_under+g_under+p_fit+fail;
Frequency = over//l_over//b_under//g_under//sev//p_fit//fail;
total = j(nrow(Frequency),1,total_metri);
create Resu_aic var {"metri_name", "Frequency", "total"};
append; 
close Resu_aic;
quit;
