/*Stat des */

PROC UNIVARIATE DATA=diabetes;
	VAR Y AGE SEX BMI BP S1 S2 S3 S4 S5 S6;
	HISTOGRAM /NORMAL;
RUN;

proc means data=diabetes skewness kurtosis;
run;


PROC CORR DATA=diabetes;
RUN;

proc sgplot data=diabetes;
vbox S5;
run;

proc reg data=diabetes;
model y=AGE SEX BMI BP S1 S2 S3 S4 S5 S6;
run;

/*diabete Stepwise*/

proc glmselect data=diabetes ;
	model Y = AGE SEX BMI BP S1 S2 S3 S4 S5 S6
	/selection=Stepwise(stop=aicc choose=aic select=sl);
run;
