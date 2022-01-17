Scalar
tcomp, texec, telapsed;
Sets
i suppliers /1*4/
u scenarios /1*16/

Parameters
c(i) selling price
/
1	48
2	50
3	49
4	58
/
pi(i) probability of each supplier's availability
/
1	0.687155
2	0.719459
3	0.913205
4	0.56801
/
Prob(u) Scenario probability
/
1	0.00329073
2	0.00432689
3	0.0346233
4	0.0455251
5	0.0084392
6	0.0110965
7	0.0887926
8	0.116751
9	0.00722799
10	0.00950388
11	0.0760489
12	0.0999946
13	0.0185364
14	0.024373
15	0.19503
16	0.25644
/
Table
R(u,i)  availability
	1	2	3	4	
1	0	0	0	0	
2	0	0	0	1	
3	0	0	1	0	
4	0	0	1	1	
5	0	1	0	0	
6	0	1	0	1	
7	0	1	1	0	
8	0	1	1	1	
9	1	0	0	0	
10	1	0	0	1	
11	1	0	1	0	
12	1	0	1	1	
13	1	1	0	0	
14	1	1	0	1	
15	1	1	1	0	
16	1	1	1	1	

Scalars
h cost of unsold item / 5 /
rc reservation cost at the backup supplier /20 /
ec execution cost at the backup supplier /71 /
mu demand mean / 1000 /
sd demand stdev / 250 /
muSP spot market price mean / 240 /
sdSP spot market price stdev / 10 /;

Parameter
zSP, pdfSP, cdfSp, delta;
zSP = (ec - muSP) / sdSP;
pdfSP =  exp(-0.5*power(zSp,2))/sqrt(2*3.14156);
cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));
If (zSP >= 0, 
         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);
else  
         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);
);

Variables
z(u) service level for scenario u
zKR(u)
expected expected profit;

Positive variables
pdf(u)    pdf for scenario u
cdf(u)    cdf for scenario u
TotalQ(u) total order quantity for scenario u
pdfKR(u)  pdf including KR for scenario u
cdfKR(u)  cdf including KR for scenario u
loss(u)   loss value for Qu scenario u
lossKR(u) loss value for Qu + KR for scenario u
cost(u)   cost scenario u
q(i)      order quantity for supplier i
KR        capacity reserved at the backup supplier;

Binary variables
T(u) determining whether z is positive or not
TKR(u);

Equations
objective        define objective function
scenarioCost(u)  scenario cost definition
calcTotalQ(u)    calculate total order for each scenarion
calcLoss(u)      define loss constraint Qu scenario u
calcLossKR(u)    define loss constraint Qu + KR scenario u
calcpdf(u)       calculate pdf of Qu for scenario u
calccdf(u)       calculate cdf of Qu for scenario u
calcpdfKR(u)     calculate pdf of Qu + KR
calccdfKR(u)     calculate cdf of Qu+KR
calcz(u)         calculate z(u)
calczKR(u)       calculate zKR(u)
ConstraintForT1(u)
ConstraintForT2(u)
ConstraintForTKR1(u)
ConstraintForTKR2(u);

objective..              expected=e= sum((u), Prob(u)*cost(u));
scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ rc*KR +
                         h*(TotalQ(u)-mu) + (h+muSP-delta)* loss(u)+ 
                         delta * lossKR(u); 
calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));
calcLoss(u)..            loss(u)=e= sd * (pdf(u) - z(u)*(T(u)*(1-cdf(u))+
                         (1-T(u))*cdf(u))); 
calcLossKR(u)..          lossKR(u) =e= sd * (pdfKR(u) - 
                         zKR(u)*(TKR(u)*(1-cdfKR(u))+ 
                         (1-TKR(u))*cdfKR(u)));
calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));
calccdf(u)..             cdf(u)=e= 0.5*(1+sqrt(1-exp(-5*power(z(u),2)*0.125)));
calcpdfKR(u)..           pdfKR(u)*sqrt(2*3.14156)=e=exp(-0.5*power(zKR(u),2));
calccdfKR(u)..           cdfKR(u)=e= 0.5*(1+sqrt(1-exp(-5*power(zKR(u),2)*0.125)));
calcz(u)..               TotalQ(u)=e=mu+sd*z(u); 
calczKR(u)..             TotalQ(u)+KR=e=mu+sd*zKR(u);
ConstraintForT1(u)..     z(u) =l= 100 * T(u);
ConstraintForT2(u)..     z(u) =g= 100 * (T(u) - 1);
ConstraintForTKR1(u)..   zKR(u) =l= 100 * TKR(u);
ConstraintForTKR2(u)..   zKR(u) =g= 100 * (TKR(u) - 1);


Model AONnew / all / ;
option optcr = 1e-5;
option minlp=baron;
option iterlim = 100000000;
option reslim = 5000;
option domlim = 20000;
solve  AONnew using minlp minimizing expected;


*display
*z.l, zKR.l, q.l, cost.l, KR.l, TotalQ.l, z.l, pdf.l, cdf.l; 


Display AONnew.modelstat;
Display AONnew.solvestat;
tcomp = TimeComp;
texec = TimeExec;
telapsed = TimeElapsed;
*Display tcomp, texec, telapsed;

file  Results storing the results / results_4_n_Instance_0_exact.txt / ; 
Results.ap=1;
put Results;
*put /'Objective,Elapsed time'/;
put /expected.l','telapsed/;
putclose Results;
