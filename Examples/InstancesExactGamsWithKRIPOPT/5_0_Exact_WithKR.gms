Sets
i suppliers /1*5/
u scenarios /1*32/

Parameters
c(i) selling price
/
1	41
2	41
3	43
4	46
5	59
/
pi(i) probability of each supplier's availability
/
1	0.609867
2	0.619785
3	0.7575
4	0.771584
5	0.802576
/
Prob(u) Scenario probability
/
1	0.00162211
2	0.00659426
3	0.00547946
4	0.0222753
5	0.00506701
6	0.0205986
7	0.0171162
8	0.0695816
9	0.00264419
10	0.0107493
11	0.00893203
12	0.0363108
13	0.00825968
14	0.0335776
15	0.0279011
16	0.113424
17	0.00253573
18	0.0103083
19	0.00856564
20	0.0348213
21	0.00792087
22	0.0322002
23	0.0267566
24	0.108772
25	0.00413347
26	0.0168035
27	0.0139628
28	0.0567619
29	0.0129118
30	0.0524893
31	0.0436157
32	0.177308
/
Table
R(u,i)  availability
	1	2	3	4	5	
1	0	0	0	0	0	
2	0	0	0	0	1	
3	0	0	0	1	0	
4	0	0	0	1	1	
5	0	0	1	0	0	
6	0	0	1	0	1	
7	0	0	1	1	0	
8	0	0	1	1	1	
9	0	1	0	0	0	
10	0	1	0	0	1	
11	0	1	0	1	0	
12	0	1	0	1	1	
13	0	1	1	0	0	
14	0	1	1	0	1	
15	0	1	1	1	0	
16	0	1	1	1	1	
17	1	0	0	0	0	
18	1	0	0	0	1	
19	1	0	0	1	0	
20	1	0	0	1	1	
21	1	0	1	0	0	
22	1	0	1	0	1	
23	1	0	1	1	0	
24	1	0	1	1	1	
25	1	1	0	0	0	
26	1	1	0	0	1	
27	1	1	0	1	0	
28	1	1	0	1	1	
29	1	1	1	0	0	
30	1	1	1	0	1	
31	1	1	1	1	0	
32	1	1	1	1	1	

Scalars
h cost of unsold item / 5 /
rc reservation cost at the backup supplier /20 /
ec execution cost at the backup supplier /31 /
mu demand mean / 800 /
sd demand stdev / 50 /
muSP spot market price mean / 240 /
sdSP spot market price stdev / 10 /
telapsed;

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
calczKR(u)       calculate zKR(u);

objective..              expected=e= sum((u), Prob(u)*cost(u));
scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ rc*KR +
                         h*(TotalQ(u)-mu) + (h+muSP-delta)* loss(u)+ 
                         delta * lossKR(u); 
calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));
calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));
calcLossKR(u)..          lossKR(u) =e= sd*(pdfKR(u) -zKR(u)*(1-cdfKR(u)));
calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));
calccdf(u)..             cdf(u)=e= errorf(z(u));
calcpdfKR(u)..           pdfKR(u)*sqrt(2*3.14156)=e=exp(-0.5*power(zKR(u),2));
calccdfKR(u)..           cdfKR(u)=e= errorf(zKR(u));
calcz(u)..               TotalQ(u)=e=mu+sd*z(u); 
calczKR(u)..             TotalQ(u)+KR=e=mu+sd*zKR(u);
cdf.up(u)=1;
cdf.lo(u)=0;
cdfKR.up(u)=1; 
cdfKR.lo(u)=0;
z.up(u)=20; 
z.lo(u)=-20;
zKR.up(u)=20;
zKR.lo(u)=-20;


Model AONnew / all / ;
option optcr = 1e-5;
option nlp=IPOPT;
option iterlim = 100000000;
option reslim = 5000;
option domlim = 2000000;
solve  AONnew using nlp minimizing expected;



Display AONnew.modelstat;
Display AONnew.solvestat;
telapsed = TimeElapsed;

file  Results storing the results / results_5_Exact_IPOPT.txt / ; 
Results.ap=1;
put Results;
*put /'Objective,Elapsed time','AONnew.solvestat','AONnew.modelstat/;
put /expected.l','telapsed','AONnew.solvestat','AONnew.modelstat/;
putclose Results;
