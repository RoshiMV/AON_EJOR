Scalar
tcomp, texec, telapsed;
Sets
i suppliers /1*5/
m SAA scenarios /1*5000/
j   coefficients  /Demand , SP/
u scenarios /1*32/


Scalars
h cost of unsold item / 5 /
rc reservation cost at the backup supplier /20 /
ec execution cost at the backup supplier /31 /
mu demand mean / 800 /
sd demand stdev / 50 /
muSP spot market price mean / 240 /
sdSP spot market price stdev / 10 /;

Parameters
oCost(i) ordering cost
/
1	0
2	0
3	0
4	0
5	0
/
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
RPrime(u,i)  availability
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
Parameters 
data(m,j) demand and SP values in each scenario 
R(m,i) availability of the suppliers in each scenario
sumOrders, sumSelectedSuppliers; 

Variables
expected expected profit;

Positive variables
TotalQ(m) total order quantity for scenario m 
Ipos(m)   amount of overage in scenarion m 
Ineg(m)   amount of underage in scenario m 
O(m)      amount of the option contract in scenario m
S(m)      spot market purchase in scenario m
cost(m)  cost scenario m
q(i) order quantity for supplier i
KR   capacity reserved at the backup supplier;

Equations
objective        define objective function
scenarioCost(m)  scenario cost definition 
calcTotalQ(m)    calculate total order for each scenarion
overage(m)  
underage(m) 
underage2(m)  
optionContract(m); 

objective..         expected=e= (1/5000) * sum((m),cost(m));
scenarioCost(m)..   cost(m)=e=  sum((i),q(i)*(c(i)*R(m,i)+oCost(i)))+ rc*KR + h*Ipos(m)+ 
                         ec* O(m) + data(m, 'SP')*S(m); 
calcTotalQ(m)..   TotalQ(m) =e=  sum((i), q(i)*R(m,i)); 
overage(m)..      Ipos(m) =g= TotalQ(m) - data(m, 'Demand');
underage(m)..    Ineg(m) =g=  data(m, 'Demand')- TotalQ(m); 
underage2(m)..   Ineg(m) =e= O(m)+S(m);  
optionContract(m)..  O(m) =l= KR;


Model AONnew / all / ;
$include demand_5_0_set_9.inc
$include SP_5_0_set_9.inc
$include Rdata_5_0_set_9.inc


option optcr = 1e-5;
option iterlim = 100000000;
option reslim = 5000;
option domlim = 20000;
solve  AONnew using lp minimizing expected;
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
Parameter
RealObj
zReal(u) service level for scenario u
zKRReal(u)  
lossReal(u)   loss value for Qu scenario u
lossKRReal(u) loss value for Qu + KR for scenario u
pdfReal(u)    pdf for scenario u
cdfReal(u)    cdf for scenario u 
TotalQReal(u) total order quantity for scenario u
pdfKRReal(u)  pdf including KR for scenario u 
cdfKRReal(u)  cdf including KR for scenario u
costReal(u)   cost scenario u  ;

Loop(u,
TotalQReal(u)=  sum((i), q.l(i)*RPrime(u,i));
zReal(u)=(TotalQReal(u)-mu)/sd;
zKRReal(u)=(TotalQReal(u)+KR.l-mu)/sd;
pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156);
cdfReal(u) = errorf(zReal(u));
lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));
pdfKRReal(u) = exp(-0.5*power(zKRReal(u), 2)) / sqrt(2 * 3.14156);
cdfKRReal(u)= errorf(zKRReal(u));

lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));
costReal(u)=  sum((i),q.l(i)*(c(i)*RPrime(u,i)+oCost(i)))+ rc*KR.l + 
                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+
                          delta * lossKRReal(u); 
);

RealObj = sum((u), Prob(u)*costReal(u));

sumOrders = sum((i), q.l(i));    
sumSelectedSuppliers = 0; 
loop(i, 
   if(q.l(i)>0,  
      sumSelectedSuppliers= sumSelectedSuppliers+1;
   );  
);

Display AONnew.modelstat;
Display AONnew.solvestat;
tcomp = TimeComp;
texec = TimeExec;
telapsed = TimeElapsed;
*Display tcomp, texec, telapsed;

file  Results storing the results / results_5_SAA.txt / ; 
Results.ap=1;
put Results;
put /expected.l','telapsed','sumSelectedSuppliers','sumOrders','KR.l','RealObj/;
