$FuncLibIn stolib stodclib
function icdfNormalTEMP     /stolib.icdfNormal/
Scalar
tcomp, texec, telapsed;
*-------------------------------------------------------------------------------
Sets
i suppliers      /1*3/
u scenarios      /1*8/

alias(i, j);
alias(i,iter,k);
set cycle(iter);
*-------------------------------------------------------------------------------
Parameters
c(i) selling price
/
1       47
2       52
3       58
/
pi(i) probability of each supplier's availability
/
1       0.504852
2       0.516907
3       0.853862
/
*cr(i) critical ratio
*/
*1        0.787755102
*2        0.767346939
*3        0.742857143
*/
cdfINV_temp(i) inverse of F for critical values
/
1        1199.664048
2        1182.534383
3        1163.044747
/

cdfINV(iter)
partX(i) X part of the equations
cdfINV_partX(i) inverse of F for fractions including X
fractileKR, cdfINVKR, A(iter,i);

*-------------------------------------------------------------------------------
Table
R(u,i)  availability
        1       2       3
1       0       0       0
2       0       0       1
3       0       1       0
4       0       1       1
5       1       0       0
6       1       0       1
7       1       1       0
8       1       1       1

parameter
Prob(u) Scenario probability
/
1       0.0349565
2       0.204246
3       0.0374033
4       0.218542
5       0.0356416
6       0.208249
7       0.0381364
8       0.222825
/;

Scalars
h cost of unsold item / 5 /
rc reservation cost at the backup supplier /20 /
ec execution cost at the backup supplier /71 /
mu demand mean / 1000 /
sd demand stdev / 250 /
muSP spot market price mean / 240 /
sdSP spot market price stdev / 10 /;

Parameter
zSP, pdfSP, cdfSp, delta, counter, check, gapAllowed, iterationCountLimit;

gapAllowed= 0.0001;
iterationCountLimit=10;

zSP = (ec - muSP) / sdSP;
pdfSP = exp(-0.5*power(zSp,2))/sqrt(2*3.14156);
cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));
If (zSP >= 0,
         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);
else
         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);
);
display delta;

*-------------------------------------------------------------------------------
Variables
z service level
expected expected cost;

Positive variables
pdf
cdf
loss      loss value for Qu scenario u
TotalQ    total order quantity for scenario u
q(i)      order quantity for supplier i;
*-------------------------------------------------------------------------------
Equations
objective        define objective function
calcTotalQ    calculate total order for each scenarion
calcLoss
calcz         calculate z
calcpdf
calccdf
LinConst;


objective..              expected=e= sum((i),c(i)*q(i)*pi(i))+h*(TotalQ-mu) + (h+muSP)* loss;
calcTotalQ..             TotalQ =e=  sum((i), pi(i)*q(i));
calcLoss..               loss =e=sd*(pdf -z*(1-cdf));
calcz..                  TotalQ=e=mu+sd*z;
calcpdf..                pdf*sqrt(2*3.14156)=e=exp(-0.5*power(z,2));
calccdf..                cdf=e= errorf(z);
LinConst(cycle)..        sum((i)$(ord(i) le counter),(A(cycle,i)*q(i))) =e= cdfINV(cycle);
*-------------------------------------------------------------------------------

Model AONnew /All/;

option limrow=100;
cycle(iter)=no;
*-------------------------------------------------------------------------------
******** Updating the left-hand-side : matrix A

loop(iter,
cycle(iter-1)$(ord(iter)<>1)=no;
cycle(iter)=yes;

A(cycle, i)$(ord(iter) <> ord(i)) = pi(i);
A(cycle, i)$(ord(iter) eq ord(i)) = 1;
Display A;
);
cycle(iter) = no;
*-------------------------------------------------------------------------------
******** Updating the right-hand-side : parameter

loop(iter,
cycle(iter - 1)$(ord(iter)<>1) = no;
cycle(iter) = yes;
         loop(i,
                 if (ord(iter) = ord(i),
                    cdfINV(cycle) = cdfINV_temp(i)
                 ));
display cdfINV;
);


*-------------------------------------------------------------------------------
********* Solving the problem

parameters

RealObj
RealObj_iter
RealObj_new
KR
gapPercent
iterationCount
zReal(u) service level for scenario u
zKRReal(u)
pdfReal(u)    pdf for scenario u
cdfReal(u)    cdf for scenario u
TotalQReal(u) total order quantity for scenario u
pdfKRReal(u)  pdf including KR for scenario u
cdfKRReal(u)  cdf including KR for scenario u
lossReal(u)   loss value for Qu scenario u
lossKRReal(u) loss value for Qu + KR for scenario u
costReal(u)   cost scenario u
temp
temp1
temp2
temp3
lowerB        for generating intitial X_i values randomly
upperB        for generating intitial X_i values randomly;


*---------- Start the iteration loop
RealObj=100000000;
counter=0;
cycle(iter)=no;

loop(iter,

RealObj_iter=100000000;
iterationCount=0;
counter=counter+1;
cycle(iter)=yes;
gapPercent=1;

*----- Find initial values of cdf inverse for the fractional values including X

fractileKR= (delta-rc) / delta;
cdfINVKR =icdfNormalTEMP(fractileKR,mu,sd);
upperB = fractileKR;

loop(i,
         lowerB= fractileKR * pi(i);
         temp1= fractileKR;
         temp2= pi(i)*(muSP-c(i))/delta;
         if (temp1<temp2,
                 upperB=temp1;
         else
                 upperB=temp2;
         );
*         display temp1,temp2, lowerB, upperB;
         partX(i)= Uniform(lowerB,upperB);
         temp = (pi(i) * (muSP - c(i)) - delta * partX(i))/(pi(i)*(h + muSP - delta));
         cdfINV_partX(i) =  icdfNormalTEMP(temp,mu,sd);
         display partX, cdfINV_partX;
);

*-------------------------------------------------------------------------------
while ((gapPercent>  gapAllowed) or (iterationCount<iterationCountLimit),
*while( aa<2,

loop(k,
         if (ord(k) le counter,
            cdfINV(cycle(k)) = cdfINV_partX(k)
         );
display cdfINV;
display counter;
);

*loop(k,
*         cycle(k - 1)$(ord(k)<>1) = no;
*         cycle(k) = Yes;
*         cdfINV(cycle) = cdfINV_partX(k);
*         display cdfINV;
*         display cycle;
*);
*        aa=aa+1 ;

*** Set the values of q_i by solving the Linear equations

q.lo(i)$(ord(i) le ord(iter))=0;
q.up(i)$(ord(i) le ord(iter))=+inf;
q.fx(i)$(ord(i) > ord(iter))=0;
option nlp=Conopt;
solve  AONnew using nlp minimizing expected;
display q.l;

KR = cdfINVKR - sum((i), q.l(i)*pi(i));

display KR,iterationCount;

*** Calc real objetive value for this solution
** Step 3 of the algorithm

TotalQReal(u)=  sum((i), q.l(i)*R(u,i));
zReal(u)=(TotalQReal(u)-mu)/sd;
zKRReal(u)=(TotalQReal(u)+KR-mu)/sd;
pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156);
cdfReal(u) = errorf(zReal(u));
Loop(u,
         lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));
);
pdfKRReal(u)= exp(-0.5*power(zKRReal(u),2))/sqrt(2*3.14156);
cdfKRReal(u)= errorf(zKRReal(u));

Loop(u,
         lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));
);
costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ rc*KR +
                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+
                         delta * lossKRReal(u);

RealObj_new= sum((u), Prob(u)*costReal(u));
Display RealObj_new;

If (RealObj_new < RealObj_iter,
         gapPercent = (RealObj_iter-RealObj_new)/RealObj_new;
         RealObj_iter=RealObj_new;
else
gapPercent=0;
);

display gapPercent;
iterationCount=iterationCount+1;

loop(i,

temp1= sum((u)$(R(u,i) eq 1), Prob(u)* cdfKRReal(u));
temp2= pi(i)*(muSP-c(i))/delta;
temp3= fractileKR;
         if (temp2<temp3,
                 upperB=temp2;
         else
                 upperB=temp3;
         );
lowerB=pi(i)*(delta-rc)/delta;

if ((temp1<upperB)and(temp1>lowerB),
         partX(i)=temp1;
elseif (temp1 ge upperB),
         partX(i)=upperB;
else
         partX(i)=lowerB;
);
*display temp1,temp2, temp3, partX;
temp = (pi(i) * (muSP - c(i)) - delta * partX(i))/(pi(i)*(h + muSP - delta));

cdfINV_partX(i) =  icdfNormalTEMP(temp,mu,sd);
);

);

If (RealObj_iter < RealObj,
        RealObj=RealObj_iter;
);
);


********* SOLVE with KR=0
loop(iter,

RealObj_iter=100000000;
counter=counter+1;
cycle(iter)=yes;
fractileKR= (delta-rc) / delta;
cdfINVKR =icdfNormalTEMP(fractileKR,mu,sd);


*** Set the values of q_i by solving the Linear equations

q.lo(i)$(ord(i) le ord(iter))=0;
q.up(i)$(ord(i) le ord(iter))=+inf;
q.fx(i)$(ord(i) > ord(iter))=0;
option nlp=Conopt;
solve  AONnew using nlp minimizing expected;
display q.l;

KR =0;
*** Calc real objetive value for this solution
** Step 3 of the algorithm

TotalQReal(u)=  sum((i), q.l(i)*R(u,i));
zReal(u)=(TotalQReal(u)-mu)/sd;
zKRReal(u)=(TotalQReal(u)+KR-mu)/sd;
pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156);
cdfReal(u) = errorf(zReal(u));
Loop(u,
         lossReal(u) = sd * (pdfReal(u) - zReal(u)* (1-cdfReal(u)));
);
pdfKRReal(u)= exp(-0.5*power(zKRReal(u),2))/sqrt(2*3.14156);
cdfKRReal(u)= errorf(zKRReal(u));

Loop(u,
         lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));
);
costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ rc*KR +
                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+
                         delta * lossKRReal(u);

RealObj_new= sum((u), Prob(u)*costReal(u));
Display RealObj_new;

If (RealObj_new < RealObj_iter,
         RealObj_iter=RealObj_new;
);

If (RealObj_iter < RealObj,
        RealObj=RealObj_iter;
);
);


*-------------------------------------------------------------------------------
Display RealObj;
tcomp = TimeComp;
texec = TimeExec;
telapsed = TimeElapsed;

file  Results storing the results / results_3_0_Heuristic.txt / ;
Results.ap=1;
put Results;
put /Realobj','telapsed/;
putclose Results;
*-----------------------------------------------------------------------------
display cycle;



