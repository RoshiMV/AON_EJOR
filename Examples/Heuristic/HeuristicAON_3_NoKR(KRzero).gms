$FuncLibIn stolib stodclib
function icdfNormalTEMP     /stolib.icdfNormal/
*-------------------------------------------------------------------------------
Sets
i suppliers      /1*3/
iter iterator    /1*3/
u scenarios      /1*8/
cycle(iter);
alias(i,j);
parameter A(iter,i);
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
cr(i) critical ratio
/
1        0.787755102
2        0.767346939
3        0.742857143
/
cdfINV_temp(i) inverse of F for critical values
/
1        1199.664048
2        1182.534383
3        1163.044747
/
cdfINV(iter);
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
1        0.03495658
2        0.204245953
3        0.037403359
4        0.218542109
5        0.035641665
6        0.208248802
7        0.038136396
8        0.222825136
/;

Scalars

h cost of unsold item / 5 /
mu demand mean / 1000 /
sd demand stdev / 250 /
muSP spot market price mean / 240 /
sdSP spot market price stdev / 10 /;

Parameter
counter, check;

*-------------------------------------------------------------------------------
Variables
z service level
zKR
expected expected cost;

Positive variables
pdf
cdf
pdfKR
cdfKR
loss      loss value for Qu scenario u
TotalQ    total order quantity for scenario u
q(i)      order quantity for supplier i
KR;

Binary variables
T determining whether z is positive or not;
*-------------------------------------------------------------------------------

Equations
objective        define objective function
calcTotalQ    calculate total order for each scenarion
calcLoss      define loss constraint Qu scenario u
calcz         calculate z
calcpdf
calccdf
calcpdfKR
calccdfKR
ConstraintForT1
ConstraintForT2
LinConst
KRvalue;



objective..              expected=e= sum((i),c(i)*q(i)*pi(i))+h*(TotalQ-mu) + (h+muSP)* loss;
calcTotalQ..             TotalQ =e=  sum((i), pi(i)*q(i));
calcLoss..               loss =e= sd * (pdf - z*(T*(1-cdf)+(1-T)*cdf));
calcz..                  TotalQ=e=mu+sd*z;
calcpdf..                pdf*sqrt(2*3.14156)=e=exp(-0.5*power(z,2));
calccdf..                cdf=e= 0.5*(1+sqrt(1-exp(-5*power(z,2)*0.125)));
calcpdfKR..              pdfKR*sqrt(2*3.14156)=e=exp(-0.5*power(zKR,2));
calccdfKR..              cdfKR=e= 0.5*(1+sqrt(1-exp(-5*power(zKR,2)*0.125)));
ConstraintForT1..        z =l= 100 * T;
ConstraintForT2..        z =g= 100 * (T - 1);
LinConst(cycle)..        sum((i)$(ord(i) le counter),(A(cycle,i)*q(i))) =e= cdfINV(cycle);
*LinConst2..              sum((i),q(i)*pi(i))+KR =g= cdfINVKR;
*LinConst3(i)$(ord(i) > counter)..       sum((j)$(ord(j) <>  ord(i)),pi(j)*q(j))+KR =g= cdfINV_temp(i);
KRvalue..                KR =e= 0;

*-------------------------------------------------------------------------------

Model AONnew /All/;

option limrow=100;
cycle(iter)=no;
*-------------------------------------------------------------------------------
******** Updating the left-hand-side : matrix A

loop(iter,
cycle(iter-1)$(ord(iter)<>1)=no;
cycle(iter)=yes;

A(cycle,i)$(ord(iter) <> ord(i))=pi(i);
A(cycle,i)$(ord(iter) eq ord(i))=1;
display A;
);

cycle(iter)=no;
*-------------------------------------------------------------------------------
******** Updating the right-hand-side : parameter

loop(iter,
cycle(iter-1)$(ord(iter)<>1)=no;
cycle(iter)=yes;
   loop(i,
      if (ord(iter)=ord(i),
         cdfINV(cycle)=cdfINV_temp(i);
     ));
display cdfINV;
);

*loop(iter,
*cycle(iter-1)$(ord(iter)<>1)=no;
*cycle(iter)=yes;
*   loop(i,
*      if (ord(iter)=ord(i),
*         cdfINV(cycle)=cdfINV_temp(i);
*     ));
*display cdfINV;
*);
*-------------------------------------------------------------------------------
********* Solving the problem
parameters

RealObj
RealObj_temp
zReal(u) service level for scenario u
pdfReal(u)    pdf for scenario u
cdfReal(u)    cdf for scenario u
TotalQReal(u) total order quantity for scenario u
lossReal(u)   loss value for Qu scenario u
costReal(u)   cost scenario u;


RealObj=100000000;
counter=0;
cycle(iter)=no;
loop(iter,
counter=counter+1;
cycle(iter)=yes;
q.lo(i)$(ord(i) le ord(iter))=0;
q.up(i)$(ord(i) le ord(iter))=+inf;
q.fx(i)$(ord(i) > ord(iter))=0;
option minlp=lindoglobal;
option optcr=0;
solve  AONnew using minlp minimizing expected;
display q.l;
TotalQReal(u)=  sum((i), q.l(i)*R(u,i));
zReal(u)=(TotalQReal(u)-mu)/sd;
pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156);
cdfReal(u) = 0.5*(1+sqrt(1-exp(-5*power(zReal(u),2)*0.125)));
Loop(u,
If (zReal(u) >= 0,
         lossReal(u) = sd * (pdfReal(u) - zReal(u)*(1-cdfReal(u)));
else
         lossReal(u) = sd * (pdfReal(u) - zReal(u)* cdfReal(u));
);
);

costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ h*(TotalQReal(u)-mu) + (h+muSP)* lossReal(u);
RealObj_temp= sum((u), Prob(u)*costReal(u));

Display RealObj_temp;

If (RealObj_temp < RealObj,
         RealObj=RealObj_temp;
);

);
*-------------------------------------------------------------------------------
Display RealObj;
*-------------------------------------------------------------------------------
display cycle;

display pdf.l;

