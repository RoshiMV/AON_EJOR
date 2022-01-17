$FuncLibIn stolib stodclib
function icdfNormalTEMP     /stolib.icdfNormal/
*-------------------------------------------------------------------------------
Sets
i suppliers      /1*3/
iter iterator    /1*4/
u scenarios      /1*8/
cycle(iter);
alias (iter, k);
alias(iter,l);
parameter A(iter,k);
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
index position of the backup supplier according to the critical value
fractileKR, cdfINVKR, piPrime(k),crPrime(k), cdfINVPrime(k), cdfINV(iter);
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
zSP, pdfSP, cdfSp, delta, counter, check;

zSP = (ec - muSP) / sdSP;
pdfSP = exp(-0.5*power(zSp,2))/sqrt(2*3.14156);
cdfSP =  0.5*(1+sqrt(1-exp(-5*power(zSP,2)*0.125)));
If (zSP >= 0,
         delta = sdSP * pdfSP -  sdSP * zSP *(1-cdfSP);
else
         delta = sdSP * pdfSP -  sdSP * zSP *(cdfSP);
);


*----------- Find the index (position) of KR among suppliers
fractileKR= (delta-rc) / delta;
cdfINVKR =icdfNormalTEMP(fractileKR,mu,sd);
index=1;
loop(i,
if (fractileKR < cr(i),
        index=index+1;
));


*----------- Update pi, cr and cdfINV vectors by inserting KR in its position
loop(k,
loop(i,
if ((ord(k) < index) and (ord(k)  eq ord(i)),
         piPrime(k) = pi(i);
         crPrime(k)= cr(i);
         cdfINVPrime(k)= cdfINV_temp(i);
elseif (ord(k) = index),
         piPrime(k)=delta/(h+muSP);
         cdfINVPrime(k)= cdfINVKR;
         crPrime(k)= fractileKR;
         cdfINVPrime(k)= cdfINVKR;
elseif ((ord(k) > index) and (ord(k) eq ord(i)+1)),
         piPrime(k) = pi(i);
         crPrime(k)= cr(i);
         cdfINVPrime(k)= cdfINV_temp(i);
);
);
);
display fractileKR;
display cdfINVKR;
display piPrime;
display crPrime;
display cdfINVPrime;
display index;
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
lossKR    loss value for Qu + KR for scenario u
TotalQ    total order quantity for scenario u
q(i)      order quantity for supplier i
qPrime(k) q for the converted problem
KR        capacity reserved at the backup supplier;

Binary variables
T determining whether z is positive or not
TKR;
*-------------------------------------------------------------------------------

Equations
objective        define objective function
calcTotalQ    calculate total order for each scenarion
calcLoss      define loss constraint Qu scenario u
calcLossKR    define loss constraint Qu + KR scenario u
calcz         calculate z
calczKR       calculate zKR
calcpdf
calccdf
calcpdfKR
calccdfKR
ConstraintForT1
ConstraintForT2
ConstraintForTKR1
ConstraintForTKR2
PrimeToReal1
PrimeToReal2
PrimeToReal3
LinConst;


objective..              expected=e= sum((i),(c(i)+h)*q(i)*pi(i))+ rc*KR -
                         h*mu + (h+muSP-delta)* loss+ delta * lossKR;
calcTotalQ..             TotalQ =e=  sum((i), pi(i)*q(i));
calcLoss..               loss =e= sd * (pdf - z*(T*(1-cdf)+(1-T)*cdf));
calcLossKR..             lossKR =e= sd * (pdfKR -zKR*(TKR*(1-cdfKR)+(1-TKR)*cdfKR));
calcz..                  TotalQ=e=mu+sd*z;
calczKR..                TotalQ+KR=e=mu+sd*zKR;
calcpdf..                pdf*sqrt(2*3.14156)=e=exp(-0.5*power(z,2));
calccdf..                cdf=e= 0.5*(1+sqrt(1-exp(-5*power(z,2)*0.125)));
calcpdfKR..              pdfKR*sqrt(2*3.14156)=e=exp(-0.5*power(zKR,2));
calccdfKR..              cdfKR=e=  0.5*(1+sqrt(1-exp(-5*power(zKR,2)*0.125)));
ConstraintForT1..        z =l= 100 * T;
ConstraintForT2..        z =g= 100 * (T - 1);
ConstraintForTKR1..      zKR =l= 100 * TKR;
ConstraintForTKR2..      zKR =g= 100 * (TKR - 1);
PrimeToReal1(i,k)$((ord(k) < index) and (ord(i) eq ord(k)))..    q(i)=e=qPrime(k);
PrimeToReal2(k)$(ord(k) eq index)..    KR=e=qPrime(k);
PrimeToReal3(i,k)$((ord(k) > index) and (ord(i)+1 eq ord(k)))..    q(i)=e=qPrime(k);
LinConst(cycle)..        sum((k)$(ord(k) le counter),(A(cycle,k)*qPrime(k))) =e= cdfINV(cycle);

*-------------------------------------------------------------------------------

Model AONnew /All/;

option limrow=100;
cycle(iter)=no;
*-------------------------------------------------------------------------------
******** Updating the left-hand-side : matrix A

loop(iter,
cycle(iter-1)$(ord(iter)<>1)=no;
cycle(iter)=yes;

A(cycle,k)$(ord(iter) <> ord(k))=piPrime(k);
A(cycle,k)$(ord(iter) eq ord(k))=1;
display A;
);

cycle(iter)=no;
*-------------------------------------------------------------------------------
******** Updating the right-hand-side : parameter

loop(iter,
cycle(iter-1)$(ord(iter)<>1)=no;
cycle(iter)=yes;
   loop(k,
      if (ord(iter)=ord(k),
         cdfINV(cycle)=cdfINVPrime(k);
     ));
display cdfINVPrime;
);
*-------------------------------------------------------------------------------
********* Solving the problem
parameters

RealObj
zReal(u) service level for scenario u
zKRReal(u)
pdfReal(u)    pdf for scenario u
cdfReal(u)    cdf for scenario u
TotalQReal(u) total order quantity for scenario u
pdfKRReal(u)  pdf including KR for scenario u
cdfKRReal(u)  cdf including KR for scenario u
lossReal(u)   loss value for Qu scenario u
lossKRReal(u) loss value for Qu + KR for scenario u
costReal(u)   cost scenario u;


counter=0;
cycle(iter)=no;
loop(iter,
counter=counter+1;
cycle(iter)=yes;
qPrime.lo(k)$(ord(k) le ord(iter))=0;
qPrime.up(k)$(ord(k) le ord(iter))=+inf;
qPrime.fx(k)$(ord(k) > ord(iter))=0;
option minlp=Lindoglobal;
solve  AONnew using minlp minimizing expected;
display qPrime.l,cdfINV;
);
*-------------------------------------------------------------------------------
TotalQReal(u)=  sum((i), q.l(i)*R(u,i));
zReal(u)=(TotalQReal(u)-mu)/sd;
zKRReal(u)=(TotalQReal(u)+KR.l-mu)/sd;
pdfReal(u) = exp(-0.5*power(zReal(u),2))/sqrt(2*3.14156);
cdfReal(u) = 0.5*(1+sqrt(1-exp(-5*power(zReal(u),2)*0.125)));
Loop(u,
If (zReal(u) >= 0,
         lossReal(u) = sd * (pdfReal(u) - zReal(u)*(1-cdfReal(u)));
else
         lossReal(u) = sd * (pdfReal(u) - zReal(u)* cdfReal(u));
);
);

pdfKRReal(u)= exp(-0.5*power(zKRReal(u),2))/sqrt(2*3.14156);
cdfKRReal(u)= 0.5*(1+sqrt(1-exp(-5*power(zKRReal(u),2)*0.125)));
Loop(u,
If (zKRReal(u) >= 0,
         lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)*(1-cdfKRReal(u)));
else
         lossKRReal(u) = sd * (pdfKRReal(u) - zKRReal(u)* cdfKRReal(u));
);
);

costReal(u)=  sum((i),c(i)*q.l(i)*R(u,i))+ rc*KR.l +
                         h*(TotalQReal(u)-mu) + (h+muSP-delta)* lossReal(u)+
                         delta * lossKRReal(u);
RealObj= sum((u), Prob(u)*costReal(u));
Display RealObj;
*-------------------------------------------------------------------------------
display cycle;



