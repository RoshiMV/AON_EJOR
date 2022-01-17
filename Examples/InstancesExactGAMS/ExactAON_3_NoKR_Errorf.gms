Sets
i suppliers /1*3/
u scenarios /1*8/

Parameters
c(i) selling price
/
1       58
2       52
3       47
/
pi(i) probability of each supplier's availability
/
1       0.853862
2       0.516907
3       0.504852
/
Prob(u) Scenario probability
/
1       0.0349565
2       0.0356416
3       0.0374033
4       0.0381364
5       0.204246
6       0.208249
7       0.218542
8       0.222825
/
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

Scalars
h cost of unsold item / 5 /
mu demand mean / 1000 /
sd demand stdev / 250 /
muSP spot market price mean / 240 /
sdSP spot market price stdev / 10 /
telapsed;

Variables
z(u) service level for scenario u
expected expected profit;

Positive variables
pdf(u)    pdf for scenario u
cdf(u)    cdf for scenario u
TotalQ(u) total order quantity for scenario u
loss(u)   loss value for Qu scenario u
cost(u)   cost scenario u
q(i)      order quantity for supplier i;


Equations
objective        define objective function
scenarioCost(u)  scenario cost definition
calcTotalQ(u)    calculate total order for each scenarion
calcLoss(u)      define loss constraint Qu scenario u
calcpdf(u)       calculate pdf of Qu for scenario u
calccdf(u)       calculate cdf of Qu for scenario u
calcz(u)         calculate z(u);


objective..              expected=e= sum((u), Prob(u)*cost(u));
scenarioCost(u)..        cost(u)=e=  sum((i),c(i)*q(i)*R(u,i))+ h*(TotalQ(u)-mu) + (h+muSP)* loss(u);
calcTotalQ(u)..          TotalQ(u) =e=  sum((i), q(i)*R(u,i));
calcLoss(u)..            loss(u)=e=sd*(pdf(u) -z(u)*(1-cdf(u)));
calcpdf(u)..             pdf(u)*sqrt(2*3.14156)=e=exp(-0.5*power(z(u),2));
calccdf(u)..             cdf(u)=e= errorf(z(u));
calcz(u)..               TotalQ(u)=e=mu+sd*z(u);

cdf.up(u)=1;
cdf.lo(u)=0;

Model AONnew / all / ;
option optcr = 1e-5;
option nlp=Conopt;
option iterlim = 100000000;
option reslim = 500000;
option domlim = 2000000;
solve  AONnew using nlp minimizing expected;


Display AONnew.modelstat;
Display AONnew.solvestat;
telapsed = TimeElapsed;
display q.l;

file  Results storing the results / results_3_n_Instance_0_exact.txt / ;
Results.ap=1;
put Results;
put /expected.l','telapsed/;
putclose Results;
display pdf.l;
