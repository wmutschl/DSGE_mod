// RESUSCITATING REAL BUSINESS CYCLES
// R.G KING AND S.T REBELO (1999)
// THIS MOD FILE REPLICATES THE THEORITICAL MOMENTS AND IRFS OF THE
// TABLE 3, FIGURE 9 AND FIGURE 10
// THIS CODE IS WRITTEN BY STEPHANE LHUISSIER (INTERN AT CEPREMAP) - NOVEMBER 2009


//********************** SPECIFICATION OF THE MODEL **********************//

// CHOICE BETWEEN MOMENTS AND IRFS : 0 for moments
//                                   1 for irfs
@#define irf  = 1

// IF IRFS, CHOICE BETWEEN A TRANSITORY SHOCK (FIGURE 9) AND A MORE 
// PERSISTENT SHOCK (FIGURE 10) : 1 for transitory shock
//                                0 for a more permanent shock
@#define transitory  = 0

//************************END OF THE SPECIFICATION ***********************//

var A C K N L W RW Y I R LAMBDA Output Consumption Investment Capital Labor
    Leisure Productivity Real_wage Real_interest_rate;

varexo EPS;

parameters bet theta eta gam alp del rho r Kss Css Nss Lss Ass Wss
           LAMBDAss Yss Iss EPSss RWss Rss;

bet   = 0.984;
eta   = 1;
gam   = 1.004;
alp   = 0.667;
del   = 0.025;
r     = gam/bet-1;
theta = ((1-0.2)/0.2)/((1-(gam-1+del)*((1-alp)/(gam/bet-1+del)))/alp);

 @#if (transitory==0)

     rho  = 0.979;

 @#else 

     rho  = 0;

  @#endif


//Compute the Steady State
Ass          = 1;
Nss          = 1/(theta/alp*(1-(gam-1+del)*((1-alp)/(gam/bet-1+del)))+1);
Lss          = 1-Nss;
Kss          = (((1-alp)*Ass)/(r+del))^(1/alp)*Nss;
Iss          = (gam-1+del)*Kss;
Wss          = theta/Lss;
Yss          = Ass*Kss^(1-alp)*Nss^(alp);
Css          = Yss - Iss;
LAMBDAss     = 1/Css;
Rss          = r;
RWss         = Wss/LAMBDAss;
EPSss        = 0;


model;

   // 1. FIRST ORDER CONDITION : C
   1/C = LAMBDA;

   // 2. FIRST ORDER CONDITION : L
   theta*L^(-eta)=W;

   // 3. FIRST ORDER CONDITION : K(+1)
   LAMBDA*A*K(-1)^(1-alp)*alp*N^(alp-1) = W;

   // 4. MARGINAL PRODUCTIVITY OF LABOR
   bet*LAMBDA(+1)*(A(+1)*(1-alp)*K^(-alp)*N(+1)^(alp)+1-del) = gam*LAMBDA;

   // 5. MARGINAL PRODUCTIVITY OF CAPITAL
   N + L = 1;

   // 6. OUTPUT
   Y = A*K(-1)^(1-alp)*N^(alp);

   // 7.  BUDGET CONSTRAINT
   C + I = Y;

   // 8. LAW OF MOTION FOR THE CAPITAL STOCK
   gam*K = (1-del)*K(-1) + I;

   // 9. TECHNOLOGY SHOCK PROCESSUS
   log(A) = rho*log(A(-1)) + EPS;

   // 10. REAL WAGE RATE
   RW = W/LAMBDA;

   // 11. REAL INTEREST RATE
   R = (1-alp)*A(+1)*K^(-alp)*N(+1)^(alp)-del;


// Measurement equations

@#if (irf==0)

Output             = log(Y/Yss)*100;
Consumption        = log(C/Css)*100;
Investment         = log(I/Iss)*100;
Capital            = log(K/Kss)*100;
Labor              = log(N/Nss)*100;
Productivity       = log(A/Ass)*100;
Leisure            = log(L/Lss)*100;
Real_wage          = log(RW/RWss)*100;
Real_interest_rate = R*100;

@#else

Output             = log(Y/Yss);
Consumption        = log(C/Css);
Investment         = log(I/Iss);
Capital            = log(K/Kss);
Labor              = log(N/Nss);
Productivity       = log(A/Ass);
Leisure            = log(L/Lss);
Real_wage          = log(RW/RWss);
Real_interest_rate = R*400;

@#endif

end;

initval;

C                  = Css;
K                  = Kss;
N                  = Nss;
W                  = Wss;
Y                  = Yss;
I                  = Iss;
A                  = Ass;
LAMBDA             = LAMBDAss;
L                  = Lss;
RW                 = RWss;
R                  = Rss;
EPS                = EPSss;
Output             = 0;
Consumption        = 0;
Investment         = 0;
Capital            = 0;
Labor              = 0;
Productivity       = 0;
Leisure            = 0;
Real_wage          = 0;

@#if (irf==0)
   Real_interest_rate = Rss*100;
   
@#else
   Real_interest_rate = Rss*400;

@#endif

end;

steady;
check;

shocks;
var EPS;

@#if (irf==1)

    stderr 1;

@#else

    stderr 0.0072;

@#endif

end;

@#if (irf==1)

     @#if (transitory==1)

stoch_simul(nomoments,order=1,irf=20)
Productivity Capital Consumption Investment Output Labor Real_wage
Real_interest_rate;

     @#else

stoch_simul(nomoments,order=1,irf=90)
Productivity Capital Consumption Investment Output Labor Real_wage
Real_interest_rate;

     @#endif
@#else

stoch_simul(hp_filter=1600,order=1,irf=0)
Output Consumption Investment  Labor Real_wage Real_interest_rate
Productivity;

@#endif
 

