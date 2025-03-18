$title Microgrid model for a research data centre



Sets
    t time period t /t0*t8759/
    k efficiency measure k /0,1,2/;

Parameter DC(t) Data-centre load ;
$CALLTOOL csvread DC.csv id=DC dimids=Index index=1 values=lastCol useHeader="Y" gdXout=DC.gdx
$GDXIN DC.gdx
$LOAD DC
$GDXIN

Parameter P(t) Performance factor of the photovoltaic system ;
$CALLTOOL csvread P.csv id=P dimids=Index index=1 values=lastCol useHeader="Y" gdXout=P.gdx
$GDXIN P.gdx
$LOAD P
$GDXIN

Parameter Vv(k) Energy reduction for efficiency measure k
                /0   0.0
                 1   0.1
                 2   0.2/;

Parameter Cc(k) Cost of adopting efficiency measure k
                /0   0
                 1   100000
                 2   200000/;


Scalars
    TT      Length of time period in hours                                              /1/
    Em      Grid emissions                                                              /0.6038/
    Epv     Photovoltaic system lifecycle emissions (per year)                          /25.86/
    Ebess   Battery lifecycle emissions (per year)                                      /12.06/
    Cm      Price for one kWh purchased from the grid                                   /11.76/
    Cpv     Amortised cost of one kW of installed photovoltaic capacity (per year)      /650/
    A       Maximum PV capacity                                                         /1000/
    Cbess   Amortised cost of one kWh of installed battery capacity (per year)          /240/
    Bb      Maximum BESS capacity                                                       /500/
    F       Battery charging efficiency factor                                          /0.9/
    D       Battery self-discharge                                                      /0.003/
    Gmax    Battery max rate of charging                                                /0.3/
    Gmin    Battery max rate of discharging                                             /0.3/                   
;

Variables
    x(t)    Grid electricity purchased at time t
    y(t)    PV electricity harnessed at time t
    r(t)    PV electricity used by the DC at time t
    h(t)    PV electricity used to charge BESS at time t
    v(k)    Binary variable for efficiency measure k
    b(t)    BESS electricity used by the DC at time t
    e(t)    BESS storage level at time t
    z       BESS installed capacity
    w       PV installed capacity
    o       total emissions
    c       total costs
;


Binary Variable v(k)    ;


Positive variables
    x
    y
    r
    h
    b
    e
    w
    z
;



Equations
        emissions           total emissions associated to the data-centre operations
        costs               total costs associated to the data-centre operations
        dc_balance          data-centre energy balance
        pv_capacity1        limits on energy sent by pv to dc-bess-grid by the photovoltaic system
        pv_capacity2        limits on the electricity produced by the photovoltaic system
        bess_balance        energy stored in the battery
        bess_limits         limits to the energy that can be stored in the battery
        bess_ramping_up     max SOC
        bess_ramping_down   min SOC
        space_restr_pv      area restrictions for PV installation
        space_restr_bess    volume restrictions for BESS installation
        eff_measures        exactly one efficiency measure may be adopted in the microgrid
;

emissions ..               o =e= sum(t, Em*x(t)) + Epv*w + Ebess*z ;
costs ..                   c =e= sum(t, Cm*x(t)) + Cpv*w + Cbess*z + sum(k, Cc(k)*v(k)) ;
dc_balance(t) ..           x(t) + r(t) + b(t) =e= DC(t)*(sum(k,(1-Vv(k))*v(k)))  ;         
pv_capacity1(t) ..         r(t) + h(t) =l= y(t) ;
pv_capacity2(t) ..         y(t) =l= w*P(t)*TT  ;   
bess_balance(t) ..         e(t) =e= (1-D)*e(t-1)$(ord(t)>1) + h(t)*F - b(t) ;        
bess_limits(t) ..          e(t) =l= z  ;        
bess_ramping_up(t) ..      h(t)*F =l= Gmax*z   ;     
bess_ramping_down(t) ..    b(t) =l= Gmin*z  ;
space_restr_pv ..          w =l= A ;
space_restr_bess ..        z =l= Bb ;
eff_measures ..            sum(k, v(k)) =e= 1 ;

option optcr=0  ;
option reslim = 3000000 ;
Option Iterlim=100000   ;
Option MIP = CPLEX  ;

Model microgrid /all/ ;
Solve microgrid using MIP minimizing o ;
Solve microgrid using MIP minimizing c ;
