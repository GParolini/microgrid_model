$title Microgrid model for a research data centre

Set
    t time period t /t0*t8759/
    k efficiency measure k /k0,k1,k2/
    o 'objective functions' / costs, emissions /
    points  'pareto points' / point1*point1000 /;

$set min -1
$set max +1

Parameter dir(o) 'direction of the objective functions 1 for max and -1 for min'
                 / costs %min%, emissions %min% /;
                 
Parameter pareto_obj(points,o) 'objective values of the pareto points';

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

Parameter Em(t) Grid carbon emissions at time t ;
$CALLTOOL csvread Em.csv id=Em dimids=Index index=1 values=lastCol useHeader="Y" gdXout=Em.gdx
$GDXIN Em.gdx
$LOAD Em
$GDXIN

Parameter Cm(t) Price for one kWh purchased from the grid at time t ;
$CALLTOOL csvread Cm.csv id=Cm dimids=Index index=1 values=lastCol useHeader="Y" gdXout=Cm.gdx
$GDXIN Cm.gdx
$LOAD Cm
$GDXIN

Parameter Vv(k) Energy reduction for efficiency measure k
                /k0   0.0
                 k1   0.15
                 k2   0.25/;

Parameter Cc(k) Amortised yearly cost of adopting efficiency measure k
                /k0       0
                 k1    5000
                 k2   10000/;


Scalars
    TT      Length of time period in hours                                              /1/
    Epv     Photovoltaic system lifecycle emissions (per year)                          /25.86/
    Ebess   Battery lifecycle emissions (per year)                                      /12.06/
    Cpv     Amortised cost of one kW of installed photovoltaic capacity (per year)      /99/
    A       Maximum PV capacity                                                         /606/
    Cbess   Amortised cost of one kWh of installed battery capacity (per year)          /40/
    Bb      Maximum BESS capacity                                                       /375/
    F       Battery charging efficiency factor                                          /0.9/
    D       Battery self-discharge                                                      /0.003/
    Sma     Battery max state-of-charge                                                 /0.8/
    Smi     Battery min state-of-charge                                                 /0/ 
    Gmax    Battery max rate of charging                                                /0.3/
    Gmin    Battery max rate of discharging                                             /0.3/                   
;

Variables
    m(o)    Objective function variables
    x(t)    Grid electricity purchased at time t
    y(t)    PV electricity harnessed at time t
    r(t)    PV electricity used by the DC at time t
    h(t)    PV electricity used to charge BESS at time t
    v(k)    Binary variable for efficiency measure k
    b(t)    BESS electricity used by the DC at time t
    e(t)    BESS storage level at time t
    z       BESS installed capacity
    w       PV installed capacity
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
        emissions           total emissions associated to the data-centre operations--OBJECTIVE
        costs               total costs associated to the data-centre operations--OBJECTIVE
        dc_balance          data-centre energy balance
        pv_capacity1        limits on energy sent by pv to dc-bess-grid by the photovoltaic system
        pv_capacity2        limits on the electricity produced by the photovoltaic system
        bess_balance        energy stored in the battery
        bess_max_soc        battery max SOC
        bess_min_soc        battery min SOC
        bess_charge         battery charging
        bess_discharge      battery discharging
        space_restr_pv      area restrictions for PV installation
        space_restr_bess    volume restrictions for BESS installation
        eff_measures        exactly one efficiency measure may be adopted in the microgrid
;

emissions ..               m('emissions') =e= sum(t, Em(t)*x(t)) + Epv*w + Ebess*z ;
costs ..                   m('costs') =e= sum(t, Cm(t)*x(t)) + Cpv*w + Cbess*z + sum(k, Cc(k)*v(k)) ;
dc_balance(t) ..           x(t) + r(t) + b(t) =e= DC(t)*(sum(k,(1-Vv(k))*v(k)))  ;         
pv_capacity1(t) ..         r(t) + h(t) =l= y(t) ;
pv_capacity2(t) ..         y(t) =l= w*P(t)*TT  ;   
bess_balance(t) ..         e(t) =e= (1-D)*e(t-1)$(ord(t)>1) + h(t)*F - b(t) ;        
bess_max_soc(t) ..         e(t) =l= Sma*z ;
bess_min_soc(t) ..         (Smi*z)$(ord(t)>1)  =l= e(t) ;      
bess_charge(t) ..          h(t)*F =l= Gmax*z   ;     
bess_discharge(t) ..       b(t) =l= Gmin*z  ;
space_restr_pv ..          w =l= A ;
space_restr_bess ..        z =l= Bb ;
eff_measures ..            sum(k, v(k)) =e= 1 ;


Model microgrid /all/ ;

Set oo(o) active objective function;
oo(o) = yes;

$onEcho > cplex.opt
threads 1
$offEcho
;

$libInclude moo EPSConstraint microgrid MIP o dir m points pareto_obj -iterations=10 -gridpoints=50 -savepoint=1 -savepoint_filename= -savepoint_dir=savepoints -solver=cplex  -optfile_init=1 -optfile_main=1 
execute 'gdxmerge savepoints%system.DirSep%*.gdx > %system.NullFile%';

display pareto_obj;
