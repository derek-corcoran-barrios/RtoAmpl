set V;   # vertex set of the spacial graph
set E within {V,V};   # edge set of spacial graph. It should contain all loops and the pair (u,v) and (v,u) for an edge uv in the graph.
param T; # time horizon
param alpha {V,0..T};

param nchains := 8;
param u {V,0..T} ; # node biomass capacity

var y {v in V,t in 0..T} >= 0; # flow on sites
var r {E,1..T} >= 0; # inner flow on transitions

###########OPTION 1###############

minimize Quad_Cost: sum {v in V, t in 0..T} y[v,t]*y[v,t];

subj to Initial_Flow:
   sum{v in V} y[v,0] = nchains;

###########OPTION 2###############

#maximize Flow: sum {v in V} y[v,0];

################################

subj to Flow_ConservationIN {w in V, t in 1..T}:
   y[w,t] = sum {(v,w) in E} r[v,w,t];

subj to Flow_ConservationOUT {w in V, t in 0..T-1}:
   y[w,t] = sum {(w,v) in E} r[w,v,t+1];

subj to Flow_Capacity {w in V, t in 0..T}:
   y[w,t] <= u[w,t];


