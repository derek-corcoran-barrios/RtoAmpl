set V;   # vertex set of the spacial graph
set E within {V,V};   # edge set of spacial graph. It should contain all loops and the pair (u,v) and (v,u) for an edge uv in the graph.
set SP; #Set of Species

param T; # time horizon

param nchains{SP} >= 0;
param u {SP, V,0..T} ; # node biomass capacity for species SP, in node V in time T

var y {s in SP,v in V,t in 0..T} >= 0; # flow on sites
var r {SP,E,1..T} >= 0; # inner flow on transitions

###########OPTION 1###############

minimize Quad_Cost: sum {s in SP, v in V, t in 0..T} y[s,v,t]*y[s,v,t];

subj to Initial_Flow{s in SP}:
   sum{v in V} y[s,v,0] = nchains[s];

###########OPTION 2###############

#maximize Flow: sum {v in V} y[v,0];

################################

subj to Flow_ConservationIN {s in SP, w in V, t in 1..T}:
   y[s,w,t] = sum {(v,w) in E} r[s,v,w,t];

subj to Flow_ConservationOUT {s in SP,w in V, t in 0..T-1}:
   y[s,w,t] = sum {(w,v) in E} r[s,w,v,t+1];

subj to Flow_Capacity {s in SP, w in V, t in 0..T}:
   y[s,w,t] <= u[s,w,t];

