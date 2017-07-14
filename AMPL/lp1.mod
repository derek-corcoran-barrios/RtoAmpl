set V;   # vertex set of the spacial graph
set E within {V,V};   # edge set of spacial graph. It should contain all loops and the pair (u,v) and (v,u) for an edge uv in the graph.
param T; # time horizon

param alpha {V,0..T}; # biomass amplification/reduction factor 
#param beta {E,0..T-1}; # biomass transition factor
param grado {V};        # grado de cada nodo
param beta {E,0..T-1};  # biomass transition factor
param bf := 400;

param b0 {V}; # initial biomass
#param bf {V}; # final (required) biomass 
param c {V}; # node cost
param u {V,0..T}; # node biomass capacity

var z {V} binary;  # buying decision for each node in V
var y {V,0..T} >= 0; # biomass flow

minimize Buy_Cost:
   sum {v in V } c[v] * z[v];

subj to Initial_Flow {w in V}:
   y[w,0] <= b0[w];

#subj to Final_Demand {w in V}:
 #  y[w,T] >= bf[w];
 
#minimum demand and the end of period T   
subj to Final_Demand_Total:
   sum {w in V} y[w,T]*alpha[w,T] >= bf;   

subj to Flow_Conservation {w in V, t in 0..T-1}:
   y[w,t+1] <= sum {(v,w) in E} beta[v,w,t]*alpha[v,t]*y[v,t];

subj to Flow_Capacity {w in V, t in 0..T}:
   y[w,t] <= u[w,t] * z[w];

#display sum {v in V} y[v,T]*alpha[v,T];
