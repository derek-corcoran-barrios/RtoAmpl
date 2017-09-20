set Sp; # set de especies
set V;   # vertex set of the spacial graph
set E within {V,V};   # edge set of spacial graph. It should contain all loops and the pair (u,v) and (v,u) for an edge uv in the graph.
param T; # time horizon

param alpha {Sp, V,0..T}; # biomass amplification/reduction factor 
#param beta {E,0..T-1}; # biomass transition factor
param grado {V};        # grado de cada nodo
param beta {Sp, E,0..T-1};  # biomass transition factor
param bf {Sp};

param b0 {Sp, V}; # initial biomass
#param bf {V}; # final (required) biomass 
param c {V}; # node cost
param u {Sp, V,0..T}; # node biomass capacity

#var z {V} binary;  # buying decision for each node in V
var z {V};
var y {Sp, V,0..T} >= 0; # biomass flow

minimize Buy_Cost:
   sum {v in V } c[v] * z[v];

subj to Initial_Flow {s in Sp,w in V}:
   y[s, w,0] <= b0[s, w];

#subj to Final_Demand {s in Sp,w in V}:
 #  y[s,w,T] >= bf[s,w];
 
#minimum demand and the end of period T   
subj to Final_Demand_Total {s in Sp}:
   sum {w in V} y[s,w,T] >= bf[s];   

subj to Flow_Conservation {s in Sp, w in V, t in 0..T-1}:
   y[s,w,t+1] <= sum {(v,w) in E} beta[s,v,w,t]*alpha[s,v,t]*y[s,v,t];

subj to Flow_Capacity {s in Sp,w in V, t in 0..T}:
   y[s,w,t] <= u[s,w,t] * z[w];

subj to zeroone {w in V}:
 0 <= z[w] <=1;
#display sum {v in V} y[v,T]*alpha[v,T];
#si z = 1, en todo, maximizar suma de bf
