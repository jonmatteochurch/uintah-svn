function plotSStructSystem(s,rfactors,titl)

figure
spy (s.A.Full)
title ([titl ': A (w/ghosts)'])

figure
plotSStructVector(s.x,rfactors,s.Ranks)
title ([titl ': x (hypre)'])

x0 = s.A0 \ s.b0;
norm (s.x0 - x0)

figure
plotMatlabVector(x0,rfactors,[s.A.StructMatrix(:).NParts],[s.A.StructMatrix(:).Parts],s.Ranks0)
title ([titl ': x (matlab)'])
