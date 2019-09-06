function plotSStructSystem(s,rfactors,titl)

figure
title ([titl ': A (w/ghosts)'])
spy (s.A.Full)

figure
title ([titl ': x (hypre)'])
plotSStructVector(s.x,rfactors,s.Ranks)

figure
title ([titl ': x (matlab)'])
title('matlab')

x0 = s.A0 \ s.b0;
norm (s.x0 - x0)
plotMatlabVector(x0,rfactors,[s.A.StructMatrix(:).NParts],[s.A.StructMatrix(:).Parts],s.Ranks0)
