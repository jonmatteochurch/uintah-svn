function s = loadSStructSystem(nlevels, fname, Aname, bname, xname, suffix)

s.A = loadSStructMatrix(nlevels,sprintf('%s/%s',fname,Aname),suffix);
s.b = loadSStructVector(nlevels,sprintf('%s/%s',fname,bname),suffix);
s.x = loadSStructVector(nlevels,sprintf('%s/%s',fname,xname),suffix);
s.Ranks = 0:length(s.x.Full)-1;

zero_rows=find(all(s.A.Full==0,2));
zero_cols=find(all(s.A.Full==0,1));

nghosts=0;
g=[0,0,0];
g(1:s.A.StructMatrix(1).Dim) = 2;
for m=1:nlevels
    for p=1:s.A.StructMatrix(m).NParts
        dim = s.A.StructMatrix(m).Parts(p).High-s.A.StructMatrix(m).Parts(p).Low+1;
        nghosts = nghosts + prod(dim+g) - prod(dim);
    end
end

if any(zero_rows~=zero_cols')
    warning('somethingwrong with ghost nodes')
end

if length(zero_rows)~=nghosts
    warning('somethingwrong with ghost nodes')
end

s.A0 = s.A.Full;
s.b0 = s.b.Full;
s.x0 = s.x.Full;
s.Ranks0 = s.Ranks;

s.A0(zero_rows,:) = [];
s.A0(:,zero_cols) = [];
s.b0(zero_rows) = [];
s.x0(zero_rows) = [];
s.Ranks0(zero_rows) = [];

end