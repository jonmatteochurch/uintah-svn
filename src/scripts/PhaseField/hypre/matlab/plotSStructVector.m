function plotSStructVector(x,rfactors,labels)

i=1;
rank0=1;
r = [1 1];
for l=1:length(x.StructVector)
    r = r .* rfactors(l,1:2);
    for p=1:x.StructVector(l).NParts
        sol(i).Low = x.StructVector(l).Parts(p).Low;
        sol(i).High = x.StructVector(l).Parts(p).High;
        sol(i).Level = l;
        dim = sol(i).High-sol(i).Low+3;
        rank1 = rank0+prod(dim)-1;
        sol(i).X = ([sol(i).Low(1):sol(i).High(1)]+0.5)/r(1);
        sol(i).Y = ([sol(i).Low(2):sol(i).High(2)]+0.5)/r(2);
        c = reshape(x.Full(rank0:rank1),dim')';
        c([1 end],:) = [];
        c(:,[1 end]) = [];
        sol(i).C = c;
        lbl = reshape(labels(rank0:rank1),dim');
        lbl([1 end],:) = [];
        lbl(:,[1 end]) = [];
        sol(i).Labels = string(reshape(lbl,[],1));
        i=i+1;
        rank0=rank1+1;
    end
end

for l=1:length(x.StructVector)
    hold on
    m=find([sol.Level]==l);
    for i=1:length(m)
        image(sol(m(i)).X,sol(m(i)).Y,sol(m(i)).C,'CDataMapping','scaled')
        [X,Y]=ndgrid(sol(m(i)).X,sol(m(i)).Y);
        text(reshape(X,[],1),reshape(Y,[],1),sol(m(i)).Labels,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
