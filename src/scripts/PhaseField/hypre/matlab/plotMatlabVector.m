function plotMatlabVector(x,rfactors,nparts,parts,labels)

b=1;
i=1;
r = [1 1];
for l=1:length(nparts)
    r = r .* rfactors(l,1:2);
    for p=1:nparts(l)
        sol(i).Low = parts(i).Low;
        sol(i).High = parts(i).High;
        sol(i).Level = l;
        dim = sol(i).High-sol(i).Low+1;
        e = b+prod(dim)-1;
        sol(i).X = ([sol(i).Low(1):sol(i).High(1)]+0.5)/r(1);
        sol(i).Y = ([sol(i).Low(2):sol(i).High(2)]+0.5)/r(2);
        sol(i).C = reshape(x(b:e),dim')';
        sol(i).Labels = string(labels(b:e));
        i=i+1;
        b = e+1;
    end
end

for l=1:length(nparts)
    hold on
    m=find([sol.Level]==l);
    for i=1:length(m)
        image(sol(m(i)).X,sol(m(i)).Y,sol(m(i)).C,'CDataMapping','scaled')
        [X,Y]=ndgrid(sol(m(i)).X,sol(m(i)).Y);
        text(reshape(X,[],1),reshape(Y,[],1),sol(m(i)).Labels,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end

end
