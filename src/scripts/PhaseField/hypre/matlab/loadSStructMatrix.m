function A = loadSStructMatrix ( nlevels, name, suffix )

global A

A.UMatrix = loadUMatrix([name '.UMatrix' suffix]);
A.Full = A.UMatrix.Full;
A.Size = A.UMatrix.Size;
A.Ghosts = [];

s = 1;
for i=1:nlevels
    A.StructMatrix(i) = loadStructMatrix(sprintf('%s.%02d.00.00%s', name, i-1, suffix));
    e = s+A.StructMatrix(i).Rank-1;
    A.Full(s:e,s:e) = A.StructMatrix(i).Full;
    s = e+1;
end

rank0 = 0;
for l=1:nlevels
    rankl = rank0;
    for p=1:A.StructMatrix(l).NParts
        low = A.StructMatrix(l).Parts(p).Low-1;
        high = A.StructMatrix(l).Parts(p).High+1;
        dim = high-low+1;
        switch A.StructMatrix(l).Dim
            case 2
                for x=[0 dim(1)-1]
                    for y=0:dim(2)-1
                        process_ghost([x;y],rank0,rankl,l,p,low,high,dim);
                    end
                end
                for y=[0 dim(2)-1]
                    for x=0:dim(1)-1
                        process_ghost([x;y],rank0,rankl,l,p,low,high,dim);
                    end
                end
            case 3
                for x=[0 dim(1)-1]
                    for y=0:dim(2)-1
                        for z=0:dim(3)-1
                            process_ghost([x;y;z],rank0,rankl,l,p,low,high,dim);
                        end
                    end
                end
                for x=0:dim(1)-1
                    for y=[0 dim(2)-1]
                        for z=0:dim(3)-1
                            process_ghost([x;y;z],rank0,rankl,l,p,low,high,dim);
                        end
                    end
                end
                for x=0:dim(1)-1
                    for y=0:dim(2)-1
                        for z=[0 dim(3)-1]
                            process_ghost([x;y;z],rank0,rankl,l,p,dim);
                        end
                    end
                end
        end
        rankl = rankl + prod(dim);
    end
    rank0 = rankl;
end

end

function process_ghost(I, rank0, rankl, l,p,low, high, dim)
global A ghosts
stride = 1;
rank = rankl;
for d=1:A.StructMatrix(l).Dim
    rank = rank + stride * (I(d));
    stride = stride * dim(d);
end
if (any(A.Full(rank+1,:)))
    error('non null ghost row!')
end
if (any(A.Full(:,rank+1)))
    ind = low+I;
    rank2 = rank0;
    I2 = [];
    for p2=1:A.StructMatrix(l).NParts
        low2 = A.StructMatrix(l).Parts(p2).Low-1;
        high2 = A.StructMatrix(l).Parts(p2).High+1;
        dim2 = high2-low2+1;
        if p2 ~= p & low2 < ind & ind < high2
            I2 = ind - low2;
            stride2 = 1;
            for d=1:A.StructMatrix(l).Dim
                rank2 = rank2 + stride2 * (I2(d));
                stride2 = stride2 * dim2(d);
            end
            break
        else
            rank2 = rank2 + prod(dim2);
        end
    end
    A.Full(:,rank2+1) = A.Full(:,rank2+1) + A.Full(:,rank+1);
    A.Full(:,rank+1) = 0;
    A.ghosts = [A.ghosts; l p I' rank p2 I2' rank2];
end
end