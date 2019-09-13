function m = loadStructVector(fname)

m.Dim = 0;
m.NParts = 0;
m.Parts = [];
m.Periodic = [];
m.Data = [];
m.Rank = 0;
m.Full = [];

fileID = fopen(fname,'r');
l = fgetl(fileID);
if ~strcmp(l,'StructVector')
    fclose(fileID);
    error ('wrong vector format');
end

while ~feof(fileID)
    l = fgetl(fileID);
    if isempty(l)
    elseif strcmp(l,'Grid:')
        m.Dim = fscanf(fileID, '%d', 1);
        m.NParts = fscanf(fileID, '%d', 1);
        parts = fscanf(fileID, ['%d:  (' repmat('%d, ',1,m.Dim-1) '%d)  x  (' repmat('%d, ',1,m.Dim-1) '%d)'], [1+2*m.Dim,m.NParts]);
        for p = 1:m.NParts
            m.Parts(p).ID = parts(1,p);
            m.Parts(p).Low = parts(2:(1+m.Dim),p);
            m.Parts(p).High = parts((2+m.Dim):(1+2*m.Dim),p);
        end
    elseif strncmp(l,'Periodic:',9)
        m.Periodic = sscanf(l,['Periodic:' repmat(' %d',1,m.Dim)]);
    elseif strcmp(l,'Data:')
        break
    else
        fclose(fileID);
        error (['unhandeled line: "' l '"']);
    end
end

if feof(fileID)
    fclose(fileID);
    error ('no data found');
end

data = fscanf(fileID, ['%d: (' repmat('%d, ',1,m.Dim-1) '%d; %d) %f'], [3+m.Dim,Inf]);
for d = 1:size(data,2)
    m.Data(d).Part = data(1,d);
    m.Data(d).Index = data(2:(1+m.Dim),d);
    m.Data(d).Stencil = data(2+m.Dim,d);
    m.Data(d).Value = data(3+m.Dim,d);
end

fclose(fileID);

m.Rank = 0;
for p=1:m.NParts
    low = m.Parts(p).Low;
    high = m.Parts(p).High;
    dim = high-low+3; % Ghost nodes
    m.Rank = m.Rank + prod(dim);
end

m.Full = zeros(1,m.Rank);

for d = 1:length(m.Data)
    rank = 0;
    for p=1:m.NParts
        low = m.Parts(p).Low;
        high = m.Parts(p).High;
        dim = high-low+3;
        if m.Parts(p).ID==m.Data(d).Part
            index = m.Data(d).Index-low;
            if ~isempty(m.Periodic)
                index = mod(index,m.Periodic);
            end
            stride = 1;
            for i=1:m.Dim
                rank = rank + stride * (index(i)+1);
                stride = stride * dim(i);
            end
            break
        else
            rank = rank + prod(dim);
        end
    end
    m.Full(rank+1) = m.Data(d).Value;
end

end
