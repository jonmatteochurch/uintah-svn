function m = loadStructMatrix(fname)

m.Symmetric = 0;
m.ConstantCoefficient = 0;
m.Dim = 0;
m.NParts = 0;
m.Parts = [];
m.Periodic = [];
m.NStencils = 0;
m.Stencils = [];
m.Data = [];
m.Rank = 0;
m.Full = [];
m.I = [];
m.J = [];

fileID = fopen(fname,'r');
l = fgetl(fileID);
if ~strcmp(l,'StructMatrix')
    fclose(fileID);
    error ('wrong matrix format');
end

while ~feof(fileID)
    l = fgetl(fileID);
    if isempty(l)
    elseif strncmp(l,'Symmetric:',10)
        m.Symmetric = sscanf(l, 'Symmetric: %d', 1);
    elseif strncmp(l,'ConstantCoefficient:',10)
        m.ConstantCoefficient = sscanf(l, 'ConstantCoefficient: %d', 1);
    elseif strcmp(l,'Grid:')
        m.Dim = fscanf(fileID, '%d', 1);
        m.NParts = fscanf(fileID, '%d', 1);
        parts = fscanf(fileID, ['%d:  (' repmat('%d, ',1,m.Dim-1) '%d)  x  (' repmat('%d, ',1,m.Dim-1) '%d)'], [1+2*m.Dim,m.NParts]);
        for p = 1:m.NParts
            m.Parts(p).ID = parts(1,p);
            m.Parts(p).Low = parts(2:(1+m.Dim),p);
            m.Parts(p).High = parts(4:(3+m.Dim),p);
        end
    elseif strncmp(l,'Periodic:',9)
        m.Periodic = sscanf(l,['Periodic:' repmat(' %d',1,m.Dim)]);
    elseif strcmp(l,'Stencil:')
        m.NStencils = fscanf(fileID, '%d', 1);
        stencils = fscanf(fileID, ['%d:' repmat(' %d',1,m.Dim)], [m.Dim+1,m.NStencils]);
        for s = 1:m.NStencils
            m.Stencils(s).ID = stencils(1,s);
            m.Stencils(s).Offset = stencils(2:(1+m.Dim),s);
        end
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
    dim = high-low+3; % Ghost nodes!!
    m.Rank = m.Rank + prod(dim);
end

m.Full = zeros(m.Rank);
m.I = zeros(m.Rank);
m.K = zeros(m.Rank);
for d = 1:length(m.Data)
    rankI = 0;
    rankJ = 0;
    for p=1:m.NParts
        low = m.Parts(p).Low-1; % Ghost nodes!!
        high = m.Parts(p).High+1; % Ghost nodes!!
        dim = high-low+1;
        if m.Parts(p).ID==m.Data(d).Part
            s = [m.Stencils.ID] == m.Data(d).Stencil;
            I = m.Data(d).Index;
            J = m.Data(d).Index+m.Stencils(s).Offset;
            if ~isempty(m.Periodic)
                I = mod(I,m.Periodic);
                J = mod(J,m.Periodic);
            end
            I = I -low;
            J = J - low;

            if any(J<-1) || any(J>dim+1)
                if m.Data(d).Value
                    msg = sprintf('ignoring data entry (%d) falling out of part\n',d);
                    msg = [msg sprintf('Part: %d\n',m.Data(d).Part')];
                    msg = [msg sprintf('Index: [%d',m.Data(d).Index(1))];
                    msg = [msg sprintf(',%d',m.Data(d).Index(2:end))];
                    msg = [msg sprintf(']\n')];
                    msg = [msg sprintf('Stencil: %d\n',m.Data(d).Stencil')];
                    msg = [msg sprintf('Value: %d\n',m.Data(d).Value')];
                    warning(msg)
                end
                break
            end
            stride = 1;
            for i=1:m.Dim
                rankI = rankI + stride * (I(i));
                rankJ = rankJ + stride * (J(i));
                stride = stride * dim(i);
            end
            break
        else
            rankI = rankI + prod(dim);
            rankJ = rankI;
        end
    end
    m.Full(rankI+1,rankJ+1) = m.Data(d).Value;
    m.I(rankI+1,rankJ+1) = I(1);
    m.J(rankI+1,rankJ+1) = I(2);
    if m.Symmetric
        m.Full(rankJ+1,rankI+1) = m.Full(rankI+1,rankJ+1);
    end
end

end
