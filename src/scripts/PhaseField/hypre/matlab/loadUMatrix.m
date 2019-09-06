function m = loadStructMatrix(fname)

m.Low = [];
m.High = [];
m.Data = [];
m.Size = [];
m.Full = [];

fileID = fopen(fname,'r');
if (fileID<0)
    error (['cannot open ' fname])
end
dim = fscanf(fileID, '%d %d %d %d', 4);
m.Low = dim([1,3]);
m.High = dim([2,4]);

data = fscanf(fileID, '%d %d %f', [3, Inf]);
for d = 1:size(data,2)
    m.Data(d).I = data(1,d);
    m.Data(d).J = data(2,d);
    m.Data(d).Value = data(3,d);
end

fclose(fileID);

m.Size = m.High-m.Low+1;
m.Full=zeros(m.Size');
for d = 1:length(m.Data)
    I = m.Data(d).I - m.Low(1) + 1;
    J = m.Data(d).J - m.Low(2) + 1;
    m.Full(I,J) = m.Data(d).Value;
end

end
