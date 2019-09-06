function v = loadSStructVector( nlevels, name, suffix )

v.Lenght = 0;
for i=1:nlevels
    v.StructVector(i) = loadStructVector(sprintf('%s.%02d.00%s', name, i-1, suffix));
    v.Lenght = v.Lenght+v.StructVector(i).Rank;
end
v.Full = zeros(v.Lenght,1);

s = 1;
for i=1:length(v.StructVector)
    e = s+v.StructVector(i).Rank-1;
    v.Full(s:e) = v.StructVector(i).Full;
    s = e+1;
end

end