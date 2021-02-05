function postprocess(symlink,lvl)
[~,name] = system(sprintf('readlink %s',symlink));
assert(~isempty(name),'empty symlink');

delim=find(name=='.',1,'last');
n=str2double(name(delim+1:end));
name=name(1:delim-1);
assert(endsWith(name,'.uda'),'non uda symlink');

pos=zeros(0,2);
vel=zeros(0,2);
cur=zeros(0,4);
amr=struct([]);

for i=n:-1:0
    dat=load(sprintf('%s.%03d/tip_position.dat',name,i));
    if ~isempty(pos)
        tend=find(dat(:,1)==pos(1,1))-1;
    else
        tend=length(dat(:,1));
    end
    pos=[dat(1:tend,:);pos];
    
    dat=load(sprintf('%s.%03d/tip_velocity.dat',name,i));
    vel=[dat(1:tend,:);vel];
    
    dat=load(sprintf('%s.%03d/tip_curvatures.dat',name,i));
    cur=[dat(1:tend,:);cur];
    
    cmd=sprintf('grep -H "DataArchiver created %s.%03d" %s.arc4.o* | cut -d : -f 1',name,i,name(1:end-4));
    [~,log]=system(cmd);
    log=sscanf(log,'%s');
    dat=loadlvls(log,lvl);
    if ~isempty(amr)
        tend=find([dat.timestep]>=amr(1).timestep)-1;
    else
        tend=length(dat);
    end
    amr=[dat(1:tend) amr];
end

dir=sprintf('%s.full',name(1:end-4));
mkdir(dir);

file=fopen(sprintf('%s/tip_position.dat',dir),'w');
fmt='%.13f\t%.13f\n';
fprintf(file,fmt,pos');
fclose(file);

file=fopen(sprintf('%s/tip_velocity.dat',dir),'w');
fmt='%.13f\t%.13f\n';
fprintf(file,fmt,vel');
fclose(file);

file=fopen(sprintf('%s/tip_curvatures.dat',dir),'w');
fmt='%.13f\t%.13f\t%.13f\t%.13f\n';
fprintf(file,fmt,cur');
fclose(file);

file=fopen(sprintf('%s/tip_curvatures.dat',dir),'w');
fmt='%.13f\t%.13f\t%.13f\t%.13f\n';
fprintf(file,fmt,cur');
fclose(file);

file=fopen(sprintf('%s/total_cells.dat',dir),'w');
fmt='%.13f\t%d\n';
T=[amr.time];
C=[amr.totalcells];
fprintf(file,fmt,[T;C]);
fclose(file);

end

function dat = loadlvls(file,lvls)
if(lvls>1)
    cmd=sprintf('grep "Timestep 0" %s > tmp; egrep "Patches|Timestep" %s | grep -B 1 Patches >> tmp', file, file);
    system(cmd);
    fid = fopen('tmp');
    dat = [];
    while ~feof(fid)
        line = fgetl(fid);
        if(line(1)=='-')         
            line = fgetl(fid);
        end
        ts = textscan(line,'Timestep %d   Time=%f     Next delT=%f       Wall Time=%f     EMA=%f     Memory Use=%f MBs');
        d.timestep = ts{1};
        d.time = ts{2};
        d.delt = ts{3};
        d.wt = ts{4};
        d.ema = ts{5};
        d.mem = ts{6};
        d.lvls = [];
        if d.timestep==0
        for m=max(2,lvls-4):lvls-1
        for l=1:m
            line = fgetl(fid);
        end
        end
        end
        for l=1:lvls
            line = fgetl(fid);
            lv = textscan(line,'  L%d      RefineRatio: [int %d, %d, %d] Patches: %d     Total Cells: %d  Mean Cells: %f stdv: %f relative stdv: %f Volume: %f');
            lvl.index = lv{1};
            lvl.ratio = [lv{2:4}];
            lvl.patches = lv{5};
            lvl.cells = lv{6};
            lvl.meancells = lv{7};
            lvl.stdv = lv{8};
            lvl.relstdv = lv{9};
            lvl.vol= lv{10};
            d.lvls = [d.lvls, lvl];
        end
        d.totalcells = sum([d.lvls.cells]);
        dat = [dat, d];
    end
else
    system(sprintf('egrep "Total" %s > tmp', file));
    fid = fopen('tmp');
    dat = [];
    d.timestep = 0;
    d.time = 0;
    d.delt = NaN;
    d.wt = NaN;
    d.ema = NaN;
    d.mem = NaN;
    d.lvls = [];
    lvl.index = 1;
    lvl.ratio = [1 1 1];
    line = fgetl(fid);
    lv = textscan(line,'Total number of cells:        %d (%f avg. per patch)');
    lvl.cells = lv{1};
    lvl.meancells = lv{2};
    line = fgetl(fid);
    lv = textscan(line,'Total patches in grid:          %d');
    lvl.patches = lv{1};
    lvl.stdv = NaN;
    lvl.relstdv = NaN;
    lvl.vol= NaN;
    d.lvls = [d.lvls, lvl];
    d.totalcells = sum(d.lvls.cells);
    dat = [dat, d];
end
	system('rm tmp');
end
