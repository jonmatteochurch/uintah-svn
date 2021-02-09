close all
ls='-:';
pb='karma3d';

VV=[0.09201 0.07088 0.09167 0.07563 0.04200 0.01321];
ww=ceil(sqrt(3)*[362 627 485 588 1058 2803]);
kk=1./[27.799 20.340 13.525 9.366 7.458 7.025];

names={};

regs=[];
vels=[];
verrs=[];
vrels=[];
crvs=[];
kerrs=[];
krels=[];
avgs=[];
mins=[];
maxs=[];

for p=1:6
	h=(p-1)/6;
	
	fprintf('%s_%d\n',pb,p);
	names{end+1,1}=sprintf('%s_%d\n',pb,p);
	dat = tmpload(sprintf('%s_%d.full/tip_position.dat',pb,p));

	b = find(dat(:,1)==0,1,'last');
	e = min([find(dat(:,1)>1800,1), size(dat,1)]);

	figure(1); hold on;
	title('tip position');
	plot(dat(b:e,1),movmean(dat(b:e,2),ww(p)),'LineStyle','-','Color',hsv2rgb([h,1,1]),'DisplayName',sprintf('%s %d',pb,p));
	xlim([0,1800])
	
	dat = tmpload(sprintf('%s_%d.full/tip_velocity.dat',pb,p));
	
	figure(2); hold on;
	title('tip velocity');
	plot(dat(b:e,1),movmean(dat(b:e,2),ww(p)),'LineStyle','-','Color',hsv2rgb([h,1,1]),'DisplayName',sprintf('%s %d',pb,p));
	plot([0,1800],[VV(p) VV(p)],'LineStyle',':','Color',hsv2rgb([h,1,1]),'DisplayName',sprintf('%s GF %d',pb,p));
	xlim([0,1800])
	ylim([0,.4]);
	
	V=mean(dat(e-ww(p):e,2));
	vels=[vels;V];
	verrs=[verrs;abs(V-VV(p))];
	vrels=[vrels;abs(V-VV(p))/VV(p)];

	dat = tmpload(sprintf('%s_%d.full/tip_curvatures.dat',pb,p));
	
	figure(3); hold on;
	title('tip curvature');
	plot(dat(b:e,1),movmean(dat(b:e,2),ww(p)),'LineStyle','-','Color',hsv2rgb([h,1,1]),'DisplayName',sprintf('%s %d',pb,p));
	plot([0,1800],[kk(p) kk(p)],'LineStyle',':','Color',hsv2rgb([h,1,1]),'DisplayName',sprintf('%s GF %d',pb,p));
	xlim([0,1800])
%	ylim([0,40]);
	
	k=mean(dat(e-ww(p):e,2));
	crvs=[crvs;k];
	kerrs=[kerrs;abs(k-kk(p))];
	krels=[krels;abs(k-kk(p))/kk(p)];

	figure(4); hold on;
	title('parabolic curvature');
	plot(dat(b:e,1),movmean(dat(b:e,4),ww(p)),'LineStyle','-','Color',hsv2rgb([h,1,1]),'DisplayName',sprintf('%s %d',pb,p));
	xlim([0,1800])
	ylim([0,.1]);
	
	tend=dat(e,1);
	dat = tmpload(sprintf('%s_%d.full/total_cells.dat',pb,p));
	e = min([find(dat(:,1)>1800,1), size(dat,1)]);
	
	T=[dat(b:e,1);tend];
	dT=T(2:end)-T(1:end-1);
	
	figure(5); hold on;
	title('total cells');
	plot(T,[dat(b:e,2);dat(e,2)],'LineStyle','-','Color',hsv2rgb([h,1,1]),'DisplayName',sprintf('%s %d',pb,p));
	set(gca,'YScale','log');
	xlim([0,1800])
	regs=[regs;e];
	avgs=[avgs;sum(dT.*dat(b:e,2))/tend];
	mins=[mins;min(dat(b:e,2))];
	maxs=[maxs;max(dat(b:e,2))];
end

for f=1:4
    figure(f);
	legend('off');
    epsprint(14,10);
end

figure(5);
legend('Location','northeastoutside');
epsprint(28,10);

varNames={'Sim','Tip Velocity','VAbsError','VRelError','Tip Curvature','KAbsError','KRelError','Regrids','AvgTotCells','MinTotCells','MaxTotCells'}
T=table(names,vels,verrs,vrels,crvs,kerrs,krels,regs,avgs,mins,maxs,'VariableNames',varNames)
writetable(T,'summary.txt');

function dat = tmpload(file)
% system(sprintf('cp %s tmp', file));
% dat = load('tmp');
dat = load(file);
end

function epsprint(width,height)
% Load file
h1=gca;
name=strrep(h1.Title.String,' ','_')

% Change size
set(gcf,'paperunits','centimeters')
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'papersize',[width,height])
set(gcf,'paperposition',[0,0,width,height])
set(gcf, 'renderer', 'painters');

% Export file
print('-depsc2',name);
end
