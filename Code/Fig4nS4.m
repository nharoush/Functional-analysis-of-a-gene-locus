%% Figure 4  -  heterozygotes
clearvars;
if ~exist('processedData', 'dir')
    mkdir('processedData')
end
fol=cd;
pn = fullfile(fol,'processedData');
gapNames={'Hb','Kr','Gt','Kni'};

% get gap at eve positions data arrays:
t4checkGap=42;
dt42WTm=nan(4,7,4);%each cell contains one line data in the form: 4 gene, n embryos, 7stripes
dt42Hetsm=nan(4,7,4);%same as dt42WT
dt42WTatHetsXm=nan(4,7,4);

 for lineID=1:4
  fn1=fullfile(pn,['PRpeakNvalley4',gapNames{lineID}]);   
  load(fn1);% these are stripe positions in individual embryos, with an 8 min mean delay
  fn2=fullfile(pn,[gapNames{lineID},'DiffData200runs5']); 
  load (fn2); 

  dtAllLines(:,:,:,:, lineID)=dt;
  dtAllLinesSTE(:,:,:,:, lineID)=squeeze(nanstd(dtA,[], 2));

    x4checkWt{lineID}=round(nanmean(EvePeakLocWT))-100;% these are the mean eve peak positions at wt
    DeltaX(lineID, :)=nanmean(EvePeakLocHets)-nanmean(EvePeakLocWT);
    deltaXposEr(lineID,:)=nanstd(EvePeakLocWT);
    deltaXposErHets(lineID,:)=nanstd(EvePeakLocHets);
    x4checkHets=round(nanmean(EvePeakLocHets))-100;% these are the mean eve peak positions at hets
    
    dt42WTm(:,:,lineID)=squeeze(dtAllLines(:, t4checkGap-9,x4checkWt{lineID},3,lineID));
    dt42Hetsm(:,:,lineID)=squeeze(dtAllLines(:, t4checkGap-9,x4checkHets,2,lineID));
    dt42WTste(:,:,lineID)=squeeze(dtAllLinesSTE(:, t4checkGap-9,x4checkWt{lineID},3,lineID));
    dt42Hetsste(:,:,lineID)=squeeze(dtAllLinesSTE(:, t4checkGap-9,x4checkHets,2,lineID));
        
 end
StripeClassNew=(abs(DeltaX)>=deltaXposEr);

StripeClassNew=StripeClassNew+1; % this way stripes that shift are class 2 and those that don't are 1

%% plot Panels A&B
% plot the Intensity comparison of gap activation levels for hets stripes :
    MS=4;% Marker size for plotting
    maxI=1.3;
    Delta4Noise=0.1;%0.08;
  Mark=['o', 's', 'd', '^'];
left=0.15;
bottom=0.15;
width=6.8;
hight=2;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
  C=colormap(jet(7));
for s=[3:5]%1:7
    for rgb=1:3
    if (C(s,rgb))~=0
        C(s,rgb)=C(s,rgb)-0.1;
    end
    end
end
  c=[0.7 ,0 ,0.6; 0.13,0.5 ,0.15];
    
    subplot(1,2,1);
    area([0,maxI], [Delta4Noise,maxI+Delta4Noise], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
    hold on
    area([0,maxI], [0,maxI], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
    hold on
    area([0,maxI], ([ -1*Delta4Noise, maxI-Delta4Noise]), 'FaceAlpha', 1,'FaceColor', [1,1,1], 'LineStyle', 'none');
    plot([0 maxI], [0,maxI], '--k');
    plot([0,maxI],[0,0.5*maxI], '--k');
    xlabel(['I@\it{2xGap}']);
    ylabel(['I@\it{1xGap}']);        
    box off
    axis  square
    xlim([0 maxI]);
    ylim([0 maxI]);
    
    for lineID=1:4
        for s=1:7
            if StripeClassNew(lineID,s)==1
                for g=1:4
                    if g==lineID
                        CC='r';
                    else
                        CC='k';
                    end
                    errorbar(dt42WTm(g,s,lineID),dt42Hetsm(g,s,lineID),...
                        dt42Hetsste(g,s,lineID),dt42Hetsste(g,s,lineID),...
                        dt42WTste(g,s,lineID),dt42WTste(g,s,lineID),...
                        'Marker', Mark(1), 'Color',CC,'LineStyle','none',...
                        'MarkerSize', MS);
                    hold on
                end
            end
        end
    end
    hold on
    for lineID=1:4
        for s=1:7
            if StripeClassNew(lineID,s)==2
                for g=1:4
                    if g==lineID
                        CC='r';
                    else
                        CC='g';
                    end
                    errorbar(dt42WTm(g,s,lineID),dt42Hetsm(g,s,lineID),...
                        dt42Hetsste(g,s,lineID),dt42Hetsste(g,s,lineID),...
                        dt42WTste(g,s,lineID),dt42WTste(g,s,lineID),...
                        'Marker', Mark(1), 'Color',CC,'LineStyle','none',...
                        'MarkerSize', MS);
                    hold on
                end
            end
        end
    end
    
    p1=errorbar(nan, nan,5, 'Marker','o','LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    hold on
    p2=errorbar(nan, nan,5, 'Marker','o','LineStyle', 'none', 'Color', 'g' , 'MarkerSize',MS);
    p3=errorbar(nan, nan,5, 'Marker','o','LineStyle', 'none', 'Color', 'r' , 'MarkerSize',MS);
    legend([p1,p2,p3],{'\DeltaX<=\deltaX_{wt}','\DeltaX>\deltaX_{wt}' ,'perturbed Gap'});
    set(gca, 'FontSize', 8);
    % now plot the relationship between DeltaI and Deltax in 1xgap: the more a nucleus increase a concentration 
    
    EvePXhets=nan(7,4);
EveTXhets=nan(6,4);
EvePXwt=nan(7,4);
EveTXwt=nan(6,4);
DeltaXP=nan(7,4);
DeltaXT=nan(6,4);
tmax=42;%42;
tmin=10;%22;

dIp=nan(4,4,7);%(linE,g,s)
dIpmin=nan(4,4,7);%(linE,g,s)
dIt=nan(4,4,6);%(linE,g,s)

subplot(1,2,2);
for linE=1:4
    fn=fullfile(pn,['PRpeakNvalley4', gapNames{linE}]);
    load(fn);
      
   EvePXhets(:,linE)=(nanmean(EvePeakLocHets));
   EveTXhets(:,linE)=(nanmean(EveValeyLocHets));
   
   EvePXwt(:,linE)=(nanmean(EvePeakLocWT));
   EveTXwt(:,linE)=(nanmean(EveValeyLocWT));
   DeltaXP(:,linE)=EvePXhets(:,linE)-EvePXwt(:,linE);%negative shift is moving anteriorly
   DeltaXT(:,linE)=EveTXhets(:,linE)-EveTXwt(:,linE);
    dt=dtAllLines(:,:,:,:, linE);
for g=1:4
    if g==linE
        continue
    end
    ix=StripeClassNew(linE,:)>0;%==2;% if option ==2 is chosen, only includes stripes wiht a significant shift
    gHetsMean=squeeze(dt(g,:,round(EveTXhets(:,linE))-100,2));
    gWtMean=squeeze(dt(g,:,round(EveTXhets(:,linE))-100,3));
    dIt(linE,g,:)=max((gHetsMean(tmin-9:tmax-9,:)-gWtMean(tmin-9:tmax-9,:)));
    
    gHetsMean=squeeze(dt(g,:,round(EvePXhets(:,linE))-100,2));
    gWtMean=squeeze(dt(g,:,round(EvePXhets(:,linE))-100,3));
    for s=find(ix)%1:7 % for each stripe find the largest |deltaI|
        a=(gHetsMean(tmin-9:tmax-9,s)-gWtMean(tmin-9:tmax-9,s));
%         a=(gHetsMean(tmax-9,s)-gWtMean(tmax-9,s));

        bb=a(a>=0);
        if ~isempty(bb)
            dIp(linE,g,s)=max(bb);
        end
        cc=a(a<=0);
        if ~isempty(cc)
            dIpmin(linE,g,s)=min(cc);
        end
    end

    for s=1:7%
        if StripeClassNew(linE,s)==1
            plot((dIp(linE,g,s)),(DeltaXP(s,linE)/1000),'Marker', 'o', 'Color',...
                'k', 'MarkerEdgeColor','k', 'MarkerSize', MS);
        else
            plot((dIp(linE,g,s)),(DeltaXP(s,linE)/1000),'Marker', 'o', 'Color',...
                'g', 'MarkerEdgeColor','g', 'MarkerSize', MS);
        end
        hold on
        plot([0 0], [-0.03 0.0], ':k')
        hold on

    end

end
ylabel('\DeltaX');
xlabel('max(\DeltaI)|_{t=10:42[min]}');
box off
end
    p1=plot(nan, nan, 'Marker','o','LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    hold on
    p2=plot(nan, nan, 'Marker','o','LineStyle', 'none', 'Color', 'g' , 'MarkerSize',MS);
    p3=plot(nan, nan, 'Marker','o','LineStyle', 'none', 'Color', 'r' , 'MarkerSize',MS);
    legend([p1,p2,p3],{'\DeltaX<=\deltaX_{wt}','\DeltaX>\deltaX_{wt}' ,'perturbed Gap'});
xlim([0.04 0.28])
set(gca, 'FontSize', 8);
axis square

%% S4E
%
cGAP=[ 0.9960,    0.8200    ,0.6040;...
    0.6120    ,0.9540    ,0.8480;...
    0.4480   , 0.5960  ,  0.9000;...
    0.8360    ,0.5640 ,   0.7120];
cGAP=cGAP*0.85;

EvePXhets=nan(7,4);
EveTXhets=nan(6,4);
EvePXwt=nan(7,4);
EveTXwt=nan(6,4);
DeltaXP=nan(7,4);
DeltaXT=nan(6,4);
Mark=['o', 's','d','^'];
Cline=['m','g','b', 'y'];
tmax=42;%42;
tmin=10;%22;
dIp=nan(4,4,7);%(linE,g,s)
dIpmin=nan(4,4,7);%(linE,g,s)
dIt=nan(4,4,6);%(linE,g,s)

left=0.15;
bottom=0.15;
width=6.8;
hight=2;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);

subplot(1,2,2);
for linE=1:4
    fn=fullfile(pn,['PRpeakNvalley4', gapNames{linE}]);
    load(fn);
      
   EvePXhets(:,linE)=(nanmean(EvePeakLocHets));
   EveTXhets(:,linE)=(nanmean(EveValeyLocHets));
   
   EvePXwt(:,linE)=(nanmean(EvePeakLocWT));
   EveTXwt(:,linE)=(nanmean(EveValeyLocWT));
   DeltaXP(:,linE)=EvePXhets(:,linE)-EvePXwt(:,linE);%negative shift is moving anteriorly
   DeltaXT(:,linE)=EveTXhets(:,linE)-EveTXwt(:,linE);
    dt=dtAllLines(:,:,:,:, linE);
for g=1:4
    if g==linE
        continue
    end
    ix=StripeClassNew(linE,:)>0;%==2;% if option ==2 is chosen, only includes stripes wiht a significant shift
    gHetsMean=squeeze(dt(g,:,round(EveTXhets(:,linE))-100,2));
    gWtMean=squeeze(dt(g,:,round(EveTXhets(:,linE))-100,3));
    dIt(linE,g,:)=max((gHetsMean(tmin-9:tmax-9,:)-gWtMean(tmin-9:tmax-9,:)));
    
    gHetsMean=squeeze(dt(g,:,round(EvePXhets(:,linE))-100,2));
    gWtMean=squeeze(dt(g,:,round(EvePXhets(:,linE))-100,3));
    for s=find(ix)%1:7 % for each stripe find the largest |deltaI|
        a=(gHetsMean(tmin-9:tmax-9,s)-gWtMean(tmin-9:tmax-9,s));
%         a=(gHetsMean(tmax-9,s)-gWtMean(tmax-9,s));

        bb=a(a>=0);
        if ~isempty(bb)
            dIp(linE,g,s)=max(bb);
        end
        cc=a(a<=0);
        if ~isempty(cc)
            dIpmin(linE,g,s)=min(cc);
        end
    end

    for s=1:7%
        plot((dIp(linE,g,s)),(DeltaXP(s,linE)/1000),'Marker', Mark(g), 'MarkerFaceColor',...
            cGAP(linE,:), 'MarkerEdgeColor',cGAP(linE,:), 'MarkerSize', 6);
        
        hold on
        plot((dIp(linE,g,s)),(DeltaXP(s,linE)/1000),'Marker', Mark(g), 'MarkerFaceColor',...
            C(s,:), 'MarkerEdgeColor',C(s,:), 'MarkerSize', 3);
        
        plot((dIpmin(linE,g,s)),(DeltaXP(s,linE)/1000),'Marker', Mark(g), 'MarkerFaceColor',...
            cGAP(linE,:), 'MarkerEdgeColor',cGAP(linE,:), 'MarkerSize', 6);
        hold on
        plot((dIpmin(linE,g,s)),(DeltaXP(s,linE)/1000),'Marker', 'o', 'MarkerFaceColor',...
            C(s,:), 'MarkerEdgeColor',C(s,:), 'MarkerSize', 2);
        
        hold on
        plot([0 0], [-0.03 0.0], ':k')
        hold on

    end

     end
ylabel('\DeltaX');
xlabel('max(\DeltaI)|_{t=10:42[min]}');
box off
end
xlim([0.04 0.28])

p1=plot(nan, nan , 'o','Color',cGAP(1,:));% ,'o','LineStyle', 'none', 'MarkerSize',6, 'Color',cGAP(1,:));%,'MarkerFaceColor',cGAP(1,:));
p2=plot(nan, nan, 'o', 'Color',cGAP(2,:));%, 'LineStyle', 'none', 'MarkerSize',6, 'Color',cGAP(2,:));%, 'MarkerFaceColor',cGAP(2,:));
p3=plot(nan, nan, 'o', 'Color',cGAP(3,:));%, 'LineStyle', 'none', 'MarkerSize',6, 'Color',cGAP(3,:));%,'MarkerFaceColor',cGAP(3,:));
p4=plot(nan, nan, 'o', 'Color',cGAP(4,:));%, 'LineStyle', 'none', 'MarkerSize',6, 'Color',cGAP(4,:));%, 'MarkerFaceColor',cGAP(4,:));
p5=plot(nan, nan, Mark(1), 'LineStyle', 'none', 'MarkerSize',3, 'Color','k','MarkerFaceColor','k');
p6=plot(nan, nan, Mark(2), 'LineStyle', 'none', 'MarkerSize',3, 'Color','k', 'MarkerFaceColor','k');
p7=plot(nan, nan, Mark(3), 'LineStyle', 'none', 'MarkerSize',3, 'Color','k','MarkerFaceColor','k');
p8=plot(nan, nan, Mark(4), 'LineStyle', 'none', 'MarkerSize',3, 'Color','k', 'MarkerFaceColor','k');

legend([p1,p2,p3,p4,p5,p6,p7,p8], {'1xhb', '1xKr', '1xgt', '1xkni', 'Hb', 'Kr','Gt', 'Kni'}, 'FontAngle', 'italic');%gapNames);
set(gca, 'FontSize', 8);
axis square

%% plot the example of 1xKr for stripe 4 and stripe 5 separately 
% plot the Intensity comparison of gap activation levels for hets stripes:
    MS=5;% Marker size for plotting
    maxI=1.3;
    Delta4Noise=0.08;
    %     c=[0.7 0 0];
Mark=['o', 's', 'd', '^'];
left=0.15;
bottom=0.15;
width=6.8;
hight=3.3;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
    C=colormap(jet(7));
for s=[3:5]%1:7
    for rgb=1:3
    if (C(s,rgb))~=0
        C(s,rgb)=C(s,rgb)-0.1;
    end
    end
end


% c=[0.7 ,0 ,0.6; 0.13,0.5 ,0.15];

% constract colorcodes for plotting:
mapPolar = [repmat([0,0,1],[32,1]);repmat([1,0,0],[32,1])];
% linear shading from min/max (colormap value) to center (white)
r = repmat(abs(linspace(1,-1,64)),[3,1])';
mapPolar = mapPolar.*r + 1 - r;
  
lineID=2;%using Kr data for example
ExampleStripe=[4,5]; % panel D shows stripe 4, and panel E shows stripe 5
group=2;% these are the hemizygots
gene=4;%Kni is the example gene in Kr hets
t=42;

for ex=1:2
subplot(2,3,2+3*(ex-1))% these are all the hets class 1
    area([0,maxI], [Delta4Noise,maxI+Delta4Noise], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
        hold on
        area([0,maxI], [0,maxI], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
        hold on
        area([0,maxI], ([ -1*Delta4Noise, maxI-Delta4Noise]), 'FaceAlpha', 1,'FaceColor', [1,1,1], 'LineStyle', 'none');
        plot([0 maxI], [0,maxI], '--k');
        plot([0,maxI],[0,0.5*maxI], '--k');
        xlabel(['I@\it{2xKr}']);
        ylabel(['I@\it{1xKr}']);        
        box off
        axis  square
    xlim([0 maxI]);
    ylim([0 maxI]);
    
        s=ExampleStripe(ex);%5
                for g=1:4
                    errorbar(dt42WTm(g,s,lineID),dt42Hetsm(g,s,lineID),...
                        dt42Hetsste(g,s,lineID),dt42Hetsste(g,s,lineID),...
                        dt42WTste(g,s,lineID),dt42WTste(g,s,lineID),...
                        'Marker', Mark(g), 'Color',C(s,:),'LineStyle','none',...
                        'MarkerSize', MS, 'MarkerFaceColor',C(s,:));

                    hold on
                end
        
    
    hold on
%     text(0,1.2, ['1x', gapName4display{lineID}] );
    set(gca, 'FontSize', 8);

end
    p1=errorbar(nan, nan,5, 'Marker',Mark(1),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    p2=errorbar(nan, nan,5, 'Marker',Mark(2),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    p3=errorbar(nan, nan,5, 'Marker',Mark(3),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    p4=errorbar(nan, nan, 'Marker',Mark(4),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
legend([p1,p2,p3,p4],gapNames, 'Position',[0.6092 0.4283 0.1 0.177]);
w=0.18;
ax1=axes('Position',[0.13 0.69 w 0.27]);
fn=fullfile(pn,'KrDiffData200runs5');
load (fn); 
x=(101:900)/1000;
y=10:56;
scalePar=0.7*y(end);

Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];


m=squeeze(nanmean(diffGap(gene,:,:,:,group)));   
imagesc(ax1,x,y,m,[-1,1]);%, 'YDir', 'normal');
                colormap(mapPolar);

hold on
plot(x,y(1)+scalePar*squeeze(nanmean(dtA(gene,:,t,:,3),2)), 'Color', [0 0 0 ], 'LineWidth',1);%'Color',[0.5,0.5,0.5], 'LineWidth',LW);
plot(x,y(1)+scalePar*squeeze(nanmean(dtB(gene,:,t,:,group),2)), 'Color',[0.6  0 0.7], 'LineWidth',1);
ylabel('t[min]');
set(gca,'XTick',[]);
box off
set(gca, 'YDir', 'normal','FontSize',8);%,'LabelFontSize',LabelFontSize );%,'XTick', [100 300 500 700],'XTickLabel', [0.2 0.4 0.6 0.8]); 
                 

fn=fullfile(pn,['PRpeakNvalley4',gapNames{2}]);
load(fn);
axes('Position',[0.13 0.4 w 0.27]);
set(gca, 'FontSize', 8);
IX=idHets;
peakLocs=EvePeakLocHets;
peakI=EvePeakIntensHets/max(mean(Eve(id2x,101:900)));
                   
SSPn=1;
for i=1:length(IX)
    plot((101:900)/1000, smooth(Eve(IX(i),101:900),SSPn)/max(mean(Eve(id2x,101:900))), 'Color',[0.9 0.9 0.9]);
    hold on
end
plot((101:900)/1000, (Eve(IX,101:900))/max(mean(Eve(id2x,101:900))), 'Color',[0.9 0.9 0.9]);
                    hold on
plot(peakLocs/1000, peakI+0., 'Marker','o','MarkerSize', 1.2,...
'LineStyle', 'none', 'Color',[0.7 0.7 0.7] );
plot((101:900)/1000, mean(Eve(IX,101:900))/max(mean(Eve(id2x,101:900))),...
    'Color',Col4Dose(group,:),'LineWidth', 0.6 );
nEm(lineID,group)=length(IX);

box off
ylabel('<I>');
set(gca, 'XTick', []);
ylim([0 2.1])
xlim([0.1 0.9])
plot((101:900)/1000, mean(Eve(id2x,101:900))/max(mean(Eve(id2x,101:900))),...
 'Color',Col4Dose(3,:),'LineWidth', 0.6);
hold on
plot((101:900)/1000, mean(Gap(id2x,101:900))/max(mean(Gap(id2x,101:900))), 'Color',[0.5 0.5 0.5],...
    'LineStyle', '--');
             
text(0.1,1.2 ,[num2str(group-1),'x', gapName4display{lineID}],...
    'FontAngle', 'italic', 'FontSize', 8, 'FontWeight' ,'bold');
text(0.1,1.8 ,['n=',num2str(nEm(lineID,group))],...
    'FontAngle', 'italic', 'FontSize', 8, 'FontWeight' ,'bold');
   
axes('Position',[0.13 0.15 w 0.2]);
xlabel('x/l');
set(gca,'XTick', 0.2:0.2:0.8, 'FontSize', 8);
%         
dltx=nanstd(peakLocs/1000);
   for s=1:7

%        DeltaX1xALL=[DeltaX1xALL;(peakLocs(:,s)-nanmean(EvePeakLocWT(:,s)))/1000];
        errorbar(nanmean(peakLocs(:,s))/1000,...
        nanmean(peakLocs(:,s)-nanmean(EvePeakLocWT(:,s)))/1000,...
        dltx(s),...
        'Marker', 's', 'MarkerFaceColor',C(s,:), 'MarkerEdgeColor',...
    C(s,:),'Color', C(s,:), 'MarkerSize', 4);  
    hold on
      end
hold on
plot([0.1,0.9], [0,0], ':k');
plot([0.1,0.9], [0.01,0.01], '--k');
plot([0.1,0.9], [-0.01,-0.01], '--k');
ylim([-0.04,0.02])
xlim([0.1 0.9])
ylabel('\Deltax');
xlabel('x/l');
box off

annotation('line',[0.232 0.232],...
    [0.150579 0.9621],'Color',[0 1 0],'LineWidth',1,...
    'LineStyle','--');
annotation('line',[0.24428125 0.24428125],[0.151 0.962],'Color',[1 1 0.0667],...
    'LineWidth',1, 'LineStyle','--');


%% Panel 4F (pca for wt and hemizygots stripes overlay)
% get data from individual embryos first:

gapNAMES={'Hb', 'Kr', 'Gt', 'Kni'};
t4checkGap=42;
dt42WT=cell(1,4);
dt42Hets=cell(1,4);
dt42WTThrufs=cell(1,4);
dt42HetsThrufs=cell(1,4);

 for lineID=1:4
     fn=fullfile(pn,[gapNAMES{lineID},'LineWithGenotypeKmeans.mat']);
    load(fn);
    clear dt
    dt(1,:,:)=Hb-min(Hb(:,101:900), [],2);
    dt(2,:,:)=Kr-min(Kr(:,101:900), [],2);
    dt(3,:,:)=Gt-min(Gt(:,101:900), [],2);
    dt(4,:,:)=Kni-min(Kni(:,101:900), [],2);
    %normalize to wt max
    for g=1:4
        dt(g,:,:)=dt(g,:,:)/max(nanmean(dt(g,Genotype==2 & Age>=42-4 & Age<42+4,:)));
    end
     
    % get the gap data per embryo as a seperate matrix for each genotype:
    dt2x=dt(:,Genotype==2 & Age>=42-4 & Age<42+4,:);
    dt1x=dt(:,Genotype==1 & Age>=42-4 & Age<42+4,:); 
   
    % for normalization prior to PCA, get data have unit variance relative
    % to the full wt profile
    for g=1:4
        StdPerGene(lineID, g)=nanstd(dt2x(g,:,101:900),0, 'all');
        MeanPerGene(lineID, g)=nanmean(dt2x(g,:,101:900), 'all');
    end
    fn2=fullfile(pn,['PRpeakNvalley4',gapNAMES{lineID}]);
    load(fn2);
    x4checkWt{lineID}=round(nanmean(EvePeakLocWT));% these are the mean eve peak positions at wt
    x4checkHets=round(nanmean(EvePeakLocHets));% these are the mean eve peak positions at hets

    % compose the data array containing the 4 gap levels at the mean stripe
    % position for *individual* embryos per genotype:
    dt42WT{lineID}=dt2x(:, :,x4checkWt{lineID});
    dt42Hets{lineID}=dt1x(:,:,x4checkHets);
    dt42WTThrufs{lineID}=dt2x(:, :,round(nanmean(EveValeyLocWT)));%+thrufDelta);
    dt42HetsThrufs{lineID}=dt1x(:,:,round(nanmean(EveValeyLocHets)));
 end        
%  run PCA for all stripes only wt PEAKS & plot the pca overlay hets by each line 1&2 pcs (overlay all lines)
% plot 1xpeaks and overlay 2x peaks for comparison
% countEperStripe=[];
clear ProjectedDT ProjectedD
left=0.15;
bottom=0.15;
width=6.8;
hight=3;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
 %
MS=3;
 for lineID=1:4
     sallStripes=[];
    for stripe =1:7
        sallStripes=[sallStripes;squeeze(dt42WT{lineID}(:,:,stripe))'];
        sallStripes=[sallStripes;...
              squeeze(dt42WT{lineID}(:,:,stripe))'];...
     
    end
    standLev1{lineID}=nanstd(sallStripes);
    sallStripes=sallStripes./nanstd(sallStripes);

    meanGAllLines{lineID}=nanmean(sallStripes)';
     sallStripes=sallStripes-nanmean(sallStripes);

[coeff{lineID},score{lineID},latent,tsquared,explained{lineID},mu]=pca(sallStripes);%,'Centered',false);
    for stripe=1:7
        ProjectedData=squeeze(dt42Hets{lineID}(:,:,stripe));
        ProjectedDataT=squeeze(dt42WT{lineID}(:,:,stripe));
         for g=1:4
              ProjectedData(g,:)=ProjectedData(g,:)/standLev1{lineID}(g);
              ProjectedDataT(g,:)=ProjectedDataT(g,:)/standLev1{lineID}(g);
         end
         ProjectedData=ProjectedData-meanGAllLines{lineID};
         ProjectedDataT=ProjectedDataT-meanGAllLines{lineID};

       ProjectedD{lineID}(stripe,:,:)= ProjectedData'*coeff{lineID};   
       ProjectedDT{lineID}(stripe,:,:)= ProjectedDataT'*coeff{lineID};   
    end
end
% now plot:
pc4show=[1,2];
% f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
subplot(1,2,1);
Stripe4Overlay={[1:7];[1:7]; [1: 7]; [1:7]};

for lineID=[1:4]%1:4
    if lineID==2
        sc=-1;
    else
        sc=1;
    end
    for stripe=Stripe4Overlay{lineID}%1:7
         plot(ProjectedD{lineID}(stripe,:,pc4show(1)),...
            sc*ProjectedD{lineID}(stripe,:,pc4show(2)), 'Marker','S', 'LineStyle', 'none',...
        'MarkerEdgeColor','k', 'MarkerFaceColor', C(stripe,:),'MarkerSize', MS-0.7 );%
        hold on
        plot(ProjectedDT{lineID}(stripe,:,pc4show(1)),...
           sc* ProjectedDT{lineID}(stripe,:,pc4show(2)), 'Marker','^', 'LineStyle', 'none',...
        'MarkerEdgeColor','k', 'MarkerFaceColor', 'k','MarkerSize', MS-0.7 );%
        hold on
    end
end
p1=plot(nan, nan,'LineStyle','none','Marker', 's', 'Color','k','MarkerSize',MS );
p2=plot(nan, nan,'LineStyle','none','Marker', '^', 'Color','k' ,'MarkerSize',MS);
legend([p1,p2], {'1xGapPeaks', '2xGapPeaks'});
xlabel(['pc',num2str(pc4show(1))']);
ylabel(['pc',num2str(pc4show(2))']);
axis square
xlim([-2.9 3.3])
ylim([-2.9 3.3])

grid on
box off
set(gca, 'FontSize', 8);

%% Figure S4 
%% plot S4 A-D 
%each genotype Intensity comparison is presented seperately :
% plot the Intensity comparison of gap activation levels for hets stripes:
    MS=4;% Marker size for plotting
    maxI=1.3;
    Delta4Noise=0.1;%0.08;
    C=colormap(jet(7));
for s=[3:5]%1:7
    for rgb=1:3
    if (C(s,rgb))~=0
        C(s,rgb)=C(s,rgb)-0.1;
    end
    end
end
  c=[0.7 ,0 ,0.6; 0.13,0.5 ,0.15];
    %     c=[0.7 0 0];
Mark=['o', 's', 'd', '^'];
left=0.15;
bottom=0.15;
width=6.8;
hight=4;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
for lineID=1:4    
    subplot(2,2,lineID)% these are all the hets class 1
    area([0,maxI], [Delta4Noise,maxI+Delta4Noise], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
        hold on
        area([0,maxI], [0,maxI], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
        hold on
        area([0,maxI], ([ -1*Delta4Noise, maxI-Delta4Noise]), 'FaceAlpha', 1,'FaceColor', [1,1,1], 'LineStyle', 'none');
        plot([0 maxI], [0,maxI], '--k');
        plot([0,maxI],[0,0.5*maxI], '--k');
        xlabel('I@\it{2xGap}');
        ylabel('I@\it{1xGap}');        
        box off
        axis  square
        xlim([0 maxI]);
        ylim([0 maxI]);
    
        for s=1:7
                for g=1:4
                    errorbar(dt42WTm(g,s,lineID),dt42Hetsm(g,s,lineID),...
                        dt42Hetsste(g,s,lineID),dt42Hetsste(g,s,lineID),...
                        dt42WTste(g,s,lineID),dt42WTste(g,s,lineID),...
                        'Marker', Mark(g), 'Color',C(s,:),'MarkerFaceColor',C(s,:),'LineStyle','none',...
                        'MarkerSize', MS);

                    hold on
                end
        end
    
    hold on
    text(0,1.2, ['1x', gapName4display{lineID}] );
    set(gca, 'FontSize', 8);

end
    
    p1=errorbar(nan, nan,1,1,1,1, 'Marker',Mark(1),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    p2=errorbar(nan, nan,1,1,1,1, 'Marker',Mark(2),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    p3=errorbar(nan, nan,1,1,1,1, 'Marker',Mark(3),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
    p4=errorbar(nan, nan,1,1,1,1, 'Marker',Mark(4),'LineStyle', 'none', 'Color', 'k' , 'MarkerSize',MS);
legend([p1,p2,p3,p4],gapNames,'Position',[0.448 0.50 0.0996 0.146]);

set(gca, 'FontSize', 8);

%% Panel E in S4
% % % % %numbers for fig2 & trend over time of "the nuclear response"
for lineID=[2,1,3,4]%4:-1:1
    fn=fullfile(pn, [gapNames{lineID},'DiffData200runs5']);
      load(fn);
      
    for gene=1:4
           for group=[2,3]%1:3% if only hets fig is ploted then run only group 2 &3
                if group~=3
                    m=squeeze(nanmean(diffGap(gene,:,:,:,group)));
                    DeviatPerPos(gene,:,group, lineID)=sqrt(sum(m.^2)/47);
                    Deviat(gene,group, lineID)=sqrt(sum(m.^2,'all')/numel(m));
                else
                     m=squeeze(diffGap(gene,:,:,:,group)).^2; % this includes all bootstraped versions
                     m1=[];
                       for n=1:200 
                            m1=[m1;squeeze(m(n,:,:))];
                       end
                    DeviatPerPos(gene,:,group, lineID)=squeeze(sqrt((nanmean(m1))));
                    Deviat(gene,group, lineID)=sqrt(nanmean(m,'all'));
                    m2=[];
                        for n=1:200 
                                m2=[m2,squeeze(m(n,:,:))];
                        end

                    DeviatAlongT(gene,:,group, lineID)=sqrt(nanmean(m2,2));
                end
            end
         end
end

% plot trend over egg length
cGAP=[ 0.9960,    0.8200    ,0.6040;...
    0.6120    ,0.9540    ,0.8480;...
    0.4480   , 0.5960  ,  0.9000;...
    0.8360    ,0.5640 ,   0.7120];

cEve=colormap(jet(7));
gapNamesforshow={'hb', 'Kr', 'gt', 'kni'};

Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];

left=0.05;
bottom=0.05;
width=6.850394;
hight=2;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
for lineID=1:4
for group=2
%     subplot(2,4,lineID+(4*(group-1)));
    subplot(1,4,lineID);
    
    for gene=1:4
            plot((101:900)/1000, DeviatPerPos(gene,:,group, lineID), 'Color', cGAP(gene,:), 'LineWidth',0.7);
            hold on
        ylim([0 0.5])
        xlim([0.1 0.9])
        box off
        xlabel('x/l');
        ylabel('(\Sigma_{t}{\DeltaI(x,t)}^2/n)^{1/2}')
    %     end
    text(0.6,0.45, [num2str(group-1),'x',gapNamesforshow{lineID}], 'FontAngle', 'italic');

    end
    for g=1:4
        plot((101:900)/1000, DeviatPerPos(g,:,3, lineID), 'Color', cGAP(g,:), 'LineStyle', '--', 'LineWidth',0.7);
        hold on
        
    end
    
end
end
legend(gapNames,'Position',[0.000622 0.00605 0.0988 0.314]);

%% Panel S4G (pca for hemizygots GAP at EVE PEAKS & TROUGHS)
clear ProjectedDT ProjectedD ProjectedDataT ProjectedData
Stripe4Overlay={[1:7];[1:7]; [1:7]; [1:7]};

left=0.15; 
bottom=0.15;
width=6.8;
hight=3;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
for lineID=1:4
    for stripe=1:7
        ProjectedData=squeeze(dt42Hets{lineID}(:,:,stripe));
        if stripe~=7
        ProjectedDataT=squeeze(dt42HetsThrufs{lineID}(:,:,stripe));
        end
         for g=1:4
              ProjectedData(g,:)=ProjectedData(g,:)/standLev1{lineID}(g);
              if stripe~=7
              ProjectedDataT(g,:)=ProjectedDataT(g,:)/standLev1{lineID}(g);
              end
         end
         ProjectedData=ProjectedData-meanGAllLines{lineID};
         if stripe~=7
            ProjectedDataT=ProjectedDataT-meanGAllLines{lineID};
         end

       ProjectedD{lineID}(stripe,:,:)= ProjectedData'*coeff{lineID};   
       if stripe~=7
        ProjectedDT{lineID}(stripe,:,:)= ProjectedDataT'*coeff{lineID};   
       end
    end
end
% now plot:
pc4show=[1,2];
% subplot(1,2,1);

for lineID=1:4
    if lineID==2
        sc=-1;
    else 
        sc=1;
    end
        
    for stripe=Stripe4Overlay{lineID}%1:7
    
        plot(ProjectedD{lineID}(stripe,:,pc4show(1)),...
            sc*ProjectedD{lineID}(stripe,:,pc4show(2)), 'Marker','S', 'LineStyle', 'none',...
        'MarkerEdgeColor','k', 'MarkerFaceColor', C(stripe,:),'MarkerSize', MS-0.7 );%
        hold on
        if stripe~=7
        plot(ProjectedDT{lineID}(stripe,:,pc4show(1)),...
            sc*ProjectedDT{lineID}(stripe,:,pc4show(2)), 'Marker','^', 'LineStyle', 'none',...
        'MarkerEdgeColor','k', 'MarkerFaceColor', 'k','MarkerSize', MS-0.7 );%
        hold on
        end
    end
end
p1=plot(nan, nan,'LineStyle','none','Marker', 's', 'Color','k','MarkerSize',MS );
p2=plot(nan, nan,'LineStyle','none','Marker', '^', 'Color','k' ,'MarkerSize',MS);
legend([p1,p2], {'1xGapPeaks', '1xGapTroughs'}, 'Location', 'southeast');
xlabel(['pc',num2str(pc4show(1))']);
ylabel(['pc',num2str(pc4show(2))']);
axis square
xlim([-2.9 3.3])
ylim([-2.9 3.3])

grid on
box off
set(gca, 'FontSize', 8);

