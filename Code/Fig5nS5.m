%% Plot figure 5 and S5  (nulls stripes PCA and responsiveness)
%% get data arrays from individual embryos for pca:
clearvars;
if ~exist('processedData', 'dir')
    mkdir('processedData')
end
fol=cd;
pn = fullfile(fol,'processedData');

gapNAMES={'Hb', 'Kr', 'Gt', 'Kni'};
t4checkGap=42;
dt42WT=cell(1,4);%each cell contains one line data in the form: 4 gene, n embryos, 7stripes
dt42Hets=cell(1,4);%same as dt42WT
dt42NullByPhysX=cell(1,4);% 4 genes, n embryos, 7 possible stripes according to physical order, not identity

 for lineID=1:4
     fn=fullfile(pn,[gapNAMES{lineID},'LineWithGenotypeKmeans.mat']);
    load(fn);

    clear dt
    dt(1,:,:)=Hb-min(Hb(:,101:900), [],2);
    dt(2,:,:)=Kr-min(Kr(:,101:900), [],2);
    dt(3,:,:)=Gt-min(Gt(:,101:900), [],2);
    dt(4,:,:)=Kni-min(Kni(:,101:900), [],2);

    %normalize to wt max like they do in MP & JD work:
    for g=1:4
        dt(g,:,:)=dt(g,:,:)/max(nanmean(dt(g,Genotype==2 & Age>=42-4 & Age<42+4,:)));
    end
    % get the gap data per embryo as a seperate matrix for each genotype:
    dt2x=dt(:,Genotype==2 & Age>=42-4 & Age<42+4,:);
    dt1x=dt(:,Genotype==1 & Age>=42-4 & Age<42+4,:); 
    dt0x=dt(:,Genotype==0 & Age>=42-4 & Age<42+4,:);
   
%     get eve positions:  
    fn1=fullfile(pn,['PRpeakNvalley4',gapNAMES{lineID}]);
    load(fn1);% these are stripe positions in individual embryos, with an 8 min mean delay
    x4checkWt{lineID}=round(nanmean(EvePeakLocWT));% these are the mean eve peak positions at wt
    x4checkHets{lineID}=round(nanmean(EvePeakLocHets));% these are the mean eve peak positions at hets
    x4checkNullsPhys{lineID}=round(nanmean(EvePeakLocNulls));% this is the mean eve peak positions at hets
    if lineID==2
        x4checkNullsPhys{lineID}(4)=[];% exclude the very rare stripe in Kr nulls
    end
    % compose the data array containing the 4 gap levels at the mean stripe
    % position for individual embryos per genotype:
    dt42WT{lineID}=dt2x(:, :,x4checkWt{lineID});
    dt42Hets{lineID}=dt1x(:,:,x4checkHets{lineID});
    idx=x4checkNullsPhys{lineID};
    idx=idx(~isnan(idx));
    dt42NullByPhysX{lineID}(:,:,1:length(idx))=dt0x(:,:,idx);
 end
%% plot eve stripes insets for each panel
GAPNames={'Hb','Kr', 'Gt', 'Kni'};
gapNamesforshow={'hb','Kr', 'gt', 'kni'};


left=0.15;
bottom=0.15;
width=6.850394;
hight=1.8;
Tag={'A','B','C', 'D', 'E', 'F','G','H','I','J','K', 'L'};
roman={'I', 'II', 'III', 'IV','V','VI', 'VII'};
Place=[4 1 3 2 ];
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);

Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];

C=colormap(jet(7));
for s=[3:5]%1:7
    for rgb=1:3
    if (C(s,rgb))~=0
        C(s,rgb)=C(s,rgb)-0.1;
    end
    end
end

for lineID=1:4
    fn=fullfile(pn,['PRpeakNvalley4',GAPNames{lineID}]);
    load(fn);

    subplot(2,2,Place(lineID))
    plot((101:900)/1000, mean(Gap(id2x,101:900))/max(mean(Gap(id2x,101:900))),...
        'Color',[0.2 0.2 0.2],'LineWidth', 0.6, 'LineStyle','--' );
    hold on
    IX=idNulls;
    peakLocs=EvePeakLocNulls;
    peakI=EvePeakIntensNulls/max(mean(Eve(id2x,101:900)));
    if lineID==2
        peakLocs(:,4)=[];
        peakI(:,4)=[];
        
    end
                for i=1:length(IX)
                    plot((101:900)/1000, smooth(Eve(IX(i),101:900),20)/max(mean(Eve(id2x,101:900))), 'Color',[0.9 0.9 0.9]);
                    hold on
                end
                    plot(peakLocs/1000, peakI+0., 'Marker','o','MarkerSize', 1.2,...
                 'LineStyle', 'none', 'Color',[0.7 0.7 0.7] );
             
             plot((101:900)/1000, mean(Eve(id2x,101:900))/max(mean(Eve(id2x,101:900))),...
                    'Color',Col4Dose(3,:),'LineWidth', 0.6 );
                hold on
                    plot((101:900)/1000, mean(Eve(IX,101:900))/max(mean(Eve(id2x,101:900))),...
                    'Color',Col4Dose(1,:),'LineWidth', 0.6 );
             for k=1:size(peakLocs,2)
                text(nanmean(peakLocs(:,k))/1000-0.05, nanmean(peakI(:,k))+0.3,roman{k},...
                    'FontSize', 7, 'FontAngle','italic');
                hold on
                
             end   
                    xlabel('x/l');
                    ylabel('I_{eve}');
                    box off
                    xlim([0.1 0.9])
                    ylim([0 2]);
            text(-0.1,2.3,Tag{Place(lineID)}, 'FontWeight','bold', 'FontSize', 10);
            set(gca, 'FontSize', 8);

end

%% now plot the pca for wt in each line with null stripes projected on top:
% NULL stripes are projected onto the wt pc plane and presented in black
% markers according to their serial rank along the egg (serial rank annotatino needs manual adjustment for its placing):

left=0.15;
bottom=0.15;
width=6.850394;
hight=6;
Tag={'A','B','C', 'D', 'E', 'F','G','H','I','J','K', 'L'};
roman={'I', 'II', 'III', 'IV','V','VI', 'VII'};
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
Mark4Serial=['x', '.', 'd', '*', '^', '+', 's'];
num2rom={'I', 'II', 'III', 'IV','V','VI', 'VII'};
MS=5;
p=cell(1,7);
pc4show=[1,2];
Place=[4 1 3 2 ];

for lineID=1:4
    subplot(2,2,Place(lineID))
    sallStripes=[];
    for stripe =1:7
        sall=squeeze(dt42WT{lineID}(:,:,stripe));
        sallStripes=[sallStripes;sall'];
    end
    TotEmNum=size(dt42WT{lineID}(:,:,stripe),2);

    standLev{lineID}=nanstd(sallStripes);
    sallStripes=sallStripes./nanstd(sallStripes);

    meanGPerLine(lineID,:)=nanmean(sallStripes);
    sallStripes=sallStripes-nanmean(sallStripes);
    [coeff,score,latent,tsquared,explained,mu]=pca(sallStripes);
% [coeff,score,latent,tsquared,explained,mu]=pca(sallStripes,'Centered',false);

    for stripe=1:7
        ix=1+(stripe-1)*TotEmNum:stripe*TotEmNum; 

        plot(score(ix,pc4show(1)), score(ix,pc4show(2)), 'Marker','o', 'LineStyle', 'none',...
            'Color',C(stripe,:), 'MarkerSize', MS );
        hold on

    end
       
for stripe =1:size(dt42NullByPhysX{lineID},3)
        NullData=squeeze(dt42NullByPhysX{lineID}(:,:,stripe));

         for g=1:4
              NullData(g,:)=NullData(g,:)/standLev{lineID}(g);
                NullData(g,:)=NullData(g,:)-meanGPerLine(lineID,g);
         end
           ProjectedNull= NullData'*coeff;
if lineID==3
    p{stripe}=plot(ProjectedNull(:,pc4show(1)), ProjectedNull(:,pc4show(2)),...
        'Marker',Mark4Serial(stripe), 'LineStyle', 'none',...
        'Color','k', 'MarkerSize', MS-2);
else
    plot(ProjectedNull(:,pc4show(1)), ProjectedNull(:,pc4show(2)),...
        'Marker',Mark4Serial(stripe), 'LineStyle', 'none',...
        'Color','k', 'MarkerSize', MS-2);
end    

text(nanmean(ProjectedNull(:,pc4show(1)))-0.9,nanmean( ProjectedNull(:,pc4show(2)))-0.9,num2rom{stripe},...
        'Color','k', 'FontAngle', 'italic', 'FontWeight','bold', 'FontSize',7 );

end
xlabel(['pc',num2str(pc4show(1)), '(', num2str(explained(pc4show(1)),3), '%)']);
ylabel(['pc',num2str(pc4show(2)), '(', num2str(explained(pc4show(2)),3), '%)']);
axis square
xlim([-3.5 3.5])
ylim([-3.5 3.5])

set(gca, 'FontSize', 8);
grid on
text(1.5, 2.5,['0x',gapName4display{lineID}], 'FontAngle', 'italic', 'FontWeight','bold');

if lineID==3
    legend([p{1}, p{2}, p{3}, p{4},p{5},p{6}, p{7}], num2rom,'FontSize', 6,...
        'Position',[0.46377 0.439 0.0928 0.188]);
end
end
% % % % %%
% % % % f=openfig ('Fig5Main');
% % % % exportgraphics(f,'Fig5pca.jpg', 'Resolution',300)
%% Figure S5 (Former Panel A is now part of Figure 6,chi square analysis)
%% Plot Panel S5B

% numbers for fig2 & trend over time of "the nuclear response"
for lineID=[2,1,3,4]%4:-1:1
    fn=fullfile(pn,[gapNames{lineID},'DiffData200runs5']);
      load(fn);
      
    for gene=1:4
           for group=[1,3]%1:3
                if group~=3
                    % dtA and dtB dim: 4(genes),n(randomSplits), 800(ap bins), 47(t),3genotypes
                    m=squeeze(nanmean(diffGap(gene,:,:,:,group)));
                    DeviatPerPos(gene,:,group, lineID)=sqrt(sum(m.^2)/47);
                    Deviat(gene,group, lineID)=sqrt(sum(m.^2,'all')/numel(m));
%                     DeviatAlongT(gene,:,group, lineID)=sqrt(sum(m.^2,2)/800);
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
%     subplot(1,4,lineID);
for group=1%:2
%     subplot(2,4,lineID+(4*(group-1)));
    subplot(1,4,lineID);
    for gene=1:4
            plot((101:900)/1000, DeviatPerPos(gene,:,group, lineID), 'Color', cGAP(gene,:), 'LineWidth',0.7);
            hold on
        ylim([0 0.9])
        xlim([0.1 0.9])
        box off
        xlabel('x/l');
        ylabel('(\Sigma_{t}{\DeltaI(x,t)}^2/n)^{1/2}')
    %     end
    text(0.6,0.8, [num2str(group-1),'x',gapNamesforshow{lineID}], 'FontAngle', 'italic');

    end
    for g=1:4
        plot((101:900)/1000, DeviatPerPos(g,:,3, lineID), 'Color', cGAP(g,:), 'LineStyle', '--', 'LineWidth',0.7);
        hold on
        
    end
                      
end
end
legend(gapNames,'Position',[0.000622 0.00605 0.0988 0.314]);
%% Plot Panel S5C
% plot the example of 0xKr 
% Generate the mean data array for comarison
w=0.2;
dt42WTm=nan(4,7);%mean expression level at: 4 gene, 7stripes
dt42NullByPhysXm=nan(4,7);%  7 possible stripes according to physical order, not identity
sID=[1 2 7 6 7];% these are the Kr nulls stripe identities estimated from the PCA

% for lineID=1:4
lineID=2;
fn=fullfile(pn,['PRpeakNvalley4',gapNAMES{lineID}]);
load(fn);% these are stripe positions in individual embryos, with an 8 min mean delay
EvePeakLocNulls(:,4)=nan;% exclude the very rare stripe that do not show up in the mean eve profile
EvePeakIntensNulls(:,4)=nan;        
x4checkWt=round(nanmean(EvePeakLocWT))-100;% these are the mean eve peak positions at wt
x4checkNullsPhys=round(nanmean(EvePeakLocNulls))-100;% this is the index for eve peak positions at nulls
x4checkNullsPhys=x4checkNullsPhys(~isnan(x4checkNullsPhys));
%     x4checkNullsPhys(4)=[];% exclude the very rare stripe in Kr nulls

% compose the data array containing the 4 gap levels at the mean stripe
% position for individual embryos per genotype:
fn2=fullfile(pn,[gapNAMES{lineID},'DiffData200runs5']);
load(fn2);
dt42WTm=nanmean(squeeze(dtA(:,:,t4checkGap-9,x4checkWt,3)),2);
dt42NullByPhysXm(:,1:length(x4checkNullsPhys))=...
    nanmean(squeeze(dtA(:,:, t4checkGap-9,x4checkNullsPhys,1)), 2);
% and standard errors:
dt42WTste=nanstd(squeeze(dtA(:,:,t4checkGap-9,x4checkWt,3)),[],2);
dt42NullByPhysXste(:,1:length(x4checkNullsPhys))=...
nanstd(squeeze(dtA(:,:,t4checkGap-9,x4checkNullsPhys,1)), [],2);
% end
 % constract colorcodes for plotting:
mapPolar = [repmat([0,0,1],[32,1]);repmat([1,0,0],[32,1])];
% linear shading from min/max (colormap value) to center (white)
r = repmat(abs(linspace(1,-1,64)),[3,1])';
mapPolar = mapPolar.*r + 1 - r;
 
MS=4;% Marker size for plotting
maxI=1.3;
Delta4Noise=0.1;%0.08;

Mark=['o', 's', 'd', '^'];
left=0.15;
bottom=0.15;
width=6.8;
hight=5;
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

subplot(4,2,[4,6])% 

area([0,maxI], [Delta4Noise,maxI+Delta4Noise], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
        hold on
        area([0,maxI], [0,maxI], 'FaceAlpha', 0.1,'FaceColor', [0.7,0.7,0.7], 'LineStyle', 'none');
        hold on
        area([0,maxI], ([ -1*Delta4Noise, maxI-Delta4Noise]), 'FaceAlpha', 1,'FaceColor', [1,1,1], 'LineStyle', 'none');
        plot([0 maxI], [0,maxI], '--k');
        plot([0,maxI],[0,0.5*maxI], '--k');
        xlabel(['I@\it{2xKr}']);
        ylabel(['I@\it{0xKr}']);        
        box off
        axis  square
    xlim([0 maxI]);
    ylim([0 maxI]);

 for s=1:5
        for g=1:4
            errorbar(dt42WTm(g,sID(s)),dt42NullByPhysXm(g,s),...
            dt42NullByPhysXste(g,s),dt42NullByPhysXste(g,s),...
            dt42WTste(g,sID(s)),dt42WTste(g,sID(s)),...
            'Marker', Mark(g), 'Color',C(sID(s),:), 'MarkerFaceColor',C(sID(s),:),'LineStyle','none',...
            'MarkerSize', MS);
        
                hold on
        end
        if s==5
            for g=1:4
            errorbar(dt42WTm(g,sID(s)),dt42NullByPhysXm(g,s),...
            dt42NullByPhysXste(g,s),dt42NullByPhysXste(g,s),...
            dt42WTste(g,sID(s)),dt42WTste(g,sID(s)),...
            'Marker', Mark(g), 'Color','k', 'MarkerFaceColor',C(sID(s),:),...
            'MarkerEdgeColor', 'k', 'LineStyle','none',...
            'MarkerSize', MS);
        hold on
            end
        end
 end
    
                 
group=1;
w=0.17;

axes('Position',[0.13 0.8 w 0.124]);
set(gca, 'FontSize', 8);
IX=idNulls;
peakLocs=EvePeakLocNulls;
peakI=EvePeakIntensNulls/max(mean(Eve(id2x,101:900)));

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
%         POSBot=[Left(group), bottoM(1)+2*panelH*(lineID-1)+spacer*(lineID-1), panelW,panelH];
   
axes('Position',[0.13 0.67 w 0.124]);

xlabel('x/l');
set(gca,'XTick', 0.2:0.2:0.8, 'FontSize', 8);
%         
peakLocs(:,4)=[];
dltx=nanstd(peakLocs/1000);
for s=1:length(sID)
        errorbar(nanmean(peakLocs(:,s))/1000,...
        nanmean(peakLocs(:,s)-nanmean(EvePeakLocWT(:,sID(s))))/1000,...
        dltx(s),...
        'Marker', 's', 'MarkerFaceColor',C(sID(s),:), 'MarkerEdgeColor',...
    C(sID(s),:),'Color', C(sID(s),:), 'MarkerSize', 4);  
    hold on
    if s==5
         errorbar(nanmean(peakLocs(:,s))/1000,...
        nanmean(peakLocs(:,s)-nanmean(EvePeakLocWT(:,sID(s))))/1000,...
        dltx(s),...
        'Marker', 's', 'MarkerFaceColor',C(sID(s),:), 'MarkerEdgeColor',...
    'k','Color','k', 'MarkerSize', 4);  
    end
end
hold on
plot([0.1,0.9], [0,0], ':k');
plot([0.1,0.9], [0.01,0.01], '--k');
plot([0.1,0.9], [-0.01,-0.01], '--k');
ylim([-0.3,0.26])
xlim([0.1 0.9])
ylabel('\Deltax');
xlabel('x/l');
box off

t=42;
x=(101:900)/1000;
y=10:56;
scalePar=0.7*y(end);

g=[4,3,1];
for gene=1:length(g)
    axes('Position',[0.13 0.45-(gene-1)*0.124 w 0.124]);

    m=squeeze(nanmean(diffGap(g(gene),:,:,:,group)));   
    imagesc(x,y,m,[-1,1]);%, 'YDir', 'normal');
                    colormap(mapPolar)
  
    hold on
    plot(x,y(1)+scalePar*squeeze(nanmean(dtA(g(gene),:,t,:,3),2)), 'Color', [0 0 0 ], 'LineWidth',1);%'Color',[0.5,0.5,0.5], 'LineWidth',LW);
    plot(x,y(1)+scalePar*squeeze(nanmean(dtB(g(gene),:,t,:,group),2)), 'Color',[0.6  0 0.7], 'LineWidth',1);
    ylabel('t[min]');
    set(gca,'XTick',[]);
    box off
    set(gca, 'YDir', 'normal','FontSize',8);%,'LabelFontSize',LabelFontSize );%,'XTick', [100 300 500 700],'XTickLabel', [0.2 0.4 0.6 0.8]); 
end
    xlabel('x/l');
    set(gca,'XTick',0.2:0.2:0.8);

    xp=[nan 0.2029,0.216,0.24172,0.2678];
    
for s=2:length(sID)
    annotation('line',[xp(s) xp(s)],...
        [0.20 0.925],'Color', C(sID(s),:),'LineWidth',1,...
    'LineStyle','--');

end

