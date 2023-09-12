
% This file calls other files:
% * sortGAP2dosek2 - sort embryos to genotype according to manipulated gene protein level
% * RcnstrctGapTseriesWithStd4RandSplit - bootstrap data and reconstruct gap dynamics 
% * RcnstrctGapTseriesWithStd - reconstruct gap dynamics
% * compute_SigmaX - coputes positional error for gene expression along the egg (from JD et. al 2013)
%% Make sure to add the following folders (with subfolders) to your path: "GapData", "EveData", "Code"
%% load gap data and sort to genotype:
% gap data in the files "gap+gap-line4gap.mat" includes:
% protein levels (raw fluorecent intensity of dorsal profile): 'Hb', 'Kr', 'Gt', 'Kni',
% temporal stamp of an embryo (in minutes into nuclear cycle 14): 'Age'
% embryo index in the raw data (the sequential acquisition index of embryos within a given hemizygot line), 'ix');

clearvars
fol=cd;
pn = fullfile(fol,'processedData');

lineNames={'hb+hb-', 'Kr+Kr-', 'gt+gt-', 'kni+kni-'};
GapNames={'Hb', 'Kr', 'Gt', 'Kni'};
    
    winSpan=5;% time window for genotype asignement, in minutes
    progress=1; % in minutes, when smaller than winSpan, there are ovelaps in embryo time groups
    StartT=6; % seems like there are not enough embryos for reasonable clustering prior to this t for hb line.
    endT=60;

    ndxMaxEdges=[200, 450;...
                500,550;...
                200,400;...
                600,650]; % this is the egg section used to identify nulls
   Thresh4NullsAll=[460,55,90,245]; % this is the threshold for nulls per line at the specified egg section
 for LineID=1:4% the index of lines as listed in lineNames below
  
    load([lineNames{LineID},'line4gap']);
    GapName=GapNames{LineID};

    Hb=(Hb-min(Hb(:,101:900),[],2));
    Kr=(Kr-min(Kr(:,101:900),[],2));
    Gt=(Gt-min(Gt(:,101:900),[],2));
    Kni=(Kni-min(Kni(:,101:900),[],2));

    if LineID==1
        Gap4sort=Hb;
        Kni=Kni-0.1*Hb-0.065*Gt; % remove the leakage of Hb and Gt into Kni channel (due to cross talk between the 2 reds for Gt and due to joint scann for heocst and Kni channels (which can also activate Hb blue))
    elseif LineID==2
        Gap4sort=Kr;
         Kni=Kni-0.13*Hb-0.07*Gt;
    elseif LineID==3
        Gap4sort=Gt;
        Kni=Kni-0.1*Hb-0.12*Gt;

    else
        Kni=Kni-0.06*Hb-0.16*Gt;
        Gap4sort=Kni;
    end
    
    ndxMax=ndxMaxEdges(LineID,1):ndxMaxEdges(LineID,2);

    % sort embryos with Threshold and kmeans for 1x vs. 2x 
 
    Thresh4Nulls=Thresh4NullsAll(LineID);
    k=2;
    genotype=NaN(size(Age,2),(endT-winSpan));
    genStatClass2=NaN(size(Age,2),3);% this will keep the record of how many embryos are at any genotype for each classification cycle 

    for i=StartT:progress:(endT-winSpan)
%         i
        idx = find(Age>=i & Age<=i+winSpan);
        idNulls=find(max(Gap4sort(idx,ndxMax),[],2)<Thresh4Nulls);
        idx1=idx; idx1(idNulls)=[];
        if LineID==4
             G4sort=Gap4sort(idx1,25:970);%For Kni 
        else
            G4sort=Gap4sort(idx1,101:900);%for all but Kni;
        end

        [idHets, id2x]=sortGAP2dosek2(G4sort,k);
        genotype(idx1(id2x),i)=2;
        genotype(idx1(idHets),i)=1;
        genotype(idx(idNulls),i)=0;


        genStatClass2(i,1)=length(idNulls);
        genStatClass2(i,2)=length(idHets);
        genStatClass2(i,3)=length(id2x);
%%%%%%%%visualize genotype assignment per time frame:
% % %         figure;
% % %         plot(Gap4sort(idx(idNulls),:)', 'g'); hold on;  
% % %         plot(Gap4sort(idx1(id2x),:)', 'b');ylim([0,4000]);
% % %         plot(Gap4sort(idx1(idHets),:)', 'k');
% % % 
% % %         title(['t=' , num2str(i), 'min'])
    end
    
% % % % % visualize how many embryos were asigned per genotype at any time point   
% %     figure;
% %     plot ( genStatClass2)
% %     legend({'0x','1x','2x'})
% %     xlabel('t')
% %     ylabel('#')
% %     title([GapName,' line , classification by thresh for 0x & kmeans for 1x/2x']);
% %     axis tight

    % Get the most likely genotype for null threshold and kmeans:
    Genotype=NaN(size(Age));
     notClass2Gen=[];
     g2=nan(size(Genotype));
     for i=1:length(Genotype)
        g=genotype(i,~isnan(genotype(i,:)));
        if isempty(g)
           notClass2Gen=[notClass2Gen, i];% these embryos are outside the time window for genotype assignent
        end
        g2(i)=mean(g);
        Genotype(i)=round(mean(g));
     end

     dg=g2-Genotype;
     ixHard2Class=find(dg>0.3 & dg<0.7) 
     Genotype(ixHard2Class)=NaN;% excluded


    g2(ixHard2Class)
    % get numbers per line & genotype
    display([GapName ,' n2x=', num2str(sum(Genotype==2))]);
    display([GapName ,' n1x=', num2str(sum(Genotype==1))]);
    display([GapName ,' n0x=', num2str(sum(Genotype==0))]);
    display([GapName ,' nAll=', num2str(sum(~isnan(Genotype)))]);

    fn=fullfile(pn,[GapName, 'LineWithGenotypeKmeans']);
    save(fn, 'Hb', 'Kr','Gt','Kni', 'Genotype', 'ix', 'Age',...
        'StartT','progress','endT', 'winSpan');
 end

%% Bootstrap gap data:
GAPNamesProtein={'Hb','Kr', 'Gt', 'Kni'};
n=200;

for lineID=1:4
    fn=fullfile(pn,[GAPNamesProtein{lineID}, 'LineWithGenotypeKmeans']);
    load (fn);

            %subtruct background from individual embryos
            Hb=(Hb-min(Hb(:,101:900),[],2));
            Kr=(Kr-min(Kr(:,101:900),[],2));
            Gt=(Gt-min(Gt(:,101:900),[],2));
            Kni=(Kni-min(Kni(:,101:900),[],2));
            
        StartT=6;
        winSpan=8;
        progress=1;
        ageSpan=StartT:progress:(60-winSpan);

       [MeanHbFine,MeanKrFine,MeanGtFine,MeanKniFine,stdHbFine,stdKrFine,...
            stdGtFineB,stdKniFineB]=RcnstrctGapTseriesWithStd(Hb,Kr,...
            Gt,Kni,Age,Genotype,...
            progress,winSpan,StartT);
% %         make sure  dt runs between 

        MeanHbFine=MeanHbFine/max(max(MeanHbFine(:,:,3)));
        MeanKrFine=MeanKrFine/max(max(MeanKrFine(:,:,3)));
        MeanGtFine=MeanGtFine/max(max(MeanGtFine(:,:,3)));
        MeanKniFine=MeanKniFine/max(max(MeanKniFine(:,:,3)));

            dt(1,:,:,:)=MeanHbFine;
            dt(2,:,:,:)=MeanKrFine;
            dt(3,:,:,:)=MeanGtFine;
            dt(4,:,:,:)=MeanKniFine;
        
    diffGap=nan(4,n,length(ageSpan),800,3);

    dtA=nan(4,n,length(ageSpan),800,3);
    dtB=nan(4,n,length(ageSpan),800,3);
    for i=1:n % run multiple random split of each dataset and reconstruct 2 time series & their I_diff
    n=i
    line=lineID
    idxA=[];
    idxB=[];
        for idg=1:3
     %   compose a vector of indexes such that each genotype is samples wiht replacements seperately 
            idx=find(Genotype==idg-1);%find(Genotype==idg);
            emN=length(idx);
            for resamp=1:emN
                   rs=randperm(emN,1);
                   idxA=[idxA,idx(rs)] ;
                   rs=randperm(emN,1);
                   idxB=[idxB,idx(rs)];
            end
               
        end
    [MeanHbFineA,MeanKrFineA,MeanGtFineA,MeanKniFineA,stdHbFineA,stdKrFineA,...
        stdGtFineA,stdKniFineA]=RcnstrctGapTseriesWithStd4RandSplit(Hb(idxA,:),Kr(idxA,:),...
        Gt(idxA,:),Kni(idxA,:),Age(idxA),Genotype(idxA),...
        progress,winSpan,StartT);

    [MeanHbFineB,MeanKrFineB,MeanGtFineB,MeanKniFineB,stdHbFineB,stdKrFineB,...
        stdGtFineB,stdKniFineB]=RcnstrctGapTseriesWithStd4RandSplit(Hb(idxB,:),Kr(idxB,:),...
        Gt(idxB,:),Kni(idxB,:),Age(idxB),Genotype(idxB),...
        progress,winSpan,StartT);

    % % % normalize each time series to its mean and max,such that each runs between 0-1 for the wt :
    % % 
%     MeanHbFineA=MeanHbFineA-min(min(MeanHbFineA(:,:,3))); 
    MeanHbFineA=MeanHbFineA/max(max(MeanHbFineA(:,:,3)));
%     MeanKrFineA=MeanKrFineA-min(min(MeanKrFineA(:,:,3))); 
    MeanKrFineA=MeanKrFineA/max(max(MeanKrFineA(:,:,3)));
%     MeanGtFineA=MeanGtFineA-min(min(MeanGtFineA(:,:,3))); 
    MeanGtFineA=MeanGtFineA/max(max(MeanGtFineA(:,:,3)));
%     MeanKniFineA=MeanKniFineA-min(min(MeanKniFineA(:,:,3))); 
    MeanKniFineA=MeanKniFineA/max(max(MeanKniFineA(:,:,3)));

    MeanHbFineB=MeanHbFineB/max(max(MeanHbFineB(:,:,3)));
    MeanKrFineB=MeanKrFineB/max(max(MeanKrFineB(:,:,3)));
    MeanGtFineB=MeanGtFineB/max(max(MeanGtFineB(:,:,3)));
    MeanKniFineB=MeanKniFineB/max(max(MeanKniFineB(:,:,3)));

    % dtA and dtB dim: 4(genes),n(randomSplits), 800(ap bins), 47(t),3genotypes
    dtA(1,i,:,:,:)=MeanHbFineA;
    dtA(2,i,:,:,:)=MeanKrFineA;
    dtA(3,i,:,:,:)=MeanGtFineA;
    dtA(4,i,:,:,:)=MeanKniFineA;

    dtB(1,i,:,:,:)=MeanHbFineB;
    dtB(2,i,:,:,:)=MeanKrFineB;
    dtB(3,i,:,:,:)=MeanGtFineB;
    dtB(4,i,:,:,:)=MeanKniFineB;

        % dtA and dtB dim: 4(genes),n(randomSplits), 800(ap bins), 47(t),3genotypes
        for gene=1:4
            for group =1:3
              % subtract wt from mutants (or wt in the ctrl case)
                diffGap(gene,i,:,:,group)=dtB(gene,i,:,:,group)-dtA(gene,i,:,:,3);
            end
        end
    end
  fn=fullfile(pn,[GAPNamesProtein{lineID},'DiffData',num2str(n),'runs5']);
  save(fn,'diffGap','dtA','dtB', 'dt'); 
end

%% F2 Panels A-E:
clearvars
fol=cd;
pn=fullfile(fol, 'processedData');
GAPNamesProtein={'Hb','Kr', 'Gt', 'Kni'};
gapNamesGenotype={'hb','Kr', 'gt', 'kni'};

% constract colorcodes for plotting:
mapPolar = [repmat([0,0,1],[32,1]);repmat([1,0,0],[32,1])];
% linear shading from min/max (colormap value) to center (white)
r = repmat(abs(linspace(1,-1,64)),[3,1])';
mapPolar = mapPolar.*r + 1 - r;

Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];
Col4Dose2=[0.9 ,0., 0.55;...
0.83, 0.59, 0.03;...
0.11, 0.7, 0.41];

load([pn,'KrLineWithGenotypeKmeans.mat']);
Kr=Kr-min(Kr(:,101:900),[],2);
Gt=Gt-min(Gt(:,101:900),[],2);
Hb=Hb-min(Hb(:,101:900),[],2);
Kni=Kni-min(Kni(:,101:900),[],2);


ix4show=find(Age>38 & Age<46);
ix0x=ix4show(Genotype(ix4show)==0);
ix1x=ix4show(Genotype(ix4show)==1);
ix2x=ix4show(Genotype(ix4show)==2);
Kr=Kr/max(mean(Kr(ix2x,101:900)));
Gt=Gt/max(mean(Gt(ix2x,101:900)));
Hb=Hb/max(mean(Hb(ix2x,101:900)));
Kni=Kni/max(mean(Kni(ix2x,101:900)));

left=0.05;
bottom=0.05;
width=6.850394;
hight=2.5;

f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
b=0.185;
h=0.266;
w=0.15665;
ll=0.0705;
spacer=0.005;
axes('Position',[ll b w h])
hold on
plot((101:900)/1000,(Kr(ix0x,101:900)), 'Color',0.7*Col4Dose(1,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Kr(ix1x,101:900)), 'Color',0.7*Col4Dose(2,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Kr(ix2x,101:900)), 'Color',0.7*Col4Dose(3,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
% 
plot((101:900)/1000,mean(Kr(ix0x,101:900)), 'Color',Col4Dose(1,:), ...
    'LineStyle', '--', 'LineWidth', 1.);
plot((101:900)/1000,mean(Kr(ix1x,101:900)), 'Color', Col4Dose(2,:), 'LineStyle', '--', 'LineWidth', 0.8);
plot((101:900)/1000,mean(Kr(ix2x,101:900)), 'Color', Col4Dose(3,:), 'LineStyle', '--', 'LineWidth', 0.8);
xlabel('x/l');
ylabel('I^{Kr}');

box off 
axis tight
axis square

axes('Position',[ll+(w+spacer)*1 b w h])
hold on
plot((101:900)/1000,(Kni(ix0x,101:900)), 'Color',0.7*Col4Dose(1,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Kni(ix1x,101:900)), 'Color',0.7*Col4Dose(2,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Kni(ix2x,101:900)), 'Color',0.7*Col4Dose(3,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);

hold on
plot((101:900)/1000,mean(Kni(ix0x,101:900)), 'Color',Col4Dose(1,:), ...
    'LineStyle', '--', 'LineWidth', 1.);
plot((101:900)/1000,mean(Kni(ix1x,101:900)), 'Color', Col4Dose(2,:), 'LineStyle', '--', 'LineWidth', 0.8);
plot((101:900)/1000,mean(Kni(ix2x,101:900)), 'Color', Col4Dose(3,:), 'LineStyle', '--', 'LineWidth', 0.8);
xlabel('x/l');
ylabel('I^{Kni}');
box off 
axis tight
axis square

axes('Position',[ll+(w+spacer)*3, b, w, h])
hold on
plot((101:900)/1000,(Gt(ix0x,101:900)), 'Color',0.7*Col4Dose(1,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Gt(ix1x,101:900)), 'Color',0.7*Col4Dose(2,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Gt(ix2x,101:900)), 'Color',0.7*Col4Dose(3,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);

hold on
plot((101:900)/1000,mean(Gt(ix0x,101:900)), 'Color',Col4Dose(1,:), ...
    'LineStyle', '--', 'LineWidth', 1.);
plot((101:900)/1000,mean(Gt(ix1x,101:900)), 'Color', Col4Dose(2,:), 'LineStyle', '--', 'LineWidth', 0.8);
plot((101:900)/1000,mean(Gt(ix2x,101:900)), 'Color', Col4Dose(3,:), 'LineStyle', '--', 'LineWidth', 0.8);
xlabel('x/l');
ylabel('I^{Gt}');
box off 
axis tight
axis square
axes('Position',[ll+(w+spacer)*2, b w h])
hold on
plot((101:900)/1000,(Hb(ix0x,101:900)), 'Color',0.7*Col4Dose(1,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Hb(ix2x,101:900)), 'Color',0.7*Col4Dose(3,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);
plot((101:900)/1000,(Hb(ix1x,101:900)), 'Color',0.7*Col4Dose(2,:), ...
    'LineStyle', '-', 'LineWidth', 0.5);

hold on
plot((101:900)/1000,mean(Hb(ix0x,101:900)), 'Color',Col4Dose(1,:), ...
    'LineStyle', '--', 'LineWidth', 1.);
plot((101:900)/1000,mean(Hb(ix1x,101:900)), 'Color', Col4Dose(2,:), 'LineStyle', '--', 'LineWidth', 0.8);
plot((101:900)/1000,mean(Hb(ix2x,101:900)), 'Color', Col4Dose(3,:), 'LineStyle', '--', 'LineWidth', 0.8);
xlabel('x/l');
ylabel('I^{Hb}');
box off 
axis tight
axis square

axes('Position',[0.74835 0.1604 0.157 0.0635]);
imagesc((101:900)/1000,1,mean(Gt(ix0x,101:900))-mean(Gt(ix2x,101:900)), [-1 1]);
colormap(mapPolar)
xlabel('x/l');

cb=colorbar;
cb.Label.String='\Delta I_{\it 0xKr}^{Gt}';
set(gca, 'YTick', []);
cb.Position= [0.91541 0.15278 0.01034 0.0871];

axes('Position',[0.74835 0.406 0.157 0.307])
d=mean(Gt(ix0x,101:900))-mean(Gt(ix2x,101:900));
xp=(101:900)/1000;
plot(xp(d>0),d(d>0), 'Marker','.', 'MarkerSize', 1,'Color', 'r', 'LineStyle', 'none');
hold on
plot(xp(d<0),d(d<0),'Marker','.', 'MarkerSize', 1,'Color', 'b', 'LineStyle', 'none');
plot(xp,zeros(size(xp)), '--k');
axis tight
box off
ylabel('\Delta I_{\it 0xKr}^{Gt}');
xlabel('x/l')
%% F2 panels F-I
GAPNamesProtein={'Hb','Kr', 'Gt', 'Kni'};
gapNamesGenotype={'hb','Kr', 'gt', 'kni'};

Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];

steLevel=0;% setting steLevel>0 will color in white pixles for which expression difference is too low
StartT=6;% in min into NC 14
winSpan=8;% in minutes
progress=1;% in minutes
ageSpan=StartT:progress:(60-winSpan);% this is the time vector for dynamic reconstruction 

left=0.05;
bottom=0.05;
width=6.850394;
hight=10;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
t=42;
LW=0.5;
x=(101:900)/1000;
y=ageSpan+winSpan/2;
scalePar=0.7*y(end);
panSizeL=0.107;
panSizeW=0.14;

h=[0.05, 0.55,0.05, 0.55];
l=[0.05,0.05,0.55,0.55];
AxesFontSize=5;
LabelFontSize=8;
n=200;
Gr=[3,2,1];% to order the genotypes from 2x to 0x (left to right)
for lineID=[2,1,3,4]%
    fn=fullfile(pn,[GAPNamesProtein{lineID},'DiffData',num2str(n),'runs5']);
    load(fn);
      
     PickN=10;% this is the "wt" noise example from the bootstrap. everytime that bootstrap is rerun, this image would change.

    for group=1:3
        for gene=1:4
            POS=[l(lineID)+(Gr(group)-1)*panSizeW, h(lineID)+(gene-1)*panSizeL,panSizeW,panSizeL ];
            if group ==3
                 m=squeeze((diffGap(gene,PickN,:,:,group)));
            else
                m=squeeze(nanmean(diffGap(gene,:,:,:,group)));
            end             
                mdtSTD=squeeze(nanstd(diffGap(gene,:,:,:,group)));
                ix=abs(m)>steLevel*mdtSTD; 
                m(~ix)=0;
                ax1=axes('Position',POS);
                imagesc(ax1,x,y,m,[-1,1]);%, 'YDir', 'normal');
                colormap(mapPolar);

                if gene~=1
                  set(gca,'XTick',[]);
                 else
                  xlabel('x/l', 'FontSize', LabelFontSize);
                  set(gca,'XTick', 0.2:0.2:0.8);
                end
               if group~=3
                 set(gca,'YTick',[]);
               else
                ylabel('t[min]', 'FontSize', LabelFontSize);

                annotation('textbox',...
                [l(lineID)-0.05 h(lineID)-0.00+(gene)*panSizeL 0.03 0.01],...
                'String',GAPNamesProtein{gene},...
                'LineStyle','none',...
                'FitBoxToText','off', 'FontSize',6,'FontWeight','bold');
               end
        if lineID==1 && gene==4 && group==1
            cb=colorbar('Position', [0.4885 0.51085 0.02731 0.09830]);
            cb.Label.String='\DeltaI';%'mean(diff)';
            cb.Label.FontSize=6;
        end
    
    hold on
    plot(x,y(1)+scalePar*squeeze(nanmean(dtA(gene,:,t,:,3),2)), 'Color', [0 0 0 ], 'LineWidth',LW);%'Color',[0.5,0.5,0.5], 'LineWidth',LW);
    plot(x,y(1)+scalePar*squeeze(nanmean(dtB(gene,:,t,:,group),2)), 'Color',[0.6  0 0.7], 'LineWidth',LW);
  if gene==4
      title([num2str(group-1),'X',gapNamesGenotype{lineID}], 'FontSize', LabelFontSize,...
          'FontAngle', 'italic');
  end
    box off
    set(gca, 'YDir', 'normal','FontSize',AxesFontSize);%,'LabelFontSize',LabelFontSize );%,'XTick', [100 300 500 700],'XTickLabel', [0.2 0.4 0.6 0.8]);
   
        end
        Lef=l(lineID)+(group-1)*panSizeW;
        Hi= h(lineID);
        if group==3
            annotation('rectangle','Position', [Lef,Hi,panSizeW,panSizeL*4 ],...
                        'Color', Col4Dose(4-group,:), 'LineWidth', 1.5);
        else
            annotation('rectangle','Position', [Lef,Hi,panSizeW-0.003,panSizeL*4 ],...
                        'Color', Col4Dose(4-group,:), 'LineWidth', 1.5);
        end
    end
        Lef=l(lineID);
        Hi= h(lineID)+panSizeL*(lineID-1);
        annotation('rectangle','Position', [Lef,Hi,panSizeW*3,panSizeL ],...
                        'Color', 'k', 'LineWidth', 1.5, 'LineStyle', '--');

end
% % for saving high resolution figure  vesions:
% print(gcf,'Fig2.pdf','-dpdf','-r300');         %*// 300 dpi
 
%% Panel S2B 
%display the noise distribution from the wt bootstrap
S3=figure('Units', 'inches', 'Position',[0.2, 0.2,width,5]);
Edges=-0.1:0.001:0.1;
STE=nan(4,4); % this will contain the standard error values per line per gene (for the bootstrap 2xGAP)
for lineID=1:4%2
    fn=fullfile(pn,[GAPNamesProtein{lineID},'DiffData',num2str(n),'runs5']);
    load(fn);
      for gene=1:4%3
           subplot(4,4,4*(lineID-1)+gene )
           histogram(reshape(diffGap(gene,:,:,:,3),1,[]),'Normalization', 'pdf','BinEdges', Edges,...
               'FaceColor', Col4Dose(3,:), 'EdgeColor',  Col4Dose(3,:));
           box off
           xlabel('\deltaI_{wt}');
           ylabel('pdf');
           ylim([0 220]);
           title([GAPNamesProtein{gene}, '@\it 2x', (gapNamesGenotype{lineID})]);
           STE(lineID,gene)=nanstd(reshape(diffGap(gene,:,:,:,3),1,[]));
      end
    
end
maxNoiseSingleSigma=max(STE, [],'all')
steGtAtKr=STE(2,3)
%% Panel S2C
% a version of Figure 2 F-I with blackened pixles where data is lower than noise limit
steLevel=0;
StartT=6;
winSpan=8;
progress=1;
ageSpan=StartT:progress:(60-winSpan);
left=0.05;
bottom=0.05;
width=6.850394;
hight=10;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);

t=42;
x=(101:900)/1000;
y=ageSpan+winSpan/2;

scalePar=0.7*y(end);
panSizeL=0.107;
panSizeW=0.14;

h=[0.05, 0.05,0.55, 0.55];
l=[0.05,0.55,0.05,0.55];
AxesFontSize=5;
LabelFontSize=8;
LW=0.5;
n=200;

Gr=1:3;%[3,2,1];% to order the genotypes from 2x to 0x (left to right)
for lineID=2%[2,1,3,4]%
    fn=fullfile(pn,[GAPNamesProtein{lineID},'DiffData',num2str(n),'runs5']);
    load(fn);
      
    PickN=10;
    for group=1:3
        for gene=3%1:4
            POS=[l(lineID)+(Gr(group)-1)*panSizeW, h(lineID)+(gene-1)*panSizeL,panSizeW,panSizeL ];
            if group ==3
                 m=squeeze((diffGap(gene,PickN,:,:,group)));
            else
                m=squeeze(nanmean(diffGap(gene,:,:,:,group)));
            end
            
                mdtSTD=squeeze(nanstd(diffGap(gene,:,:,:,3),[],'all'));
                ix=abs(m)>2*mdtSTD;
                ax1=axes('Position',POS);
                ax2=axes('Position',POS);
                coverIX=zeros(size(m));
                coverIX(ix)=1;% gray colormap is indicating zero with black
                 
                imagesc(ax1,x,y,m,[-1,1]);%, 'YDir', 'normal');
                colormap(ax1,mapPolar);
                hold on
                plot(x,y(1)+scalePar*squeeze(nanmean(dtA(gene,:,t,:,3),2)), 'Color', [0 0 0 ], 'LineWidth',LW);%'Color',[0.5,0.5,0.5], 'LineWidth',LW);
                plot(x,y(1)+scalePar*squeeze(nanmean(dtB(gene,:,t,:,group),2)), 'Color',[0.6  0 0.7], 'LineWidth',LW);
                xlim([0.1 0.9]);
                ylim([10, 56]);
                if lineID==1 && gene==4 && group==1
                cb=colorbar(ax1,'Position', [0.49 0.497 0.0162 0.04]);
                cb.Label.String='\DeltaI';%'mean(diff)';
                text('FontSize',6,'String','\DeltaI',...
                'Position',[-0.0999 1.30 0]);

                end
    
                imagesc(ax2,x,y,coverIX,[0,1]);%, 'YDir', 'normal');  
                 colormap(ax2,'gray');
                
                linkaxes([ax1,ax2])
%                 %%Hide the top axes
                ax2.Visible = 'off';
                ax2.XTick = [];
                ax2.YTick = [];
               
                alpha(ax2, 0.5);
                if gene~=1
                  set(ax1,'XTick',[]);
                 else
                  xlabel(ax1, 'x/l', 'FontSize', LabelFontSize);
                  set(ax1,'XTick', 0.2:0.2:0.8);
                end
               if group~=1
                 set(ax1,'YTick',[]);
               else
                ylabel(ax1, 't[min]', 'FontSize', LabelFontSize);
                annotation('textbox',...
                [l(lineID)-0.05 h(lineID)-0.00+(gene)*panSizeL 0.03 0.01],...
                'String',GAPNamesProtein{gene},...
                'LineStyle','none',...
                'FitBoxToText','off', 'FontSize',6,'FontWeight','bold');
               end
                
        
  if gene==4
      title(ax1,[num2str(group-1),'x',gapNamesGenotype{lineID}], 'FontSize', LabelFontSize,...
          'FontAngle', 'italic');
  end
    box off
    set(ax1, 'YDir', 'normal','FontSize',AxesFontSize);%,'LabelFontSize',LabelFontSize );%,'XTick', [100 300 500 700],'XTickLabel', [0.2 0.4 0.6 0.8]);
   
        end
        Lef=l(lineID)+(group-1)*panSizeW;
        Hi= h(lineID);
        if group==3
            annotation('rectangle','Position', [Lef,Hi,panSizeW,panSizeL*4 ],...
                        'Color', Col4Dose(group,:), 'LineWidth', 1.5);
        else
            annotation('rectangle','Position', [Lef,Hi,panSizeW-0.003,panSizeL*4 ],...
                        'Color', Col4Dose(group,:), 'LineWidth', 1.5);
        end
    end
    Lef=l(lineID);
        Hi= h(lineID)+panSizeL*(lineID-1);
            annotation('rectangle','Position', [Lef,Hi,panSizeW*3,panSizeL ],...
                        'Color', [0.1 0 0.8], 'LineWidth', 1.5, 'LineStyle', '--');
end
%%  Panel S2D 
tRange=[38 46];% The time window of interest (in minutes into NC14)
DT=cell(4,3);% this will contain the data at trange per lineID (n=4) and genotype(n=3)
% prepare data:
for lineID=1:4
    fn=fullfile(pn, [GAPNamesProtein{lineID},'LineWithGenotypeKmeans']);
    load (fn);

    % normalize data and identify 2x gap embryos in the 38-48 min window:
    %remove background from individual embryos and normalize to mean(wt(38-48))
    Hb=Hb-min(Hb(:,101:900), [],2);
    Kr=Kr-min(Kr(:,101:900), [],2);
    Gt=Gt-min(Gt(:,101:900), [],2);
    Kni=Kni-min(Kni(:,101:900), [],2);

    for idg=0:2 % take each genotype: 0x, 1x and 2x gap
        idxx=find(Genotype==idg);
        id2x=find(Genotype==2);
        idAge=find(Age>tRange(1) & Age<tRange(2));
        ix=intersect(idxx,idAge);
        ixWT=intersect(id2x,idAge);
        HbMut=Hb(ix,:);
        KrMut=Kr(ix,:);
        GtMut=Gt(ix,:);
        KniMut=Kni(ix,:);

        HbWT=Hb(ixWT,:);
        KrWT=Kr(ixWT,:);
        GtWT=Gt(ixWT,:);
        KniWT=Kni(ixWT,:);


        HbMut=HbMut/max(mean(HbWT(:,101:900)));
        KrMut=KrMut/max(mean(KrWT(:,101:900)));
        GtMut=GtMut/max(mean(GtWT(:,101:900)));
        KniMut=KniMut/max(mean(KniWT(:,101:900)));

        HbWT=HbWT/max(mean(HbWT(:,101:900)));
        KrWT=KrWT/max(mean(KrWT(:,101:900)));
        GtWT=GtWT/max(mean(GtWT(:,101:900)));
        KniWT=KniWT/max(mean(KniWT(:,101:900)));

        
        nEmbryos=size(KrMut,1);
        dt=nan(4,nEmbryos,1000);
        for n=1:nEmbryos
            dt(1,n,101:900)=HbMut(n,101:900);
            dt(2,n,101:900)=KrMut(n,101:900);
            dt(3,n,101:900)=GtMut(n,101:900);
            dt(4,n,101:900)=KniMut(n,101:900);
        end
        DT{lineID,idg+1}=dt;
    end
end
% get sigma I
SigmaI=nan(4, 1000,4,3);
gI=nan(4, 1000,4,3);
for lineID=1:4
    for idg=1:3
        for g=1:4
           if idg==1 && g==lineID
               continue
           end
           SigmaI(g,101:900,lineID,idg)=squeeze(nanstd(DT{lineID,idg}(g,:,101:900)));
           gI(g,101:900,lineID,idg)=squeeze(nanmean(DT{lineID,idg}(g,:,101:900)));

        end
      
    end
end
% plot Panel S2D  

left=0.15;
bottom=0.15;
width=6.8;
hight=3;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
for g=1:4
    subplot(2,4,g);% plot SigmaI along the egg per gene
    plot((101:900)/1000,squeeze(nanmean(dt(g,:,101:900))), 'Color',[0.5 0.5 0.5]);
    hold on
    for lineID=1:4
        for idg=1:2
             if idg==1 && lineID==g
                        continue
             end
              plot((1:1000)/1000,SigmaI(g,:,lineID,idg)./gI(g,:,lineID,idg), 'Color', Col4Dose(idg,:),...
                  'Marker', '.', 'LineStyle','none', 'MarkerSize', 0.5);
              hold on
              plot((1:1000)/1000,SigmaI(g,:,lineID,3)./gI(g,:,lineID,3), 'Color',Col4Dose(3,:),...
                  'Marker', '.', 'LineStyle','none', 'MarkerSize', 0.5);
              hold on   
        end
            box off
            xlabel('x/l');
            ylabel('\sigma_i/g_i');
            xlim([0.1 0.9]);
            ylim([0 2]);
            title(GAPNamesProtein{g});
    end
end

for g=1:4
    subplot(2,4,g+4);% plot the mean expression intensity vs. its std
    for lineID=1:4
        for idg=1:2    
            if idg==1 && lineID==g
                continue
            end
              plot(gI(g,:,lineID,idg),SigmaI(g,:,lineID,idg), 'Color', Col4Dose(idg,:),...
                  'Marker','.', 'LineStyle','none', 'MarkerSize', 0.5);
              hold on
              plot(gI(g,:,lineID,3),SigmaI(g,:,lineID,3), 'Color',Col4Dose(3,:),...
                  'Marker','.', 'LineStyle','none', 'MarkerSize', 0.5);
              hold on

        end
        box off
        xlabel('g_i');
        ylabel('\sigma_i');
        ylim([0 0.35]);
    end
end
p0=plot(nan, nan, 'Color', Col4Dose(1,:), 'Marker','.', 'LineStyle', 'none');
p1=plot(nan, nan, 'Color', Col4Dose(2,:), 'Marker','.', 'LineStyle', 'none');
p2=plot(nan, nan, 'Color', Col4Dose(3,:), 'Marker','.', 'LineStyle', 'none');
legend([p0,p1,p2],'0xGap','1xGap','2xGap', 'Marker','.', 'LineStyle', 'none',...
    'Position',[0.009 0.43 0.12 0.15]);
%% Panels E&F in S2
% For sup figure of hets sigmai - show sigmax for 4g and best 1g? get the 1g positional error for the hets for the genes that increase their levels:
tRange=[38 46];
zsmooth=30; %
nbins=1000;
% use data from DT=cell(4,3);% this will contain the data at tRange per lineID (n=4) and genotype(n=3)

cGAP=[ 0.9960,    0.8200    ,0.6040;...
    0.6120    ,0.9540    ,0.8480;...
    0.4480   , 0.5960  ,  0.9000;...
    0.8360    ,0.5640 ,   0.7120];

cEve=colormap(jet(7));

Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];
GAPNamesProtein={'Hb','Kr', 'Gt', 'Kni'};
sigma_4Gm=nan(10, nbins,4);%10 replicas, nbins along the egg, 4 lines
sigma_2Gm=nan(10,nbins,6,4);%10 replicas, nbins along the egg, 6 combinations of pairs, 4 lines
sigma_1Gm=nan(10,nbins,4,4);%10 replicas, nbins along the egg, 4 genes, 4 lines
sigma_3Gm=nan(10,nbins,4,4);%10 replicas, nbins along the egg, 4 combinations of triplets, 4 lines
for lineID=1:4
     % test sigma of wt against Previous datasets 
    dt=DT{lineID,3};
    counter=0;
    [sigma_4Gm(:,:,lineID) mgm dgm cgm]  = compute_SigmaX(dt, nbins, zsmooth);
    for g1=1:4
        for g2=(g1+1):4
            counter=counter+1;
     [sigma_2Gm(:,:,counter,lineID) mgm dgm cgm]  = compute_SigmaX(dt([g1,g2],:,:), nbins, zsmooth);
        end
    end
    for g=1:4
        Genes=1:4;
        Genes(g)=[];
        [sigma_1Gm(:,:,g,lineID) mgm dgm cgm]  = compute_SigmaX(dt(g,:,:), nbins, zsmooth);
        [sigma_3Gm(:,:,g,lineID) mgm dgm cgm]  = compute_SigmaX(dt(Genes,:,:), nbins, zsmooth);
    end
end

% now get the positional error for the gap mutants:
sigma_1GmMut=nan(10,1000,4,4,2);
sigma_2GmMut=nan(10,1000,6,4,2);
sigma_3GmMut=nan(10,1000,4,4,2);
sigma_4GmMut=nan(10,1000,4,2);
for lineID=1:4
     for group=1:2
          dt=DT{lineID,group};
          [sigma_4GmMut(:,:,lineID,group) mgm dgm cgm]=...
              compute_SigmaX(dt, nbins, zsmooth);
            for g1=1:4           
                for g2=(g1+1):4
                    if (group==1 && g1==lineID) || (group==1 && g2==lineID)
                        continue
                    end

                    counter=counter+1;
                    [sigma_2GmMut(:,:,counter,lineID,group) mgm dgm cgm]=...
                        compute_SigmaX(dt([g1,g2],:,:), nbins, zsmooth);
                end
            end
            
            for g=1:4
                if group==1 && g==lineID
                    sigma_1GmMut(:,:,g,lineID,group)=nan;
                else
                    [sigma_1GmMut(:,:,g,lineID,group) mgm dgm cgm]  = compute_SigmaX(dt(g,:,:), nbins, zsmooth);

                end    
                Genes=1:4;
                Genes(g)=[];
                [sigma_3GmMut(:,:,g,lineID,group) mgm dgm cgm]  = compute_SigmaX(dt(Genes,:,:), nbins, zsmooth);
            end
     end
end

% get the single gene that provides the lowest sigmax1g per position for
% each of the mutants:
Sigma1xBest=nan(900,4,2);%(x,lineID, group);
BestGene=nan(900,4,2);%(x,lineID, group);
    for group=1:2
        for lineID=1:4
            ixG=1:4;
            if group==1 
                ixG(lineID)=[];
            end
              Ms1x=squeeze(mean(sigma_1GmMut(:,:,:,lineID,group)));
                for x=101:900
                    [Sigma1xBest(x,lineID, group),ix]=min(Ms1x(x,ixG));
                    BestGene(x,lineID, group)=ixG(ix);
                end
        end
    
    end
% plot sigmax for both hets and nulls, just for sigma1xBestMut and sigma4xmut
left=0.15;
bottom=0.15;
width=6.8;
hight=4;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
for lineID=1:4
    fn=fullfile(pn,['PRpeakNvalley4',GAPNamesProtein{lineID}]);
    load (fn, 'EvePeakLocNulls' );
    for group=1:2
        if group==2
            subplot(2,4,lineID);
        else
            subplot(2,4,lineID+4*(group));
        end
        p3=plot((1:1000)/1000,0.06*squeeze(nanmean(dt(lineID,:,:))), 'Color',[0.5 0.5 0.5]);
        hold on
        p2x=plot((1:1000)/1000,mean(sigma_4Gm(:,:,lineID)), 'Color', Col4Dose(3,:),...
             'LineWidth', 1.2);
        if group==1
            p0x=plot((1:1000)/1000,mean(sigma_3GmMut(:,:,lineID,lineID,group)),...
                'Color', Col4Dose(group,:));
            hold on
             p2x3=plot((1:1000)/1000,mean(sigma_3Gm(:,:,lineID, lineID)), 'Color', 'b');
                 
            hold on
        else
            p1x=plot((1:1000)/1000,mean(sigma_4GmMut(:,:,lineID,group)),...
                'Color', Col4Dose(group,:));
            hold on
            p0x=plot(nan, nan, 'Color', Col4Dose(1,:));
        end
        hold on
   
        p4=plot ((101:900)/1000,Sigma1xBest(101:900,lineID, group),'--k');
        hold on
        plot([0.1 0.9 ], [0.015, 0.015], '--k');

        xlim([0.1,0.9]);
        ylim([0 0.06])
        box off
        xlabel('x/l');
        ylabel('\sigma_x')

        if lineID==1 && group==2
            legend([p0x,p1x, p2x, p2x3,p4],{'\sigma_{x3g}0xGap', '\sigma_{x4g}1xGap', '\sigma_{x4g}2xGap',...
            '\sigma_{x3g}2xGap', 'min|_g(\sigma_{x1g}mut)'},'Position',[-0.0387 0.429 0.176 0.233]);

        end
        axis square
        set(gca, 'FontSize', 8);

    end
end
