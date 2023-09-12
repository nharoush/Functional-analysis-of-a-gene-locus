clearvars;
% if ~exist('processedData', 'dir')
%     mkdir('processedData')
% end
fol=cd;
pn = fullfile(fol,'processedData');

cGAP=[ 0.9960,    0.8200    ,0.6040;...
    0.6120    ,0.9540    ,0.8480;...
    0.4480   , 0.5960  ,  0.9000;...
    0.8360    ,0.5640 ,   0.7120];

Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];
GapNameS={'Hb','Kr','Gt','Kni'};
gapNameS={'hb','Kr','gt','kni'};

MD0xAll=[];
MD1xAll=[];
MD2xAll=[];
minMD=cell(1,4);
minMD1x=cell(1,4);
minMD2x=cell(1,4);
minMDind=cell(1,4);
minMDind1x=cell(1,4);
minMDind2x=cell(1,4);
meanDAll=cell(1,4);

 for lineID=1:4
     fn=fullfile(pn,[GapNameS{lineID},'LineWithGenotypeKmeans.mat']);
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

    % calc chi square/mahal, only around 42 min? % mahal and chi square a-la-JD are the same thing...

    dt2x=dt(:,Genotype==2 & Age>=42-4 & Age<42+4,:);
    dt1x=dt(:,Genotype==1 & Age>=42-4 & Age<42+4,:); % for checking wt from wt
    dt0x=dt(:,Genotype==0 & Age>=42-4 & Age<42+4,:);


        MD=nan(1000,1000,size(dt0x,2));
        MD1x=nan(1000,1000,size(dt1x,2));
        MD2x=nan(1000,1000,size(dt2x,2));
        
    for x1=101:1:900
        for x2=101:1:900
            MD(x1,x2,:)=(mahal(squeeze(dt0x(:,:,x1))',squeeze(dt2x(:,:,x2))'));
          
            MD1x(x1,x2,:)=(mahal(squeeze(dt1x(:,:,x1))',squeeze(dt2x(:,:,x2))'));
            MD2x(x1,x2,:)=(mahal(squeeze(dt2x(:,:,x1))',squeeze(dt2x(:,:,x2))'));
        end
          [minMD{lineID}(x1,:),minMDind{lineID}(x1,:)]=min(squeeze(MD(x1,:,:)));
          [minMD1x{lineID}(x1,:),minMDind1x{lineID}(x1,:)]=min(squeeze(MD1x(x1,:,:)));
          [minMD2x{lineID}(x1,:),minMDind2x{lineID}(x1,:)]=min(squeeze(MD2x(x1,:,:)));

          
    end
        MD0xAll=[MD0xAll;minMD{lineID}(~isnan(minMD{lineID}))];
        MD1xAll=[MD1xAll;minMD1x{lineID}(~isnan(minMD1x{lineID}))];
        MD2xAll=[MD2xAll;minMD2x{lineID}(~isnan(minMD2x{lineID}))];
 end
 MD0xAll=MD0xAll/4;
MD2xAll=MD2xAll/4;
MD1xAll=MD1xAll/4;


 MaxWt=max(MD2xAll,[],'all');
% get eve peak poistions for nulls:
for lineID=1:4
    fn=fullfile(pn,['PRpeakNvalley4', GapNameS{lineID}]);
    load(fn);
    peakLocsNulls{lineID}=EvePeakLocNulls;
end
 
% plot:
left=0.15;
bottom=0.15;
width=6.8;
hight=4;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);

 for lineID=1:4
            subplot(2,4,lineID); 
            plot((101:10:900)/1000,(minMD{lineID}(101:10:900,:))/4, 'Color', [0.7 0.7 0.7]);
            hold on
            plot((101:10:900)/1000,(mean(minMD{lineID}(101:10:900,:),2))/4, '--k');
            plot((101:900)/1000, MaxWt*ones(1,800), '--k');

            hold on
            plot((101:900)/1000,15*squeeze(nanmean(dt2x(lineID,:,101:900))),'Color', cGAP (lineID,:), 'LineWidth', 1.2);
            hold on
            xlabel('x/l')
            ylabel('\chi^2_{per gene}');
            sp=(nanmean(peakLocsNulls{lineID}))/1000;

            for s=1:length(sp)
                plot([sp(s),sp(s)], [0,50],'Color', Col4Dose(1,:), 'LineStyle', '--');
                hold on
            end
            title(['0x',gapNameS{lineID}], 'FontAngle', 'italic');
            xlim([0.1,0.9]);
            ylim([0,15]);
            box off

            subplot(2,4,4+lineID); 
            plot((101:10:900)/1000,(minMD1x{lineID}(101:10:900,:))/4, 'Color', [0.7 0.7 0.7]);
            hold on
            plot((101:10:900)/1000,(mean(minMD1x{lineID}(101:10:900,:),2))/4, '--k');
            plot((101:900)/1000, MaxWt*ones(1,800), '--k');
            ylabel('\chi^2_{per gene}');
      
            xlabel('x/l')
            title(['1x',gapNameS{lineID}]);
    %         ylim([0,max(minMD, [],'all')]);
            ylim([0,15]);
            xlim([0.1,0.9]);
            box off
        
 end
        set (gca, 'FontSize', 8);
 %%
 Col4Dose=[1 ,0.078, 0.65;...
0.93, 0.69, 0.13;...
0.21, 0.8, 0.51];

clear MD MD2x
[counts,edges]=histcounts(MD0xAll,1:0.1:max(MD0xAll,[],'all'));
[counts2x,edges2x]=histcounts(MD2xAll,1:0.01:max(MD2xAll,[],'all'));
[counts1x,edges1x]=histcounts(MD1xAll,1:0.1:max(MD1xAll,[],'all'));
 left=0.15;
bottom=0.15;
width=4.5;
hight=2;
f=figure('Units', 'inches', 'Position',[left, bottom,width,hight]);
% subplot(1,5,2:4)
    semilogx(edges(1:end-1)+diff(edges(1:2))/2,cumsum(counts, 'reverse')/sum(counts),...
        'Color', Col4Dose(1,:));
%   plot(edges(1:end-1)+diff(edges(1:2))/2,cumsum(counts, 'reverse')/sum(counts),...
%         'Color', Col4Dose(1,:));
    hold on
     semilogx(edges1x(1:end-1)+diff(edges1x(1:2))/2,...
         cumsum(counts1x, 'reverse')/sum(counts1x),'Color', Col4Dose(2,:));
% plot(edges1x(1:end-1)+diff(edges1x(1:2))/2,...
%          cumsum(counts1x, 'reverse')/sum(counts1x),'Color', Col4Dose(2,:));
hold on
     semilogx(edges2x(1:end-1)+diff(edges2x(1:2))/2,...
         cumsum(counts2x, 'reverse')/sum(counts2x), 'Color', Col4Dose(3,:));
% plot(edges2x(1:end-1)+diff(edges2x(1:2))/2,...
%          cumsum(counts2x, 'reverse')/sum(counts2x), 'Color', Col4Dose(3,:));

     hold on
     semilogx([MaxWt,MaxWt,MaxWt], [10^(-5) 0.5 1], '--k');
%      semilogx([4,4,4], [10^(-5) 0.5 1], '--r');
% plot([MaxWt,MaxWt,MaxWt], [10^(-5) 0.5 1], '--k');
% plot([3,3,3], [10^(-5) 0.5 1], '--r');

     legend({['0xGap'], ['1xGap'],['2xGap']}, 'FontAngle', 'italic');
    ylabel('Cumulative probability');
    xlabel('\chi^2_{per gene}');
% xlabel('$\sqrt{\chi^2}$', 'Interpreter', 'latex');
set (gca, 'FontSize', 8);
box off
axis tight

%%
