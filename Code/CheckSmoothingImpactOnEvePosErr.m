
%% check what is the difference if I smooth Eve profiles by 2%, 5% or 10% of egg lenght
%%%% find eve peaks for all datasets 
nP=[6, 6, 7, 6];
% % gen=2; % gen=0 gap null mut, gen =1 gap hets or maternal mut, gen=2 wt 
% minIPI=[40,40, 40, 30]; % minimal inter peak interval, per line, adjusted manually
minIPI=40;%40;
xix=205:920;
%  xix=300:820;
SmoothSpan=[20, 50, 100];%this is the number of AP pixels for smoothing the Eve trace for finding peaks and troughs 
for ss=1:lenght(SmoothSpan)
for LineID=1:4
    for gen=0:2
        ix=GenotypeEveAll{LineID}==gen;
        Eve=EveAll{LineID}(ix,:);
        nEm=sum(ix);
        EvePeakIntens=NaN(nEm,nP(LineID));
        EvePeakLoc=NaN(nEm,nP(LineID));

        for i=1:nEm%length(idwt)/2+1:length(idwt)%1:length(idwt)/2
                 Eve(i,:)=smooth(Eve(i,:), 50);
                [pck2x1,loc2x1]=findpeaks(Eve(i,xix),...
                         'SortStr', 'descend', 'MinPeakDistance', minIPI);%, 'MinPeakHeight', -0.2

            pck2x1(nP(LineID)+1:end)=[];
            loc2x1(nP(LineID)+1:end)=[];
        %                  [pck2x1,loc2x1]=findpeaks(Eve(id2x(i),xix),...
        %                  'MinPeakProminence',mp,'MinPeakDistance',40,'NPeaks',nP);    
                figure; plot(Eve(i,xix));
                  hold on
                  plot(loc2x1,pck2x1,'*r');
                  loc2x1=loc2x1'+xix(1)-1;
                   if ~isempty(loc2x1)
                          index=sortrows([loc2x1,(1:length(loc2x1))']);    
                            index(:,1)=[];
                          EvePeakIntens(i,1:length(loc2x1))=pck2x1(index);
                           EvePeakLoc(i,1:length(loc2x1))=loc2x1(index);    

                       %                EvePeakIntensWT(i,:)=sortrows(pck2x1');
        %                EvePeakLocWT(i,:)=sortrows(loc2x1')+xix(1)-1;
                      else
                       EvePeakIntens(i,:)=nan;
                       EvePeakLoc(i,:)=nan;
                       continue 
                   end

        end
        lineID
        gen
        nanstd(EvePeakLoc/1000)
        
        if gen==0
            % now make sure that flickering stripes are assigned with their correct
            % serial identity:
            % nNulls=size(EvePeakLoc, 1);
            a=~isnan(EvePeakLoc);
            sN=max(sum(a,2));% max num of stripes identified
            %      find embryos that show all stripes:
                emLine=~isnan(EvePeakLoc(:,sN));
                if sum(emLine)==1
                    meanPeakLocNulls=(EvePeakLoc(emLine,:));

                else
                 meanPeakLocNulls=nanmean(EvePeakLoc(emLine,:));
                end
                 EvePeakIntensNulls1=nan(nEm, 7);
                 EvePeakLocNulls1=nan(nEm, 7);
                for i=1:nEm
                    locNulls1=EvePeakLoc(i,~isnan(EvePeakLoc(i,:)));
                           for k=1:length(locNulls1)
                               dP=abs(meanPeakLocNulls-locNulls1(k));
                               [m,ind]=min(dP);
                               EvePeakIntensNulls1(i,ind)=EvePeakIntens(i,k);
                               EvePeakLocNulls1(i,ind)=locNulls1(k);    
                           end
                end  
                EvePeakLocAll{lineID, gen, ss}=EvePeakLocNulls1;
                EvePeakIntensAll{lineID, gen, ss}=EvePeakIntensNulls1;
%                 EveAll{lineID, gen, ss}=Eve
            gen
            nanstd(EvePeakLocNulls1/1000)
        else
            EvePeakLocAll{lineID, gen, ss}=EvePeakLoc;
            EvePeakIntensAll{lineID, gen, ss}=EvePeakIntens;
    %                 EveAll{lineID, gen, ss}=Eve
        end

        end
        

   end
end

