%% get the mean fine time trace with the same time line for all:
function[MeanHbFine,MeanKrFine,MeanGtFine,MeanKniFine,stdHbFine,stdKrFine,...
    stdGtFine,stdKniFine]=RcnstrctGapTseriesWithStd4RandSplit(Hb,Kr,Gt,Kni,Age,Genotype,...
    progress,winSpan,StartT)
Counter=0;
    for i=StartT:progress:(60-winSpan)
        Counter=Counter+1
        idx = find(Age>=i & Age<=i+winSpan);
        id2x=find(Genotype(idx)==2);
        idHets=find(Genotype(idx)==1);
        idNulls=find(Genotype(idx)==0);
        
       MeanHbFine(Counter,:,1)=nanmean(Hb(idx(idNulls),101:900),1);
       MeanGtFine(Counter,:,1)=nanmean(Gt(idx(idNulls),101:900),1);
       MeanKrFine(Counter,:,1)=nanmean(Kr(idx(idNulls),101:900),1);
       MeanKniFine(Counter,:,1)=nanmean(Kni(idx(idNulls),101:900),1);
        
        MeanHbFine(Counter,:,2)=nanmean(Hb(idx(idHets),101:900),1);
        MeanHbFine(Counter,:,3)=nanmean(Hb(idx(id2x),101:900),1);

        MeanGtFine(Counter,:,2)=nanmean(Gt(idx(idHets),101:900),1);
        MeanGtFine(Counter,:,3)=nanmean(Gt(idx(id2x),101:900),1);
        
        MeanKrFine(Counter,:,2)=nanmean(Kr(idx(idHets),101:900),1);
        MeanKrFine(Counter,:,3)=nanmean(Kr(idx(id2x),101:900),1);

        MeanKniFine(Counter,:,2)=nanmean(Kni(idx(idHets),101:900),1);
        MeanKniFine(Counter,:,3)=nanmean(Kni(idx(id2x),101:900),1);
        
        
        stdHbFine(Counter,:,1)=nanstd(Hb(idx(idNulls),101:900),[],1);
        stdHbFine(Counter,:,2)=nanstd(Hb(idx(idHets),101:900),[],1);
        stdHbFine(Counter,:,3)=nanstd(Hb(idx(id2x),101:900),[],1);

        stdGtFine(Counter,:,1)=nanstd(Gt(idx(idNulls),101:900),[],1);
        stdGtFine(Counter,:,2)=nanstd(Gt(idx(idHets),101:900),[],1);
        stdGtFine(Counter,:,3)=nanstd(Gt(idx(id2x),101:900),[],1);

        stdKrFine(Counter,:,1)=nanstd(Kr(idx(idNulls),101:900),[],1);
        stdKrFine(Counter,:,2)=nanstd(Kr(idx(idHets),101:900),[],1);
        stdKrFine(Counter,:,3)=nanstd(Kr(idx(id2x),101:900),[],1);

        stdKniFine(Counter,:,1)=nanstd(Kni(idx(idNulls),101:900),[],1);
        stdKniFine(Counter,:,2)=nanstd(Kni(idx(idHets),101:900),[],1);
        stdKniFine(Counter,:,3)=nanstd(Kni(idx(id2x),101:900),[],1);
    end
end

