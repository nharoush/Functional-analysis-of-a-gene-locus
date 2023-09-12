%% Computes positional error carried by nG genes. Assumes that at every position, the
%% gene expression variability is well described by a multivariate gaussian. Bootstraps 
%% the result be repeatedly estimating the error with leave-one-out procedure.
%
%   
%   Julien Dubuis, Gasper Tkacik (2008, 2010, 2014)
%
%   INPUT:
%       data 
%           nG x nN x nX matrix containing samples of expression profiles
%           for nG genes from nN embroys, sampled at positions nX
%           data should be normalized such that the mean profile is between
%           0 and 1 (individual profiles can fluctuate out of these bounds)
%       nbins 
%           the new discretization to use along the position axis; nbins
%           has to divide nX without a reminder; consecutive bins in the
%           original discretization will be averaged together
%       zsmooth
%           smoothing to apply to the mean profiles before computing the
%           numerical derivative; indicates the number of consecutive bins
%           to smooth over; 1 = no smoothing
%    OUTPUT:
%       sigmaX
%           estimate of positional error, of dimension (nshuffles) x
%           (nbins), containing nshuffles leave-one-out estimates
%       meanG
%           the mean of the profiles using the new smoothing and
%           discretization; dimension nG x nbins
%       dGdX
%           numerical spatial derivative of the gene expression patterns;
%           dimension nG x nbins
%       covG
%           the covariance matrix of expression levels of dimension nG x nG
%           x nbins

function [sigmaX meanG dGdX covG] = compute_SigmaX(data, nbins, zsmooth)

    if (nargin <= 1) 
        nbins = 100; 
        zsmooth = 1;
    end;
    
    nshuffles = 10;
    
    if (numel(size(data)) == 2)
        data = reshape(data, 1, size(data,1), size(data, 2));
    end
    
    [nG nN nX] = size(data);
   
    if (mod(nX, nbins)~=0)
        error('The new number of spatial bins must divide the original number.');
    end
    
    disp(sprintf('**** COMPUTE_SIGMAX: Estimating positional error of %d genes from %d embryos.', nG, nN));
    disp(sprintf('**** COMPUTE_SIGMAX: %d spatial bins -> %d spatial bins, smoothing over %d bins.', nX, nbins, zsmooth));
    
    
    nK = nX/nbins;

    dX = 1./nbins;

    meanG = zeros(nG, nshuffles, nbins);
    
    covgg=cell(nshuffles, nbins);

    divgg=cell(nshuffles, nbins);
    
    sigmax4G=NaN(nshuffles, nbins);

    for kk=1:nshuffles

        samps=randperm(nN);
        samps=samps(1:nN-1);

        gg = data(:,samps,:);
        
        for l=1:nbins,
            
            meanG(:,kk,l) = squeeze(nanmean(nanmean(gg(:,:,(l-1)*nK+1:l*nK),2),3));
            
            for q=1:nG,
                tt(q,:)   = flatten(gg(q,:,(l-1)*nK+1:l*nK));
            end
            
            covgg{kk,l}=nancov(tt');

        end

        if (~isempty(zsmooth))
            for q=1:nG,
                meanG(q,kk,:) = smooth(meanG(q,kk,:),zsmooth);
            end
        end

        divgg{kk,1} = (meanG(:,kk,2) - meanG(:,kk,1))/dX;
        divgg{kk,nbins} = (meanG(:,kk,nbins) - meanG(:,kk,nbins-1))/dX;
       
        for l=2:nbins-1,
            divgg{kk,l} = (meanG(:,kk,l+1)-meanG(:,kk,l-1))/(2*dX);
        end

        
        %size(divgg{kk,l})
        for l=1:nbins,   
            if (sum(flatten(isnan(covgg{kk,l}))) > 0 || sum(isnan(divgg{kk,l}))>0)
                sigmaX(kk,l) = nan;
            else
                sigmaX(kk,l) = 1/sqrt(divgg{kk,l}'*covgg{kk,l}^(-1)*divgg{kk,l});
            end
        end
    end

      for l=1:nbins,
           covG(:,:,l) = covgg{1,l}; 
           dGdX(:,l) = divgg{1,l};
      end
      meanG = squeeze(mean(meanG(:,:,:),2));
end