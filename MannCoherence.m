function [k1,cocoh,quadcoh] = MannCoherence(PHI,k11,k2,k3,d2,d3,k2_log,k3_log)
% [k1,coh] = MannCoherence(PHI,k11,k2,k3,d2,d3,k2_log,k3_log) computes wind
% coherence based on Mann spectral tensor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHI : 5D Spectral tensor for the three wind components  [2Ndk1 x 2Ndk2 x 2Ndk3 x 3 x 3]
% k2 : wavenumber [2 x Ndk1, 2 x Ndk2, 2 x Ndk3]
% k3 : wavenumber [2 x Ndk1, 2 x Ndk2, 2 x Ndk3]
% k2_log : two sided wavenumber [1, 2 x Ndk2]
% k3_log : two sided wavenumber [1, 2 x Ndk3]
% d2 : Single distance along an horizontal line [1 x 1]
% d3 : Single distance along a vertical line [1 x 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k1: wavenumver in the along-wind direction : [1 x Ndk1]
% cocoh: co-coherence calculated for d2 and d3
% quadcoh: quad-coherence calculated for d2 and d3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author: Etienne Cheynet -- updated on : 16/04/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see also MannTurb.m fitMannTensor.m

%%
N1=size(PHI,1);
k1 = interp(k11,5);
cocoh = zeros(3,3,numel(k1));
quadcoh = zeros(3,3,numel(k1));
% for u, v and w
% 1: u-direction
% 2: v-direction
% 3: w-direction
for ii=1:3,
    % single point spectra
    for jj=1:3,
        FM1=squeeze(trapz(k3_log,trapz(k2_log,PHI(:,:,:,ii,ii),2),3));
        FM2=squeeze(trapz(k3_log,trapz(k2_log,PHI(:,:,:,jj,jj),2),3));
%         crossPSD = real(trapz(k3_log,trapz(k2_log,PHI(:,:,:,ii,jj).*exp(-1i.*(k2.*d2+k3.*d3)),3),2));
%         A = crossPSD(end-round(N1/2)+1:end);
%         B = sqrt(FM1(end-round(N1/2)+1:end).*FM2(end-round(N1/2)+1:end));
%         coh(ii,jj,:) = interp1(k11,A./B,k1); % cross-coherence
        
        crossPSD = trapz(k3_log,trapz(k2_log,PHI(:,:,:,ii,jj).*exp(-1i.*(k2.*d2+k3.*d3)),3),2);
        coPSD = real(crossPSD(end-round(N1/2)+1:end));
        quadPSD = imag(crossPSD(end-round(N1/2)+1:end));
        B = sqrt(FM1(end-round(N1/2)+1:end).*FM2(end-round(N1/2)+1:end));
        
        cocoh(ii,jj,:) = interp1(k11,coPSD./B,k1); % cross-coherence
        quadcoh(ii,jj,:) = interp1(k11,quadPSD./B,k1); % cross-coherence
    end
end


% test = ppval(pp,k1);
% plot(k1,test,k1,coh0);xlim([0,1])


end

