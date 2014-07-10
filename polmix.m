function [sih,siv] = polmix(tscat,sih,siv)

%   calculates the polarization mixing of the interface reflectivities
%       of each layer (taking into account the first order scattering)
%
%   [sih,siv] = polmix(tscat,sih,siv)
%       sih:   interface reflectivity at h-pol
%       siv:   interface reflectivity at v-pol
%       tscat: tau scat
%
%   Version history:
%      1.0    wi 14.10.97
%      1.1    wi  4.11.97  bug fix (layer numbering problem)
%
%   Uses: -
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

tscat = [tscat;1];

smean = 0.5 .* (sih + siv);
deltas = 0.5 .* tscat .* (sih - siv);
sih = smean + deltas;
siv = smean - deltas;
