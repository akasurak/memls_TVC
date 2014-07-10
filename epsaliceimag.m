function z = epsaliceimag(fGHz,TK,Sppt)
% MATLAB Function for calculating the imaginary part of the
% relative permittivity of impure ice in the microwave region, 
% according to C. Matzler, IET Book (2006), chapter on ice
% Input:
% TK = temperature (K), range 20 to 273.15
% fGHz = frequency in GHz, range 0.01 to 3000
% Sppt = salinity in parts per thousand	
 B1 = 0.0207;
 B2 = 1.16e-11;
 b = 335;
 deltabeta = exp(-9.963 + 0.0372.*(TK-273));
 betam = (B1./TK).* ( exp(b./TK)./ ((exp(b./TK)-1).^2) ) + B2*fGHz.^2;
 beta = betam + deltabeta;
 theta = 300./TK - 1;
 alfa = (0.00504 + 0.0062*theta).*exp(-22.1*theta);
 zpure=alfa./fGHz + beta.*fGHz;
 g0=1866*(exp(-0.317*fGHz));
 g1=72.2+6.02*fGHz;
 dinv=0.013*(g0+g1.*abs(273.16-TK));
 dz=Sppt./dinv;
 z=zpure+dz;