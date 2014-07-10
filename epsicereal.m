function z = epsicereal(TK)
% Function for calculating the real relative permittivity 
% of pure ice in the microwave region, according to 
% Matzler, C.(ed),"Thermal Microwave Radiation - Applications for Remote Sensing", 
% IET, London, UK (2006), Chapter 5.
% Input:
% TK = temperature (K), range 240 to 273.15

	z = 3.1884 + 9.1e-4*(TK-273);
