function result = loadsoil(filename)
y=load(filename);
result.t=y(1,1); %Enter Soil temperature in Kelvin
result.mv=y(1,2);    %Enter soil moisture m3/m3
result.sig=y(1,3)./100; %Enter soil roughness (sigSoil is in meters)