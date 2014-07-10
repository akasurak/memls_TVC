% tehon heijastuskerroin, horisontaalipolarisaatio
% kaava 6.7
% 5.12.89 jpk; 17.3.1993 JP

function Gammah = gammah(Epss,theta)
% kompleksinen dielektrisyysvakio Epss , kulma radiaaneina (asteina)

theta_rad = theta./180.*pi;
%theta_rad = theta;
cosTheta = cos(theta_rad);
neliojuuri = sqrt(Epss-sin(theta_rad).^2);
Gammah  = (abs((cosTheta-neliojuuri)./(cosTheta+neliojuuri))).^2;

