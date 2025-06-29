function [zii,zit,zib,zil,zir]=elementZ_CPU(zi,zt,zb,zl,zr)
%calculate harmonic averages of the bottom elevation at cell interfaces

% Arithmetic avarage of the bedrock elevations
zit=(zt+zi)./2;   zit(isnan(zit))=0; 
zib=(zb+zi)./2;   zib(isnan(zib))=0; 
zil=(zl+zi)./2;   zil(isnan(zil))=0; 
zir=(zr+zi)./2;   zir(isnan(zir))=0; 

% for the diriclet boundaries
zii=zi;



