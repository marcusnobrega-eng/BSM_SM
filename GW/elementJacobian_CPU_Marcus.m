function [ji,jt,jb,jl,jr,f0]=elementJacobian_CPU(hi,ht,hb,hl,hr,kt,kb,kl,kr,...
    zt,zb,zl,zr,Sy,Dbound,ki,zi,hOld,inflow,dt,kRi,itDo,ibDo,ilDo,irDo,min_depth)
% creates jacobian for 2d depth averaged groundwater equation
%% PROGRAM INFORMATION
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% Creations date:  March./April 2015

% minimum water table
mwt=min_depth; 
% eps for creating finite differences 
e=1e-6;

% Inflow from the top cell==============================
at0=itDo .* kt.*((hi+ht)./2-zt).*(ht-hi);
atT=itDo .* kt.*((hi+ht+e)./2-zt).*(ht+e-hi);
atI=itDo .* kt.*((hi+e+ht)./2-zt).*(ht-(hi+e));
ix=(hi+ht)./2-zt<mwt & itDo;
if any(ix)
    % If the water table too low that (h-z) is fixed
    at0(ix)=kt(ix).*mwt.*(ht(ix)-hi(ix));
    atT(ix)=kt(ix).*mwt.*(ht(ix)+e-hi(ix));
    atI(ix)=kt(ix).*mwt.*(ht(ix)-(hi(ix)+e));
end

% Inflow from the bottom cell==============================
ab0=ibDo .* kb.*((hi+hb)./2-zb).*(hb-hi);
abB=ibDo .* kb.*((hi+hb+e)./2-zb).*(hb+e-hi);
abI=ibDo .* kb.*((hi+e+hb)./2-zb).*(hb-(hi+e));
ix=(hi+hb)./2-zb<mwt & ibDo;

if any(ix)
    ab0(ix)=0*kb(ix).*mwt.*(hb(ix)-hi(ix));
    abB(ix)=0*kb(ix).*mwt.*(hb(ix)+e-hi(ix));
    abI(ix)=0*kb(ix).*mwt.*(hb(ix)-(hi(ix)+e));
end

% Inflow from the left cell==============================
al0=ilDo .* kl.*((hi+hl)./2-zl).*(hl-hi);
alL=ilDo .* kl.*((hi+hl+e)./2-zl).*(hl+e-hi);
alI=ilDo .* kl.*((hi+e+hl)./2-zl).*(hl-(hi+e));

ix=(hi+hl)./2-zl<mwt & ilDo;

if any(ix)
    al0(ix)=0*kl(ix).*mwt.*(hl(ix)-hi(ix));
    alL(ix)=0*kl(ix).*mwt.*(hl(ix)+e-hi(ix));
    alI(ix)=0*kl(ix).*mwt.*(hl(ix)-(hi(ix)+e));
end

% Inflow from the right cell==============================
ar0=irDo .* kr.*((hi+hr)./2-zr).*(hr-hi);
arR=irDo .* kr.*((hi+hr+e)./2-zr).*(hr+e-hi);
arI=irDo .* kr.*((hi+e+hr)./2-zr).*(hr-(hi+e));

ix = (hi+hr)./2-zr<mwt & irDo;
if any(ix)
    ar0(ix)=0*kr(ix).*mwt.*(hr(ix)-hi(ix));
    arR(ix)=0*kr(ix).*mwt.*(hr(ix)+e-hi(ix));
    arI(ix)=0*kr(ix).*mwt.*(hr(ix)-(hi(ix)+e));
    
end

% Inflow from the Dirichlet boundaries==============================
b0=zeros(size(ki)); bI=b0;
ix = ~isnan(Dbound);
if any(ix)
    % OBS: ONLY CONSIDERING ONE DIRICLET PER CELL.....  
    b0(ix)=ki(ix).*((Dbound(ix)+hi(ix))/2-zi(ix)).*(Dbound(ix)-hi(ix));
    bI(ix)=ki(ix).*((Dbound(ix)+hi(ix)+e)/2-zi(ix)).*(Dbound(ix)-(hi(ix)+e));     
       
    ix2=Dbound-zi<mwt;
    
    if any(ix2)
        b0(ix2)=ki(ix2).*mwt.*(Dbound(ix2)-hi(ix2));
        bI(ix2)=ki(ix2).*mwt.*(Dbound(ix2)-(hi(ix2)+e));
    end
end

% Newton-Rahpson functions
f0=Sy.*(hi-hOld)-dt.*(at0+ab0+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fI=Sy.*(hi+e-hOld)-dt.*(atI+abI+alI+arI+bI)-dt.*inflow+dt.*kRi.*(hi+e);

fT=Sy.*(hi-hOld)-dt.*(atT+ab0+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fB=Sy.*(hi-hOld)-dt.*(at0+abB+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fL=Sy.*(hi-hOld)-dt.*(at0+ab0+alL+ar0+b0)-dt.*inflow+dt.*kRi.*hi;
fR=Sy.*(hi-hOld)-dt.*(at0+ab0+al0+arR+b0)-dt.*inflow+dt.*kRi.*hi;

% Construct the Jacobian entries 
% main diagonal
ji=(fI-f0)./e;
% first positive off diagonal
jt=(fT-f0)./e;
% first negative off diagonal
jb=(fB-f0)./e;
% second negative off diagonal
jl=(fL-f0)./e;
% second positive off diagon
jr=(fR-f0)./e;





