function [ji,jt,jb,jl,jr,f0]=elementJacobian_CPU(hi,ht,hb,hl,hr,kt,kb,kl,kr,...
    zt,zb,zl,zr,Sy,Dbound,ki,zi,hOld,inflow,dt,kRi,itDo,ibDo,ilDo,irDo)
% creates jacobian for 2d depth averaged groundwater equation
%% PROGRAM INFORMATION
% Created by: Daniel Erdal (daniel.erdal@uni-tuebingen.de)
% Creations date:  March./April 2015

% minimum water table
% mwt=1e-3; 

mwt=1e-6; 

% Flow depth
flow_depth = 1*mwt; % No flow depth for cells with very little gw depth

% eps for creating finite differences 
e=1e-8;

% Inflow from the top cell==============================
at0=itDo .* kt.*((hi+ht)./2-zt).*(ht-hi);
atT=itDo .* kt.*((hi+ht+e)./2-zt).*(ht+e-hi);
atI=itDo .* kt.*((hi+e+ht)./2-zt).*(ht-(hi+e));
ix=(hi+ht)./2-zt<mwt & itDo;
if any(ix)
    % If the water table is too low that (h-z) is fixed
    at0(ix)=kt(ix).*flow_depth.*(ht(ix)-hi(ix));
    atT(ix)=kt(ix).*flow_depth.*(ht(ix)+e-hi(ix));
    atI(ix)=kt(ix).*flow_depth.*(ht(ix)-(hi(ix)+e));
end

% Inflow from the bottom cell==============================
ab0=ibDo .* kb.*((hi+hb)./2-zb).*(hb-hi);
abB=ibDo .* kb.*((hi+hb+e)./2-zb).*(hb+e-hi);
abI=ibDo .* kb.*((hi+e+hb)./2-zb).*(hb-(hi+e));
ix=(hi+hb)./2-zb<mwt & ibDo;

if any(ix)
    ab0(ix)=kb(ix).*flow_depth.*(hb(ix)-hi(ix));
    abB(ix)=kb(ix).*flow_depth.*(hb(ix)+e-hi(ix));
    abI(ix)=kb(ix).*flow_depth.*(hb(ix)-(hi(ix)+e));
end

% Inflow from the left cell==============================
al0=ilDo .* kl.*((hi+hl)./2-zl).*(hl-hi);
alL=ilDo .* kl.*((hi+hl+e)./2-zl).*(hl+e-hi);
alI=ilDo .* kl.*((hi+e+hl)./2-zl).*(hl-(hi+e));

ix=(hi+hl)./2-zl<mwt & ilDo;

if any(ix)
    al0(ix)=kl(ix).*flow_depth.*(hl(ix)-hi(ix));
    alL(ix)=kl(ix).*flow_depth.*(hl(ix)+e-hi(ix));
    alI(ix)=kl(ix).*flow_depth.*(hl(ix)-(hi(ix)+e));
end

% Inflow from the right cell==============================
ar0=irDo .* kr.*((hi+hr)./2-zr).*(hr-hi);
arR=irDo .* kr.*((hi+hr+e)./2-zr).*(hr+e-hi);
arI=irDo .* kr.*((hi+e+hr)./2-zr).*(hr-(hi+e));

ix = (hi+hr)./2-zr<mwt & irDo;
if any(ix)
    ar0(ix)=kr(ix).*flow_depth.*(hr(ix)-hi(ix));
    arR(ix)=kr(ix).*flow_depth.*(hr(ix)+e-hi(ix));
    arI(ix)=kr(ix).*flow_depth.*(hr(ix)-(hi(ix)+e));    
end

% Inflow from the Dirichlet boundaries==============================
b0=zeros(size(ki)); bI=b0;
ix = ~isnan(Dbound);
if any(ix)
    % OBS: ONLY CONSIDERING ONE DIRICLET PER CELL.....  
    % b0(ix)=ki(ix).*((Dbound(ix)+hi(ix))/2-zi(ix)).*(Dbound(ix)-hi(ix));
    % bI(ix)=ki(ix).*((Dbound(ix)+hi(ix)+e)/2-zi(ix)).*(Dbound(ix)-(hi(ix)+e));     
    
    if min(Dbound-zi)<0
        % error('Please enter a dirchlet value larger than the bottom elevation.')
    end

    % ix2 = ((Dbound + hi)/2 - zi) < 0 & ix ;
    ix2 = ix; % All cells at dirchlet B.C
    if any(ix2)
        % % Not allowing outflow of dirchlet cells 
        at0(ix2) = max(at0(ix2),0);
        ab0(ix2) = max(ab0(ix2),0);
        al0(ix2) = max(al0(ix2),0);
        ar0(ix2) = max(ar0(ix2),0);

        % Other values
        atI(ix2) = max(atI(ix2),0);
        abI(ix2) = max(abI(ix2),0);
        alI(ix2) = max(alI(ix2),0);
        arI(ix2) = max(arI(ix2),0);

        atT(ix2) = max(atT(ix2),0);
        abB(ix2) = max(abB(ix2),0);
        alL(ix2) = max(alL(ix2),0);
        arR(ix2) = max(arR(ix2),0);

        b0(ix2)= -(at0(ix2) + ab0(ix2) + al0(ix2) + ar0(ix2)); % Inflows = Outflows   
        ksat_factor = 10;
        flow_depth_dirch = b0(ix2)./(ksat_factor*ki(ix2).*(Dbound(ix2)-hi(ix2))); % Equivalent Depth
        idx = Dbound(ix2) >= hi(ix2); % Linearize around mwt
        flow_depth_dirch(idx) = mwt;
        % Derivative
        bI(ix2)=ksat_factor*ki(ix2).*flow_depth_dirch.*(Dbound(ix2)-(hi(ix2)+e));                
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