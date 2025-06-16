function [f0]=elementf0_CPU(hi,ht,hb,hl,hr,kt,kb,kl,kr,zt,zb,zl,zr,Sy,Dbound,...
    ki,zi,hOld,inflow,dt,kRi,itDo,ibDo,ilDo,irDo)
%writes 2d depth averaged groundwater equation as 0=f0

% minimum water table
mwt=1e-6; 

% flow depth
flow_depth = mwt; % No flow depth for cells with very little gw depth

% compute fluxes across the four cell interfaces
at0=itDo .* kt.*((hi+ht)./2-zt).*(ht-hi); 
ix=(hi+ht)./2-zt<mwt & itDo;
at0(ix)=kt(ix).*flow_depth.*(ht(ix)-hi(ix)); 

ab0=ibDo .* kb.*((hi+hb)./2-zb).*(hb-hi); 
ix=(hi+hb)./2-zb<mwt & ibDo;
ab0(ix)=kb(ix).*flow_depth.*(hb(ix)-hi(ix));  

al0=ilDo .* kl.*((hi+hl)./2-zl).*(hl-hi);   
ix=(hi+hl)./2-zl<mwt & ilDo;
al0(ix)=kl(ix).*flow_depth.*(hl(ix)-hi(ix));   

ar0=irDo .* kr.*((hi+hr)./2-zr).*(hr-hi); 
ix=(hi+hr)./2-zr<mwt & irDo;
ar0(ix)=kr(ix).*flow_depth.*(hr(ix)-hi(ix));  

% compute fluxes across Dirichlet boundaries
b0=zeros(size(ki));
ix=~isnan(Dbound);
if any(ix)
    % OBS: ONLY CONSIDERING ONE DIRICLET PER CELL.....
    b0(ix)=ki(ix).*((Dbound(ix)+hi(ix))/2-zi(ix)).*(Dbound(ix)-hi(ix));
    ix2=Dbound-zi<mwt;
    b0(ix2)=ki(ix2).*flow_depth.*(Dbound(ix2)-hi(ix2));       
   
end

% calculate f0
f0=Sy.*(hi-hOld)-dt.*(at0+ab0+al0+ar0+b0)-dt.*inflow+dt.*kRi.*hi;