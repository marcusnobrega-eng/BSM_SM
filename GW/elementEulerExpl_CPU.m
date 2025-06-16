function [flowB,dtMinX]=elementEulerExpl_CPU(hi,ht,hb,hl,hr,kt,kb,kl,kr,...
    zt,zb,zl,zr,Dbound,ki,zi,inflow,Sy,kRi,Emin)
%solve 2d depth averaged groundwater equation with explicit Euler scheme
error('not adjusted for low GW');
% minimum water table
mwt=0.01; 

% compute fluxes across the four cell interfaces
qt=kt.*((hi+ht)./2-zt).*(ht-hi);
ix=(hi+ht)./2-zt<mwt;
if any(ix)
    qt(ix)=kt(ix).*mwt.*(ht(ix)-hi(ix));
end

qb=kb.*((hi+hb)./2-zb).*(hb-hi); 
ix = (hi+hb)./2-zb<mwt;
if any(ix)
    qb(ix)=kb(ix).*mwt.*(hb(ix)-hi(ix));   
end

ql=kl.*((hi+hl)./2-zl).*(hl-hi); 
ix = (hi+hl)./2-zl<mwt;
if any(ix)
    ql(ix)=kl(ix).*mwt.*(hl(ix)-hi(ix));   
end

qr=kr.*((hi+hr)./2-zr).*(hr-hi);
ix = (hi+hr)./2-zr<mwt;
if any(ix)
    qr(ix)=kr(ix).*mwt.*(hr(ix)-hi(ix));  
end

% compute fluxes across Dirichlet boundaries
b0=zeros(size(hi));
ix = ~isnan(Dbound);
if any(ix)
    % OBS: ONLY CONSIDERING ONE DIRICLET PER CELL.....
    b0=ki.*((Dbound+hi)/2-zi).*(Dbound-hi);
    ix2 = Dbound-zi<mwt;
    if any(ix2)
        b0(ix2)=ki(ix2).*mwt.*(Dbound(ix2)-hi(ix2));  
    end
end

% calculate dh/dt
flowB=(qt+qb+ql+qr+b0+inflow-kRi.*hi)./Sy;

% calculte min. time step
dtMinX=Emin./abs(flowB);