function a=forNorm_CPU(ji,jt,jb,jl,jr,dxi,dxt,dxb,dxl,dxr,f0)
%calculates the norm

aa=f0+ji.*dxi+jt.*dxt+jb.*dxb+jl.*dxl+jr.*dxr;
a=abs(aa).^2;









