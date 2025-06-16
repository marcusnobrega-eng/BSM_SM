classdef GWSolver < handle
    %GWSolver Numerical solver for the 2d depth averaged groundwater
    %equation
    
    properties
        doLineSearch = true;            % use LineSearch
        how = 1;                        % how to solve linear equations ...
                                        % (1: Matlab; 2: Biconjugate gradient; ...
                                        % 3: gmres or GPU-Jacobi)
        doExplicitPredcond = false;     % precondition with explicit Euler?
        EulerMax = 1e-5;                % max dh for explicit methods
        cErr = 1e-7;                    % convergence tolerance        
        maxJacobi = 1000;               % maximum number of iteration steps
                                        % for Jacobi method
        mwt=1e-6;                      % minimum water table to allow flow
                                        % out of cell
    end
    
    methods
        %% constructor
        function obj = GWSolver()
        end
        %% solve the groundwater equation
        function h = solveGW(obj,tStart,tEnd,dt,GWmodel,R)
            
            %% Sort out indexes:
            % Welcome to an indexing nightmare; this takes into account also
            % non-rectangular systems.

            % Total index (index of each active cell)
            tmp=reshape(1:GWmodel.nx*GWmodel.ny,GWmodel.ny,GWmodel.nx);
            tmp(GWmodel.boundaries.mask==0)=nan;
            tmp=reshape(tmp,GWmodel.nx*GWmodel.ny,1);
            tmp(isnan(tmp))=[];
            indexTot=tmp;

            % Number of active cells:
            nrCells=length(indexTot);

            %===Neighbour index (top, bottom, left, right)===

            tmp=nan(GWmodel.ny,GWmodel.nx);
            tmp(GWmodel.boundaries.mask==1)=1:nrCells;
            tmp2=nan(GWmodel.ny+2,GWmodel.nx+2);
            tmp2(2:end-1,2:end-1)=tmp;

            % Actual cell index
            ime(1:nrCells,1)=1:nrCells;

            % identify cells with low GW
            lowGWmask = GWmodel.h(indexTot) - GWmodel.z0(indexTot) < obj.mwt;
            lowGWid = indexTot(lowGWmask);
            lowGWh =  GWmodel.h(indexTot);

            % Right hand neighbour
            tmp=tmp2(2:end-1,3:end);
            tmp(GWmodel.boundaries.mask==0)=nan;
            ir=tmp(indexTot);       % the neighbours index
            irDo=true(size(ir));    % index as boolean (used to create the Jacobian)
            irDo(isnan(ir))=false;
            ir(isnan(ir))=1;

            % Left hand neighbour
            tmp=tmp2(2:end-1,1:end-2);
            tmp(GWmodel.boundaries.mask==0)=nan;
            il=tmp(indexTot);
            ilDo=true(size(il));
            ilDo(isnan(il))=false;
            il(isnan(il))=1;

            % Top neighbour
            tmp=tmp2(3:end,2:end-1);
            tmp(GWmodel.boundaries.mask==0)=nan;
            it=tmp(indexTot);
            itDo=true(size(it));
            itDo(isnan(it))=false;
            it(isnan(it))=1;

            % Bottom neighbour
            tmp=tmp2(1:end-2,2:end-1);
            tmp(GWmodel.boundaries.mask==0)=nan;
            ib=tmp(indexTot);
            ibDo=true(size(ib));
            ibDo(isnan(ib))=false;
            ib(isnan(ib))=1;

            % if a cell has a low GW table, do not allow outflow
            % (i.e. flow to a neighbor with lower h)

            if ~isempty(lowGWid)
                % disp(lowGWid)
                cellh = lowGWh(lowGWid);
                righth = lowGWh(ir(lowGWid));
                lefth = lowGWh(il(lowGWid));
                toph = lowGWh(it(lowGWid));
                bottomh = lowGWh(ib(lowGWid));
                
                % if the right hand neighbour has a lower potential,
                % no outflow to that cell
                irDo(lowGWid(cellh >= righth)) =false;
                % set the left flux for the right cell to false
                % to ensure consistency
                ilDo(ir(lowGWid(cellh >= righth))) = false;

                ilDo(lowGWid(cellh >= lefth)) =false;
                irDo(il(lowGWid(cellh >= lefth))) = false;

                itDo(lowGWid(cellh >= toph)) =false;
                ibDo(it(lowGWid(cellh >= toph))) = false;

                ibDo(lowGWid(cellh >= bottomh)) =false;
                itDo(ib(lowGWid(cellh >= bottomh))) = false;
            end

            % Put the indexes on the GPU
            % ime=(ime);
            % ir=(ir);
            % il=(il);
            % it=(it);
            % ib=(ib);

            % What goes where in the Jacobian matrix
            iA1=([ime;ime(itDo);ime(ibDo);ime(ilDo);ime(irDo)]);
            iA2=([ime;it(itDo);ib(ibDo);il(ilDo);ir(irDo)]);

            % the index of the perimiter cells (boundary cells)
            iPerim=GWmodel.boundaries.perimI(indexTot);
            clear tmp tmp2
            
            %% Make sure that a boundary is only on a cell with a boundary!
            % (hence: delete from the boundary index cells that are on inner 
            % edges on a boundary)
            maskL=zeros(GWmodel.ny+2,GWmodel.nx+2);
            maskL(2:end-1,2:end-1)=GWmodel.boundaries.mask;
            b1=diff(maskL(:,2:end-1),1,1)~=0;
            b2=diff(maskL(2:end-1,:),1,2)~=0;
            clear maskL

            % find out how many boundaries each cell holds
            % (for the GPU implementation this is not used)
            perim=(b1(1:end-1,:)+b1(2:end,:)+b2(:,1:end-1)+b2(:,2:end)) .* ...
                GWmodel.boundaries.mask;
            iPerim(perim(indexTot)==0)=0;
            clear perim b1 b2
            
            %% Pre-process the river boundary
            % Asumption: river is on a finer (or same) grid as the problem 
            % AND the two grids are just different scales of each other 
            % (dx/riverDX=integer>=0)
            if GWmodel.rivers.doRivers
                r=GWmodel.dx/GWmodel.rivers.dx;
                if r==1
                    Ar=GWmodel.rivers.dx*GWmodel.rivers.dy*GWmodel.rivers.riverMask;
                else        
                    Ar=zeros(size(GWmodel.boundaries.mask));
                    % Loop through each actual grid cell
                    % (could be speeded up on the gpu.....)
                    for i=1:GWmodel.ny
                        for j=1:GWmodel.nx
                            Ar(i,j)=sum(sum(GWmodel.rivers.riverMask(1+(i-1)* ...
                                r:i*r,1+(j-1)*r:j*r)));
                        end
                    end
                    Ar=GWmodel.rivers.dx.*GWmodel.rivers.dy.*Ar;
                end
                Ar(GWmodel.boundaries.mask==0)=0;
                Ar(GWmodel.boundaries.perimI~=0)=0;

                KRiver=(GWmodel.rivers.K./GWmodel.rivers.d.*Ar./ ...
                    (GWmodel.dx*GWmodel.dy));
                Ar(Ar~=0)=1;
                GWmodel.rivers.riverMaskDomain=logical(Ar);
                clear Ar

                KRiverI=KRiver(indexTot);
            else
                KRiverI=0;
            end


            %% Wells: part 1
            % Get the index of the cells with the wells

            if GWmodel.wells.doWells
                wXY=zeros(GWmodel.ny,GWmodel.nx);
                for i=1:GWmodel.wells.nrWells
                    wXY(ceil(GWmodel.wells.y(i)/GWmodel.dy),...
                        ceil(GWmodel.wells.x(i)/GWmodel.dx))=i;
                end
            else
                wXY=0;
            end

            %% Conductivity (K) and bedrock elevation (Z)
            % Get the values at the interfaces between the cells

            GWmodel.parameter.K(GWmodel.boundaries.mask==0)=nan; % If we want to change K, is right here
            K=(nan(GWmodel.ny+2,GWmodel.nx+2));
            K(2:end-1,2:end-1)=GWmodel.parameter.K;
            % Call the GPU function for averaging the Ks
            dxF=zeros(GWmodel.ny,GWmodel.nx)+GWmodel.dx; 
            dyF=zeros(GWmodel.ny,GWmodel.nx)+GWmodel.dy; 
            [ki,kt,kb,kl,kr]=elementK_CPU(K(2:end-1,2:end-1),K(3:end,2:end-1),...
                K(1:end-2,2:end-1),K(2:end-1,1:end-2),K(2:end-1,3:end),dxF,dyF);
            clear K

            ki=ki(indexTot);    % Conductivity for the potential Diricelt boundary
            kt=kt(indexTot);    % Conductivity at top interface
            kb=kb(indexTot);    % Conductivity at bottom interface
            kl=kl(indexTot);    % Conductivity at left interface
            kr=kr(indexTot);    % Conductivity at right interface

            % Same as for K but for Z
            Z=(nan(GWmodel.ny+2,GWmodel.nx+2));
            Z(2:end-1,2:end-1)=GWmodel.z0;
            [zi,zt,zb,zl,zr]=elementZ_CPU(Z(2:end-1,2:end-1),Z(3:end,2:end-1),...
                Z(1:end-2,2:end-1),Z(2:end-1,1:end-2),Z(2:end-1,3:end));
            clear Z

            zi=zi(indexTot);
            zt=zt(indexTot);
            zb=zb(indexTot);
            zl=zl(indexTot);
            zr=zr(indexTot);

            % Porosity on the GPU
            Sy=(GWmodel.parameter.Sy(indexTot));
            
            % soil thickness
            lsurf = GWmodel.lsurf(indexTot);
            
            %% Setup for Transient Solution

            % Create the storage vector for the pressure heads
            hInit = GWmodel.h;

            % Old head = inital head
            if numel(hInit)==GWmodel.nx*GWmodel.ny
                hOld=(hInit(indexTot));
            else
                hOld=hInit;
            end
            
            % initialize solver parameters, counters etc
            zahler2=0; ti=tStart;
            forceQuit=false; 
            
            %% time loop
            
            while ti<tEnd
                [Q,boundDir,boundNeu] = GWmodel.boundaries.updateBoundaries(ti,...
                    iPerim,indexTot);
                if GWmodel.wells.doWells
                    W = GWmodel.wells.updateWells(ti,GWmodel.ny,GWmodel.nx,wXY);
                    Q=Q+W(indexTot);
                end
                if GWmodel.rivers.doRivers
                    hR = GWmodel.rivers.updateRivers(ti);
                    hRiver = hR+GWmodel.lsurf;
                    hRiver(KRiver==0) = 0;
                    Q = Q-KRiverI.*hRiver(indexTot);
                end
                Q = Q+R(indexTot);

                
                %% PRECONDITION
                if obj.doExplicitPredcond
                    % Solve the time step with explicit Euler and use the solution
                    % as initial quess for the nonlinear iterations.
                    hOldE=hOld;
                    inflow=(-Q + boundNeu);
                    tx=0;
                    c=1;
                    while tx<dt

                        [flowB,dtMin]=elementEulerExpl_CPU(hOldE(ime),hOldE(it),...
                            hOldE(ib),hOldE(il),hOldE(ir),kt,kb,kl,kr,zt,zb,zl,zr,...
                            boundDir,ki,zi,inflow,Sy,KRiverI,obj.EulerMax);
                        dtX=min(min(dtMin),dt-tx);
                        hOldE=hOldE+dtX.*flowB;
                        tx=tx+dtX;
                        c=c+1;

                        if c>10000
                            hOldE=hOld;
                            break
                        end

                    end
                    hItOld=(hOldE);
                 %   clear hOldE flowB dtMin
                else
                    hItOld=(hOld);
                end

                dtF=zeros(nrCells,1)+dt;
            
                        %% START THE SOLVER
            zahler3=0;   zMax=50;

            while 1
                
                dh=inf; zahler1=0; normf1=inf;
        
                %% Start of the nonlinear iterations
                hOld=(hOld);   % (assure it is on the GPU)
                while max(abs(dh))>1e-4 && zahler1<zMax && (normf1>1e-7 || ...
                        max(abs(dh))>1e-4)
                    
                    inflow=(-Q + boundNeu);
            
            
                  % Compute the entries of the Jacobian             
                    [ji,jt,jb,jl,jr,f0]=elementJacobian_CPU(hItOld(ime),...
                        hItOld(it),hItOld(ib),hItOld(il),hItOld(ir),kt,kb,kl,kr,...
                        zt,zb,zl,zr,Sy,boundDir,ki,zi,hOld,inflow,dtF,KRiverI,itDo,ibDo,ilDo,irDo);

                    if obj.how == 1
        %================= use the MATLAB backslash operator=======================

                        % Construct the Jacobian
                        J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));...
                            (jr(irDo))]);
                        % Solve the system
                        dh=J\-(f0);                

                    elseif obj.how == 2
        %================= use the iterative BICGSTAB METHOD=======================

                        % Construct the Jacobian
                        J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));...
                            (jr(irDo))]);
                        % Incomlete LU-factorization for preconditioning
                        [L,U]=ilu(J);
                        % Solve the system
                        [dh,out]=bicgstab(J,-(f0),obj.cErr,200,L,U);

                        % Looking like the system did not converge
                        if out~=0
                            % check if we actually have a problem
                            if max(abs(J*dh+f0))>1e-9
                                % We do: solve with \
                                disp('No covergence, solving using mldivide')
                                dh=J\-(f0);
                            end
                        end

                    elseif obj.how == 3
        %================= use the iterative JACOBI METHOD=========================                

                       % Initial guess (0 is here often best)
                       dx=zeros(nrCells,1);

                      c=1;


                      % Current setup: calculate the error in the same way as the
                      % BICGSTAB.... Rather costly but gives comparable results

                      % Norm of the right hand side
                      normf0=norm(f0,2);
                      while 1
                          % Compute the new dx
                          [dxnew]=solveEQ_CPU(ji,jt,jb,jl,jr,dx(it),dx(ib),...
                              dx(il),dx(ir),-f0);


                          % Calculate the norm (sum of a GPU-array is EXPENSIVE!)
                          eX=forNorm_CPU(ji,jt,jb,jl,jr,dxnew,dxnew(it),...
                              dxnew(ib),dxnew(il),dxnew(ir),f0);
                          err2=sqrt(sum(eX))/normf0; % norm(f0-J*dxnew)/norm(f0)


                          if err2 < obj.cErr
                              % Convergance: save iteration number and exit
                              cTot=cTot+c;
                              break

                          elseif c>obj.maxJacobi
                              % If no convergance, use \
                              disp('No covergence, solving using mldivide')
                              J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));...
                                  (jl(ilDo));(jr(irDo))]);
                              dxnew=J\-(f0);
                              break
                          end
                          % Iterate on.....
                          c=c+1;
                          dx=dxnew;
                      end
                      clear dx
                      dh=dxnew;
                      clear dxnew

                     elseif obj.how == 4
        %================= MIXTURE =======================
                    if zahler1<=1
                        % do matlab
                         % Construct the Jacobian
                        J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));...
                            (jr(irDo))]);
                        % Solve the system
                        dh=J\-(f0);                
                    else
        %                 % do bicgstab
        %                 
                        % Construct the Jacobian
                        J=sparse(iA1,iA2,[(ji);(jt(itDo));(jb(ibDo));(jl(ilDo));...
                            (jr(irDo))]);
                        % Incomlete LU-factorization for preconditioning
                        [L,U]=ilu(J);
                        % Solve the system
                        [dh,out]=bicgstab(J,-(f0),obj.cErr,200,L,U);

                        % Looking like the system did not converge
                        if out~=0
                            % check if we actually have a problem
                            if max(abs(J*dh+f0))>1e-9
                                % We do: solve with \
                                disp('No covergence, solving using mldivide')
                                dh=J\-(f0);
                            end
                        end


                    end
                    end
                 %   clear ji jt jb jl jr
        %==========================================================================
                %% DO LINESEARCH            

                    alpha=1;
                    if obj.doLineSearch

                        % original norm
                        if obj.how~=3                   
                           normf0=norm(f0,2);
                        end

                        % alpha=1 norm
                        hX=hItOld+dh;
                        f1_1=elementf0_CPU(hX(ime),hX(it),hX(ib),hX(il),hX(ir),...
                            kt,kb,kl,kr,zt,zb,zl,zr,Sy,boundDir,ki,zi,hOld,...
                            inflow,dtF,KRiverI,itDo,ibDo,ilDo,irDo);
                        normf1_1=norm(f1_1,2);
                        normf1=normf1_1;

                        % Iterate until the new soltuion has a smaller norm than
                        % the old one
                        while normf1 >= normf0 && alpha > 2*1/1000 && normf1>1e-6
                            alpha=alpha*0.5;
                            hX=hItOld+alpha*dh;
                            f1=elementf0_CPU(hX(ime),hX(it),hX(ib),hX(il),hX(ir),...
                                kt,kb,kl,kr,zt,zb,zl,zr,Sy,boundDir,ki,zi,hOld,...
                                inflow,dtF,KRiverI,itDo,ibDo,ilDo,irDo);
                            normf1=norm(f1,2);
                        end

                        if alpha <= 2*1/1000 %(max(abs(dh))*alpha) <= 0.01%
                            % hence, if nothing good was found
                            % take the largest or smallest step
                            if round(normf1/normf0) >= round(normf1_1/normf0)
                                normf1=normf1_1;
                                alpha=1;
                            end
                        end

                    else
                        % Just calculate the norm of the new soltion
                        % (if needed for judging convergance)
                        normf1=norm(f0,2);

                    end
                    %-----------------

                    % Update and go on
                   hNew=hItOld+alpha*dh;             
                  
                    zahler1=zahler1+1;
                    hItOld=hNew;                    
            
                end
                
                %% Check the solution        
        
                % Update the counter
                zahler2=zahler2+zahler1;

                % Check if we are ok or if exited for wrong reasons
                if zahler1==zMax || any(isnan(hNew))

                    % No convergence: reduce time step and retry (It does
                    % not work)
                    dt=dt/2;
                    warning('Model did not converge') % Edit Marcus
                    error('GW model not converging. Halving the time-step to check if it converges.')
                    zahler3=zahler3+1;
                    hItOld=(hOld);
                    % Forcing quit because this is not working
                    forceQuit=true;
                    break

                else
                    % Update the time
                    ti=ti+dt;
                    break
                end

                if zahler3>1
                    % after 1 times reduction of dt, quite
                    forceQuit=true;
                    break
                end
        
        
            end

            % if we broke out for real, then this is the end of the line: quit
            if forceQuit
%                 hNew = NaN*ones(size(hNew));
                % EDIT
                error('GW Model did not converge. Saving potentially wrong results.')
                ti=ti+dt;
                break
            end
            
            %% STORE AND UPDATE

            % Old is new and on we go.....
            hOld=hNew; 
    
            end
            
            h = hNew;
        
        end
    
    end

end
