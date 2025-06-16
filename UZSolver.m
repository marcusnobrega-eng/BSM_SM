classdef UZSolver < handle
    %UZSolver Numerical solver for the 1d Richards Equation
    
    properties
        eps = 1.0e-7;           % convergence criterion for iteration
        epsh = 1.0e-7;          % Delta h for the finite difference of the retention function [m]
        jacSingBreak = false;   % if the jacobian is ill conditioned/ singular should the program break out
        convBreak1 = 300;       % Convergence check 1: reduce the time step
        convBreak2 = 1200;      % Convergence check 2: break out
        doLineSearch = true;    % do line search?
        MaxIter = inf;          % maximum number of iterations
    end
    
    methods
        %% constructor
        function obj = UZSolver()
        end
        % solve the Richards Equation
        function [h,R,hu] = solveUZ(obj,tStart,tEnd,source,dt,UZmodel,hGW)
            
            % initial time
            t = tStart;
            % counter
            zTot=0;
            
            % get often used variables
            parameter = UZmodel.parameter;
            boundaries = UZmodel.boundaries;
            poro = parameter.poro;
            trans = boundaries.trans;
            roots = UZmodel.roots;
            ss = UZmodel.ss+source;
            
            % times that need to be hit exactly
            mustDos=sort([boundaries.changeT(2:end),boundaries.changeB(2:end),...
                boundaries.changeTr(2:end),tEnd]);
            a=mustDos(2:end)-mustDos(1:end-1);
            a= a~=0;
            mustDos=mustDos(a);
            if isempty(mustDos)
                mustDos=tEnd;
            end
            
            % get old state vector
            h_old = UZmodel.h;
            headinter = UZmodel.h;
            
            % set the grid
            delzNN=UZmodel.delzNN;
            delzCS=UZmodel.delzCS;
            
            % calculate the water content
            [sat, se]=calcSat(UZmodel.h,parameter);
            wc_new=poro.*sat;
            nZ = UZmodel.nZ;
            nz = nZ+2;
            
            % calculate the harmonic average of the saturated
            % permeabilities
            parameter.getHarmonicAverages();
            
            %% --------------------------Time loop --------------------------------
            while t < tEnd
                
                % update boundary conditions
                boundaries.updateBoundaries(t);
                if boundaries.noFluxB
                    boundaries.currHeadB = UZmodel.h(1)+UZmodel.delzCS(1)/2;
                end
                if ~boundaries.doFreeDrain
                    parameter.calculateBoundaryKR(boundaries.currHeadB);
                end
                
                % calculate the relative permeability
                old.kr = newtonKrel(headinter,se,nz,parameter,boundaries);
                
                % Check for the maximum possible time step size
                a=mustDos-t;
                a(a<=0)=inf;
                dt=min(dt,min(a));
                
                % Set the water content at the old timestep
                wc_old=wc_new;
                old.wc=wc_old;
                
                %% ----------------------- Start of Nonlinear Iteration --------
                % Initialize a counter
                zaehler1=0;
                zaehler2=0;
                
                %-------------------------------------------------------
                % While loop, do as many Newton-Raphson steps until convercence 
                % is reached. Check every solver.convBreak1 steps if convergance 
                % problems and break out if solver.convBreak2 is reached.
                % The Newton-Raphson Iterations are done to find the minimum for
                % (S Storage h^n+1-h^n)/dt + nf S^n+1-S^n/dt) - div flux^n+1) =
                % f(h^n+1).

                % Initialize some variables
                rhs=1;
                rhsFirst=1;
                incrOrg=1;
                old.rhs=[];
                doBreak=false;
                
                % Keep iterating until the defect between current and last iteration
                % step is small enough. Second break of criteria is present only 
                % avoid a crash when reaching steady state (hence no defect decrease)
                while (sqrt(sum(rhs.^2))/sqrt(sum(rhsFirst.^2))>obj.eps && ...
                        max(abs(incrOrg))>1e-5) || any(isnan(incrOrg)) 
                    
                    %%  Do the iteration steps
        
                    % Calculate the Jacobian and the right hand side for the 
                    % Newton iterations
                    [J, rhs]=newtonJacobian(headinter,old,wc_old,h_old,parameter,...
                        dt,delzCS,delzNN,boundaries,obj.epsh,nz,roots,ss);
                    
                    % solve for the increment
                    incr=J\(-rhs);
                    incrOrg=incr;
                    
                    [~, msgid] = lastwarn;
                    lastwarn('')
                    
                    if any(isinf(incr)) || any(isnan(incr)) || ...
                            strcmp(msgid,'MATLAB:singularMatrix') || ...
                            strcmp(msgid,'MATLAB:nearlySingularMatrix')
                        % Jacobian is singular
                        if obj.jacSingBreak
                           % Break out!
                            doBreak=true;
                            headinter = ones(nZ,1)*inf;

                            disp('---------------Jacobian singular, breaking out!')
                            break

                        else               
                            % restart (i.e. make the system go to restart)
                            incr=zeros(nZ,1);
                            rhs=ones(nZ,1)*1e-20;
                            zaehler1=obj.convBreak1-1;
                        end
                    end

                    % set the first defect to the current if first iteration step
                    if zaehler1 == 0; rhsFirst=rhs; end
                    
                    %% Line search
                    % decrease the increment if the solution deteriorates
                    % implemeted after the numerics script of O. Ippisch p.140
                    if obj.doLineSearch           
                        % initialize alpha
                        alpha=1;

                        % original norm
                        normf0=norm(rhs,inf);

                        % alpha=1 norm
                        [sat,se]=calcSat((headinter+incr),parameter);
                        wc_new=sat.*poro;
                        storeX=store((headinter+incr),wc_new,wc_old,parameter,dt,...
                            delzCS,h_old);
                        kr = newtonKrel((headinter+incr),se,nz,parameter,boundaries);
                        sour = sourceRichards((headinter+incr),roots,trans,delzCS)...
                            +ss;
                        rhs=richardsNewton((headinter+incr),kr,storeX,delzNN,...
                            parameter,boundaries,sour);            
                        normf1=norm(rhs,inf);

                        % line search loop
                        while normf1 >= normf0 && alpha >= 2*1/10000 && normf1>1e-5
                            alpha=alpha*0.5;
                            [sat,se]=calcSat((headinter+alpha*incr),parameter);
                            wc_new=sat.*poro;
                            storeX=store((headinter+alpha*incr),wc_new,wc_old,...
                                parameter,dt,delzCS,h_old);
                            kr = newtonKrel((headinter+alpha*incr),se,nz,parameter,...
                                boundaries);
                            sour = sourceRichards((headinter+alpha*incr),roots,...
                                trans,delzCS)+ss;
                            rhs=richardsNewton((headinter+alpha*incr),kr,storeX,...
                                delzNN,parameter,boundaries,sour);               
                            normf1=norm(rhs,inf);
                        end

                        % decrease the increment
                        incr=alpha*incr;
                    else
                        % if no line search, just calculate the defect
                        [sat,se]=calcSat((headinter+incr),parameter);
                        wc_new=sat.*poro;
                        storeX=store((headinter+incr),wc_new,wc_old,parameter,dt,...
                            delzCS,h_old);
                        kr = newtonKrel((headinter+incr),se,nz,parameter,boundaries);
                        sour = sourceRichards((headinter+incr),roots,trans,delzCS)...
                            +ss;
                        rhs=richardsNewton((headinter+incr),kr,storeX,delzNN,...
                            parameter,boundaries,sour);
                    end
                    
                    %% back to the original solver        
              
                    % Save same variables for later re-use
                    old.rhs=rhs;
                    old.kr=kr;
                    old.store=storeX;

                    % Update the head
                    headinter=headinter+incr;

                    % Increase the counter by one
                    zaehler1=zaehler1+1;

                    % Check all X iteration for convergence problem. If such, 
                    % decrease timestep size and restart the timestep  

                    if(mod(zaehler1,obj.convBreak1)==0)

                        % reduce with 60%
                        dt=max(UZmodel.dtMin,0.4*dt);

                        % reset the counters
                        zaehler2=zaehler2+zaehler1;
                        zaehler1=0;

                        if zaehler2 >= obj.convBreak2
                            % Quite and bail out  
                            headinter = ones(nZ,1)*inf;
                            disp('---------------Conv. probl., breaking out!')

                            doBreak=true;
                            break
                        end

                        % otherwise restart and restart
                        headinter=h_old;
                        [sat,se]=calcSat(headinter,parameter);
                        clear rhs kr storeX
                        rhs=1000000;
                        wc_new=sat.*poro;
                        wc_old=wc_new;
                        old.wc=wc_new;
                        old.store=store(headinter,wc_new,wc_old,parameter,dt,...
                            delzCS,h_old);
                        old.kr = newtonKrel(headinter,se,nz,parameter,boundaries);
                        old.rhs=[];
                    end

                    %-------------- end of Newton iteration ---------------
                end
                %%
    
                if doBreak
                    break
                end

                % Set time to the new time
                t=t+dt;

                zTot=zTot+zaehler1+zaehler2;
                if zTot>obj.MaxIter && t<tEnd
                    headinter = ones(nZ,1)*inf;
                    disp('---------------Too many iterations, breaking out')
                    break
                end
                
                %Check if the next timestep should be different
                if UZmodel.doAdaptiveTime
                    dt=timecheck(dt,zaehler1,sat-wc_old./parameter.poro);
                    % and check the boundaries
                    dt=min(UZmodel.dtMax,dt);
                    dt=max(UZmodel.dtMin,dt);
                end

                % Set old time to the actual time to start the new time.
                h_old=headinter;

                %----------------------- End of time loop -----------------------------
            end
            
            %% calculate groundwater table position and recharge
            full = false;
            if any(isinf(headinter))
                R = 0;
                h = headinter;
                hu = nan;
            else
                headinter(headinter>UZmodel.H-UZmodel.z) = UZmodel.H - ...
                    UZmodel.z(headinter>UZmodel.H-UZmodel.z);
                if headinter(1)>0
                    hi = find(headinter<0,1,'first');
                    if isempty(hi)
                        full = true;
                    else
                        hs = headinter(hi-1);
                        zs = UZmodel.z(hi-1);
                        hu = headinter(hi);
                        zu = UZmodel.z(hi);
                    end
                else
                    hs = boundaries.currHeadB;
                    zs = UZmodel.z(1)-delzCS(1)/2;
                    hu = headinter(1);
                    zu = UZmodel.z(1);
                end
                if ~full
                    hNeu = zs+(zu-zs)*(0-hs)/(hu-hs);
                else
                    hNeu = UZmodel.z(end)+headinter(end);
                end

                R = (hNeu-hGW);
                hu = hNeu;
                h=headinter;

            end
        end
    end
end

