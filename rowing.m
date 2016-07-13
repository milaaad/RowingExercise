function result = rowing(N, amplitude, duration, varargin)
% Finds optimal motion of a rowing machine.  Task is to move from one
% static position to reach to a specific point in 1st phase, and go back to
% zero position from that specific point in phase 2. Duration of the phases
% is unknown, but total duration is fixed.

% author: Milad Zare m.zarei@csuohio.edu

% ________         _______________     ____________       ___ 
%|Flywhell|       |Spring & Damper|   |Machine Mass|     |Arm|
%|________|-------|_______________|---|____________|--->>|___|
%                      ______                     |
%---------------------|Spring|--------------------|
%                     |______|                 


% __________       _______________     ____       ____ 
%|m1, b1, c1|     |    k1 & b2    |   | m2 |     | m3 |
%|__________|-----|_______________|---|____|--->>|____|
%                      _______                     |
%---------------------|  k2   |--------------------|
%                     |_______|                 
% Inputs:
%		N			(scalar) Number of collocation points to use in each phase of the rowing cycle
%		amplitude   It defines the displacement of the machine
%       duration    The duration of phases is given
%       varargin{1}	It can be either off or on. 'on' adds muscle dynamics
%                   to the model and 'off' use linear equations to
%                   calculate the force
%       varargin{2} To use old results from previous run, varargin{2} needs
%                   to be defined. Old results from previous run can be used as an initial guess
%                   To activate this option, you need to use "result" as a 5th input. 
%                   i.e. rowing(100, 0.2, 2, 'off', 'result')
% Method: Direct Collocation
    clc
	tic
    global model

	% initializations
 	model.N          = N;                   % number of collocation nodes in each phase
    model.amplitude  = amplitude;           % amplitude of the desired motion (m)
    model.duration   = duration;            % duration of the motion (s)
    model.muscle_dyn = varargin{1};         % having muscle dynamics on or off
   	model.Ncon       = 6*(2*N-1)+7;         % 2N-1 dynamics constraints and 7 task constraints(periodicityof all states and 
                                            % controls except for flywheel position) and T1+T2=duration)
    
    % muscle parameters
    model.L          = 0.60;                % Length of forearm + hand(m)
    model.d          = 0.052;               % Moment arm at elbow(m) at 90 degrees (m)
    model.Tact       = 10e-3;               % Activation time(s)
    model.Tdeact     = 40e-3;               % Deactivation time(s)
    model.Lm0        = 0.288;               % Muscle-tendon length at 90 deg elbow angle (m)  
    
    % muscle properties of the biceps brachii
	model.Fmax      = 495;                  % (N) maximal isometric force (sum of both parts of the muscle)
	model.umax      = 0.033;                % (dimensionless) strain in the series elastic element at load of Fmax
	model.W         = 0.63;                 % (dimensionless) width parameter of the force-length relationship of the contractile element
	model.AHill     = 0.25;                 % (dimensionless) Hill parameter of the force-velocity relationship
	model.FVmax     = 1.5;                  % (dimensionless) maximal eccentric force
	model.Lceopt    = 0.1106;               % (m) optimal CE length (MUST BE AVERAGE OF THE TWO PARTS OF THE BICEPS)
	model.b         = 0.01;                 % (s/m) damping coefficient of damper parallel to the CE (normalized to Fmax)
   	model.SEELslack = 0.192;                % (m) slack length of the series elastic element
    model.PEELslack = 1;                    % (dimensionless) slack length of the parallel elastic element, divided by Lceopt

    % constants derived from the muscle properties
	model.Vmax  = 10*model.Lceopt;                                         % Maximum shortening velocity (m/s) is 10 fiber lengths per sec
	model.d1	= model.Vmax*model.AHill*(model.FVmax-1)/(model.AHill+1);  % parameter in the eccentric force-velocity equation
    
    % Since the model parameters are obtained from parameters
    % identification of a Concep-2 rowing machine, they have to be scaled
    % down to use in this model
    % So the model parameters should be multiplied by scale factor
    
    % the force scale factor(SF) = maximal pull force based on biceps
    % muscle strength / maximal force from Kleshnev paper
    % maximum force in Kleshnev data is 800 N
    FSF = ((model.d/model.L)*model.Fmax)/(800);     % force scale factor
    
    % The displacement scale factor(SD) = amplitude / maximum displacement
    % in Kleshnev test
    % maximum displacement in Kleshnev data is 1.35 m
    DSF = 0.2/(1.35);         % distance scale factor
    
    % Scale Factor = force scale factor/distance scale factor
    SF = FSF/DSF;
    
    % So stiffness parameters, linear damping stiffness and masses should
    % be multiplied by (SF)
    % The quadratic damping parameter should be multiplied by (FSF/DSF^2)
    
    % model parameters (obtained from parameters identification code)
    model.m1 = SF*(66.12/2);          % mass of the flywheel (kg)
    model.m2 = SF*0.0;              % mass of rowing machine (kg)
    model.m3 = 1.1285;              % efective mass of arm (kg)
    model.k1 = SF*(3046.07);        % spring stiffness of rewinding spring (N/m)
    model.k2 = SF*(31.98);          % spring stiffness of resistive spring (N/m)
    model.b1 = SF*(0);              % linear damping coefficient of flywheel (Ns/m)
    model.b2 = SF*(59.23);          % linear damping of shock cord (Ns/m)
    model.c1 = SF*(70.13)/(DSF);    % quadratic damping coefficient of flywheel (Ns^2/m^2)
    
    model.discretization = 'ME';    % it can be either Backward Euler(BE) or Midpoint Euler(ME)
   
    % state variable are x1, v, x2, v2, Lce, and a 
	% control variable is u
    % duration of phases is unknown, but the total duration is prescribed

	% if oldresult was provided, use it as initial guess, otherwise use zero initial guess
    if (nargin == 5)
        load('result')
		oldtime = result.time;         % sample times of old result, from 0 to 1
		h   = duration/(2*N-1);
        newtime = h*(0:2*N-1)';
        x1  = interp1(oldtime, result.x1, newtime);
        v1  = interp1(oldtime, result.v1, newtime);
        x2  = interp1(oldtime, result.x2, newtime);
        v2  = interp1(oldtime, result.v2, newtime);
        Lce = interp1(oldtime, result.Lce, newtime);
        a   = interp1(oldtime, result.a, newtime);
        u   = interp1(oldtime, result.u, newtime);
        T1  = result.t1;
        T2  = result.t2;
	else
		x1  = zeros(2*model.N,1);           % initial guess for flywheel position
        v1  = zeros(2*model.N,1);           % initial guess for flywheel velocity
        x2  = zeros(2*model.N,1);           % initial guess for rowing machine position
        v2  = zeros(2*model.N,1);           % initial guess for rowing machine velocity
        Lce = zeros(2*model.N,1);           % initial guess for flywheel velocity
        a   = zeros(2*model.N,1);           % initial guess for rowing machine position
        u   = zeros(2*model.N,1);           % initial guess for muscle neural control
        T1  = model.duration/2;             % initial guess for duration of phase 1
        T2  = model.duration/2;             % initial guess for duration of phase 2
    end

    % encode initial guess of unknowns into a long column vector X
    % there are N nodes in each phase. since there are 2 phases, totally there should be 2*N nodes
	X0 = [x1 ; v1; x2; v2; Lce; a; u; T1; T2];
    
    model.ix1  = (1:2*model.N);                     % index to elements in X where flywheel positions are stored
	model.iv1  =   2*model.N + (1:2*model.N);       % index to elements in X where velocities of flywheel are stored
    model.ix2  = 2*2*model.N + (1:2*model.N);       % index to elements in X where positions of rowing machine are stored
    model.iv2  = 3*2*model.N + (1:2*model.N);       % index to elements in X where velocities of rowing machine are stored
    model.iLce = 4*2*model.N + (1:2*model.N);       % index to elements in X where fiber lentghs are stored
    model.ia   = 5*2*model.N + (1:2*model.N);       % index to elements in X where velocities of activation are stored
    model.iu   = 6*2*model.N + (1:2*model.N);       % index to elements in X where controls u are stored
    model.iT1  = 7*2*model.N+1;                     % index to elements in X where duration of phase 1 is stored
    model.iT2  = 7*2*model.N+2;                     % index to elements in X where duration of phase 2 is stored  
	model.NX   = size(X0,1);                        % number of unknowns

    % checking the derivations that function calculates with the numerical derivations
    if (model.N == 10)
        checkderiv(X0); 
    end
    
    % determine Jacobian structure and number of nonzeros
    model.Jstructure = double(conjac ~= 0.0);
    
    % make bounds for the unknowns
    options.lb(model.ix1)    = -1.0;         % lower bound of flywheel position
    options.ub(model.ix1)    = 1.0;          % upper bound of flywheel position
    options.lb(model.ix1(1)) = 0;            % lower bound of flywheel position
    options.ub(model.ix1(1)) = 0;            % upper bound of flywheel position
    options.lb(model.iv1)    = 0;            % lower bound of flywheel velocity
    options.ub(model.iv1)    = 3.0;          % upper bound of flywheel velocity
    options.lb(model.ix2)    = 0;            % lower bound of handle position
    options.ub(model.ix2)    = amplitude;    % upper bound of handle position
    options.ub(model.ix2(1)) = 0;            % at t=0, handle must be at 0
    options.lb(model.ix2(N)) = amplitude;    % lower bound of handle position (at the end of phase 1, handle must be at amplitude)
    options.ub(model.ix2(N)) = amplitude;    % upper bound of handle position (at the end of phase 1, handle must be at amplitude)
    options.lb(model.iv2)    = -3.0;         % lower bound of handle velocity
    options.ub(model.iv2)    = 3.0;          % upper bound of handle velocity
    options.lb(model.iLce)   = -3;           % lower bound of fiber length
    options.ub(model.iLce)   = 3;            % upper bound of fiber length
    options.lb(model.ia)     = 0.0;          % lower bound of activation
    options.ub(model.ia)     = 1.0;          % upper bound of activation
    options.lb(model.iu)     = 0.0;          % lower bound of control
    options.ub(model.iu)     = 1.0;          % upper bound of control
    options.lb(model.iT1)    = 0.1;          % lower bound of time
    options.ub(model.iT1)    = duration;     % upper bound of time
    options.lb(model.iT2)    = 0.1;          % lower bound of time
    options.ub(model.iT2)    = duration;     % upper bound of time
    
    % solve the NLP with IPOPT
    funcs.objective         = @objfun;
    funcs.gradient          = @objgrad;
    funcs.constraints       = @confun;
    funcs.jacobian          = @conjac;
    funcs.jacobianstructure = @conjacstructure;
    options.cl              = zeros(model.Ncon,1);
    options.cu              = zeros(model.Ncon,1);	
    options.ipopt.max_iter  = 1000;
    options.ipopt.hessian_approximation = 'limited-memory';
    
    [X, info] = ipopt(X0,funcs,options);
    
    comp.t = toc;
    
    % plot the results
    show(X);    
    
    % calculating time step for 2 phases
    h1 = X(model.iT1)/(N-1);
    h2 = X(model.iT2)/N;
    
    % calculating time nodes in each phase
    times1 = h1*(0:N-1);
    times2 = X(model.iT1)+h2*(1:N);
	
    % store results
    result.x1              = X(1:2*model.N);
    result.v1              = X(2*model.N+(1:2*model.N));
    result.x2              = X(2*2*model.N+(1:2*model.N));
    result.v2              = X(2*3*model.N+(1:2*model.N));
    result.Lce             = X(2*4*model.N+(1:2*model.N));
    result.a               = X(2*5*model.N+(1:2*model.N));
    result.u               = X(2*6*model.N+(1:2*model.N));
    result.time            = [times1'; times2'];
    result.t1              = X(2*7*model.N+1);
    result.t2              = X(2*7*model.N+2);
    result.costfunction    = objfun(X);
    result.MuscleForce     = model.MuscleForce;
    result.CableForce      = model.CableForce;
    result.status          = ipoptmessage(info.status);
    result.info            = info;
    result.ComputationTime = comp.t;
    result.modelparam      = model;
    % save result
    filename = 'result.mat';
    save(filename,'result')
    
end

%==========================================================================
function F = objfun(X)
% objective function: mean of squared controls
    global model
    F = mean(X(model.iu).^2);
end

%==========================================================================
function G = objgrad(X)
% gradient of the objective function
    global model
    G  = zeros(size(X));
    G(model.iu) = X(model.iu)/model.N;
end

%==========================================================================
function c = confun(X)
% constraint function (dynamics constraints)
	global model
	N = model.N;
	% size of constraint vector
	c = zeros(model.Ncon,1);
	% dynamics constraints
	ic = 1:6;		% index for constraint violations c (which depends on # of the equations)
                    % since there are 6 equations, there are 6 possible constraints violations
                    
    Fmuscle = zeros(1,2*N);
    Fcable  = zeros(1,2*N);
    [~,~,~,Fmuscle(1),Fcable(1)] = dynfun(X(1:(2*N):(14*N)),zeros(1,7),1);
    for i=1:2*N-1
        if i < N
            phase = 1;
            h = X(end-1)/(N-1);
        else
            phase = 2;
            h = X(end)/N;
        end
        iy  = i:(2*N):(14*N);
        iyn = iy+1;
        y   = X(iy);
        yn  = X(iyn);
        

        % use Backward Euler formula or Midpoint Euler formula as dynamics constraint
        if strcmp(model.discretization, 'BE')
            [c(ic),~,~,Fmuscle(i+1), Fcable(i+1)] = dynfun(yn , (yn-y)/h, phase);
        else
            [c(ic),~,~,Fmuscle(i+1),Fcable(i+1)] = dynfun((y+yn)/2 , (yn-y)/h, phase);
        end
        ic = ic + 6;
    end

	% periodicity constraints for all states except flywheel position
    c(end-6) = X(model.iu(end)) - X(model.iu(1));           % final and intial control must be equal
	c(end-5) = X(model.iv1(end)) - X(model.iv1(1));         % final and initial velocity of the flywheel should be the same
	c(end-4) = X(model.ix2(end)) - X(model.ix2(1));         % final and initial position of the handle should be the same            
	c(end-3) = X(model.iv2(end)) - X(model.iv2(1));         % final and initial velocity of the handle should be the same
	c(end-2) = X(model.iLce(end)) - X(model.iLce(1));       % final and initial value of the fiber length should be the same  
	c(end-1) = X(model.ia(end)) - X(model.ia(1));           % final and initial value of the muscle activation should be the same  
    
    % constraint on total duration
    c(end)   = X(end-1) + X(end) - model.duration;          % total duration of motion is given
    
    % store muscle force and handle force history
    model.MuscleForce   = Fmuscle';
    model.CableForce    = Fcable';
    
end
%==========================================================================
function J = conjac(X)
    % constraint Jacobian
	global model
    N = model.N;
	% Jacobian matrix dc/dX of the constraints coded in confun
	
	% if no input given, generate a random input to determine Jacobian structure
    if nargin == 0
        X = rand(model.NX,1);
    end
	    
	% size of Jacobian
	J = spalloc(model.Ncon, numel(X), 1000);
	
	% Jacobian of dynamics constraints
	ic = 1:6;		% index for constraint violations c (which depends on # of the equations)
                    % since there are 6 equations, there are 6 possible constraints violations
for i = 1:2*N-1
         if i < N
            phase = 1;
            iT = model.NX - 1;        % T1 index
            h = X(iT) / (N-1);
            dh_dT = 1/(N-1);
         else
            phase = 2;
            iT = model.NX;            % T2 index
            h = X(iT) / (N);
            dh_dT = 1/(N);
        end   
        iy  = i:(2*N):(14*N);
        iyn = iy+1;
        y   = X(iy);
        yn  = X(iyn);
		
        % use Backward Euler formula or Midpoint Euler formula as dynamics constraint
        if strcmp(model.discretization, 'BE')
    		% c(ic) = dynfun(yn , (yn-y)/h, phase);
            [~, dfdy, dfdydot] = dynfun(yn , (yn-y)/h, phase);
            J(ic, iy)  = -dfdydot/(h);            % h1 = T1/N1
            J(ic, iyn) = dfdy + dfdydot/(h);
            
        else
            % c(ic) = dynfun((y+yn)/2 , (yn-y)/h, phase);
            [~, dfdy, dfdydot] = dynfun((y+yn)/2, (yn-y)/h, phase);
            J(ic, iy)  = dfdy/2 - dfdydot/h;
            J(ic, iyn) = dfdy/2 + dfdydot/h;
        end
        % df/dt = (df/dh)*(dh/dT) ==> df/dh = (df/dydot)*(dydot/dh)
        % ydot = (yn-y)/h    dydot/dh = (-1/h^2)*(yn-y)
        J(ic, iT) = dh_dT * dfdydot * (-1/h^2) * (yn-y);
        ic = ic + 6;
end	

% periodicity constraints for all states except flywheel position
% c(end-6) = X(model.iu(end)) - X(model.iu(1));
    J(end-6, model.iu(end)) = 1;
    J(end-6, model.iu(1))   = -1;
% c(end-5) = X(iv1(end)) - X(iv1(1));
    J(end-5, model.iv1(end)) = 1;
    J(end-5, model.iv1(1))   = -1;
% c(end-4) = X(ix2(end)) - X(ix2(1));         
    J(end-4, model.ix2(end)) = 1;
    J(end-4, model.ix2(1))   = -1;
% c(end-3) = X(iv2(end)) - X(iv2(1));         
    J(end-3, model.iv2(end)) = 1;
    J(end-3, model.iv2(1))   = -1;
% c(end-2) = X(iLce(end)) - X(iLce(1));         
    J(end-2, model.iLce(end)) = 1;
    J(end-2, model.iLce(1))   = -1;
% c(end-1) = X(ia(end)) - X(ia(1));       
    J(end-1, model.ia(end)) = 1;
    J(end-1, model.ia(1))   = -1;
    
% constraint on total duration
% c(end) = X(end-1) + X(end) - model.duration;
    J(end, end-1) = 1;
    J(end, end)   = 1;
	
    if (nargin == 0)
        model.Jnnz = nnz(J);
        perc = 100*model.Jnnz / (model.Ncon * numel(X));
        fprintf('Jacobian has %d non-zero elements (%.3f%%)\n', model.Jnnz, perc);
    end
end

%==========================================================================
    function [f, dfdy, dfdydot, Fsee, Fpos] = dynfun(y,ydot,phase)
	% rowing machine dynamics using model parameters p
	% f(y,ydot, Fsee) = 0       [note: Fsee is a function of y!]
	% y = (x1, v1, x2, v2, Lce, a, u)
	% Lce is dimensionless, it is the muscle fiber length divided by Lceopt
	global model

	% copy parameter values into variables with meaningful names
	b1 = model.b1;
	c1 = model.c1;
	m1 = model.m1;
	k1 = model.k1;
    k2 = model.k2;
	b2 = model.b2;
	m2 = model.m2;
    m3 = model.m3;
    
    d   = model.d;
    L   = model.L;
    Lm0 = model.Lm0;
    
    SEELslack = model.SEELslack;
    PEELslack = model.PEELslack;
    Fmax      = model.Fmax;
    Lceopt    = model.Lceopt;
    W         = model.W;
    AHill     = model.AHill;
    Vmax      = model.Vmax;
    FVmax     = model.FVmax;
    
    kPEE2 = 1/W^2;                       % PEE quadratic stiffness, so Fpee = Fmax when Lce = Lce*(1+W)
    kSEE2 = 1/(SEELslack*model.umax)^2;  % SEE quadratic stiffness, so Fsee = Fmax at strain of umax

    % we need to determine the number of nonzeros in Jacobians dfdy and dfdydot
    % number of nonzeros differs in phase 1 and phase 2
	% this also depends on the muscle model used
    % so there are 4 different cases in which nnz differs
       
    % muscle dynamics
    f = zeros(6,1);
    if strcmp(model.muscle_dyn, 'off')          % when the muscle dynamics is off, linear equations are used to get Lcedot, actdot, and Fsee
       f(5) = ydot(5) + y(5);                   % Lcedot = -Lce
       f(6) = ydot(6) + y(6);                   % actdot = -act
       Fsee = y(7)*model.Fmax;                  % Fsee = u*Fmax
        
       dFsee_dy  = [0;0;0;0;0;0;model.Fmax];    % Fsee is a function of y(7)
    if (nargout > 1)
            % before start calculating dfdy and dfdydot we need to
            % determine the number of non-zero elements in Jacobian matrix
            if phase == 1;
            nnz_dfdy    = 13;
            nnz_dfdydot = 6;
            elseif phase == 2;
            nnz_dfdy    = 7;
            nnz_dfdydot = 6;
            end
            % pre-allocate the Jacobian matrix
            dfdy    = spalloc(6,7,nnz_dfdy);            
            dfdydot = spalloc(6,7,nnz_dfdydot);
           % derivative of equations wrt y(5), y(6)
    dfdy(5,5) = 1;
    dfdy(6,6) = 1;
    dfdydot(5,5) = 1;
    dfdydot(6,6) = 1;
    end
    elseif strcmp(model.muscle_dyn, 'on')       % when the muscle dynamics is on, a nonlinear Hill-based muscle model is added to the model to calculate Lcedot, actdot, Fsee
		Lm  = Lm0 - y(3)*(d/L); 	% muscle+tendon length
		Lce = y(5);
		a   = y(6);
		Lcedot = ydot(5);
        
        if (nargout > 1)
        % before start calculating dfdy and dfdydot we need to
        % determine the number of non-zero elements in Jacobian matrices (dfdy & dfdydot)
        if phase == 1;
            nnz_dfdy    = 16;
            nnz_dfdydot = 6;
        elseif phase == 2;
            nnz_dfdy    = 9;
            nnz_dfdydot = 6;
        end
        % pre-allocate the Jacobian matrix
        dfdy    = spalloc(6,7,nnz_dfdy);            
        dfdydot = spalloc(6,7,nnz_dfdydot);
        end
            
		% F1 is the normalized isometric force-length relationship at
		% maximum activation
		ff       = (Lce - 1.0)/W;   % [dimensionless]
		F1       = exp(-ff^2);		% Gaussian force-length curve
		dF1_dLce = -2.0*ff*F1 / W;
		
		% F2 is the dimensionless force-velocity relationship
		if (Lcedot < 0)
			% concentric contraction
			ff          = Vmax - Lcedot/AHill;
			F2          = (Vmax + Lcedot)/ff;
			dF2_dLcedot = (1.0 + F2/AHill)/ff;
		else
			% eccentric contraction
			c           = Vmax * AHill * (FVmax - 1.0) / (AHill + 1.0); % parameter in the eccentric force-velocity equation
			ff          = Lcedot + c;
			F2          = (FVmax*Lcedot + c) / ff;
			dF2_dLcedot = (FVmax - F2)/ff;
		end
	
		% F3 is the dimensionless PEE force (in units of Fmax)
		kPEE     = 1.0/Fmax*Lceopt;             % stiffness of the linear term is 1 N/m, convert to Fmax/Lceopt units	
		ff       = (Lce - PEELslack);           % elongation of PEE, relative to Lceopt
		F3       = kPEE*ff;                     % low stiffness linear term
		dF3_dLce = kPEE;
		if (ff>0)                               % add quadratic term for positive elongation						
			F3 = F3 + kPEE2*ff^2; 
			dF3_dLce = dF3_dLce + 2*kPEE2*ff;
		end
	
		% F4 is the dimensionless SEE force (in units of Fmax)
		kSEE     = 1.0/Fmax;                     % stiffness of the linear term is 1 N/m, convert to Fmax/m
		ff       = Lm - Lce*Lceopt - SEELslack;  % elongation of SEE, in meters
		F4       = kSEE*ff;                      %  low stiffness linear term
		dF4_dLce = -kSEE*Lceopt;
		dF4_dLm  = kSEE;
		if (ff>0)                                % add quadratic term for positive deformation
			F4 = F4 + kSEE2*ff^2;
			dF4_dLce = dF4_dLce - 2 * kSEE2 * Lceopt * ff;
			dF4_dLm  = dF4_dLm + 2 * kSEE2 * ff;
		end

		% F5 is viscous damping parallel to the CE (0.001 of Fmax at 1 Lceopt/s) to
		% ensure that df/dLcedot is never zero
		F5 = 0.001*Lcedot;
		dF5_dLcedot = 0.001;
		
		% The muscle dynamics equation: f = Fsee - Fce - Fpee - Fdamping = 0
		f(5) = F4 - a*F1*F2 - F3 - F5;
		if (nargout > 1)
			dfdy(5,6)    = -F1*F2;                              % df/da
			dfdy(5,5)    = dF4_dLce - a*dF1_dLce*F2 - dF3_dLce;	% df/dLce
			dfdydot(5,5) = -a*F1*dF2_dLcedot - dF5_dLcedot;		% df/dLcedot
			df_dLm       = dF4_dLm;
			dfdy(5,3)    = -(d/L)* df_dLm;                      % derivative of f with respect to hand position
		end
		
		% Force in SEE is Fmax*F4
		Fsee = Fmax*F4;
		if (nargout > 1)
			dFsee_dLm  = Fmax * dF4_dLm;
			dFsee_dLce = Fmax * dF4_dLce;
			dFsee_dy   = [0;0;(-d/L)*dFsee_dLm; 0;dFsee_dLce;0;0];  % Fsee is a function of y(3) & y(5)
		end
		
		% Activation dynamics equation
		f(6) = ydot(6)-(y(7)-y(6))*(y(7)/model.Tact + (1-y(7))/model.Tdeact);
		if (nargout > 1)
			dfdy(6,6) = (y(7)/model.Tact + (1-y(7))/model.Tdeact);
			dfdy(6,7) = -(y(7)/model.Tact + (1-y(7))/model.Tdeact) ...
						-(y(7)-y(6))*(1/model.Tact - 1/model.Tdeact);
			dfdydot(6,6) = 1;
		end
	else
		error('muscle_dyn option must be either on or off');
    end

	% rowing machine dynamics f(1,2,3,4)=0
    if phase == 1;
        epsilon = 0.5;
        F    =  k1*(y(3)-y(1)) + b2*(y(4)-y(2));
        Fpos = 0.5*(F+sqrt(F^2 + epsilon^2)) - 0.5*epsilon;
        
		f(1) = ydot(1) - y(2);
		f(2) = m1*ydot(2) + b1*y(2) + c1*y(2)^2 - Fpos;
		f(3) = ydot(3) - y(4);
		f(4) = (m2 + m3)*ydot(4) - d/L*Fsee + Fpos + k2*y(3);
        
        if (nargout > 1)
            D         = 0.5*(1 + F/sqrt(F^2 + epsilon^2));
            % derivative of equations wrt y(1), y(2), y(3), y(4)
            dfdy(1,2) = -1;
            dfdy(2,1) = k1*D;
            dfdy(2,2) = b2*D + b1 + 2*c1*y(2);
            dfdy(2,3) = -k1*D;
            dfdy(2,4) = -b2*D;
            dfdy(3,4) = -1;
            dfdy(4,1) = -k1*D;
            dfdy(4,2) = -b2*D; 
            dfdy(4,3) = k1*D + k2;
            dfdy(4,4) = b2*D;
            dfdy(4,:) = dfdy(4,:) - d/L*dFsee_dy';
            
            % derivative of equations wrt dy(1), dy(2), dy(3), dy(4)
            dfdydot(1,1) = 1;
            dfdydot(2,2) = m1;
            dfdydot(3,3) = 1;
            dfdydot(4,4) = m2 + m3;
        end
    elseif phase == 2;
        Fpos = 0;
        f(1) = ydot(1) - y(2);
        f(2) = m1*ydot(2) + b1*y(2) + c1*y(2).^2;
        f(3) = ydot(3) - y(4);
        f(4) = (m2 + m3)*ydot(4) - d/L*Fsee + k2*y(3);
    
        if (nargout > 1)
            % derivative of equations wrt y(1), y(2), y(3), y(4)
            dfdy(1,2) = -1;
            dfdy(2,2) = b1 + 2*c1*y(2);
            dfdy(3,4) = -1;
            dfdy(4,3) = k2;
            dfdy(4,:) = dfdy(4,:) - d/L*dFsee_dy';
            
            % derivative of equations wrt dy(1), dy(2), dy(3), dy(4)
            dfdydot(1,1) = 1;
            dfdydot(2,2) = m1;
            dfdydot(3,3) = 1;
            dfdydot(4,4) = (m2 + m3);

        end
    end

end
%==========================================================================

	function J = conjacstructure(X)
    global model

      J = model.Jstructure;
			
    end
    %======================================================================
    function checkderiv(X)
	% use finite differences to check that the code in objgrad and conjac is correct
	global model
    Ncon = model.Ncon;
    NX   = model.NX;
    
    X        = randn(NX,1);
    hh       = 1e-6;
	f        = objfun(X);
	grad     = objgrad(X);
	c        = confun(X);
	cjac     = conjac(X);
	cjac_num = spalloc(Ncon, NX,10000);
	grad_num = zeros(NX,1);
	
    for i = 1:NX
        fprintf('checking derivatives for unknown %4d of %4d\n',i,NX);
        Xisave        = X(i);
        X(i)          = X(i) + hh;
        cjac_num(:,i) = sparse(confun(X) - c)/hh;
        grad_num(i)   = (objfun(X) - f)/hh;
        X(i)          = Xisave;
    end
	
	% report maximal differences between analytical derivatives and numerical results
	fprintf('Max. error in constraint jacobian: ');
	matcompare(cjac, cjac_num);
	fprintf('Max. error in objective gradient: ');
	matcompare(grad, grad_num);
	disp('Hit ENTER to continue');
    pause;
    end

    %======================================================================
    function matcompare(a,b)
	% compares two matrices and prints element that has greatest difference
	[maxerr,irow] = max(abs(a-b));
	[maxerr,icol] = max(maxerr);
	irow = irow(icol);
	fprintf('%9.6f at %d %d (%9.6f vs. %9.6f)\n', full(maxerr), irow, icol, full(a(irow,icol)), full(b(irow,icol)));
    end

    %======================================================================
	function show(X)
    global model
    N = model.N;

    % plot the current solution
    x1  = X(model.ix1);
    v1  = X(model.iv1);
    x2  = X(model.ix2);
    v2  = X(model.iv2);
    Lce = X(model.iLce);
    a   = X(model.ia);
    u   = X(model.iu);

    % calculating time step for 2 phases
    h1 = X(model.iT1)/(N-1);
    h2 = X(model.iT2)/N;
    
    % calculating time nodes in each phase
    times1 = h1*(0:N-1);
    times2 = X(model.iT1)+h2*(1:N);
	time   = [times1'; times2'];
   
    figure(1)
    set(gcf, 'Color', 'w');
  	set(gcf,'Position',[10 10 800 900]);

    subplot(3,3,1);set(gca,'fontsize',14);h1 = plot(time(1:N),x1(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);h2 = plot(time(N:end),x1(N:end),'r','linewidth',2);title('Flywheel Position(m)');
    legend([h1 h2], 'Phase 1', 'Phase 2')
       
    subplot(3,3,2);set(gca,'fontsize',14);plot(time(1:N),v1(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end),v1(N:end),'r','linewidth',2);title('Flywheel Velocity(m/s)');
    subplot(3,3,3);set(gca,'fontsize',14);plot(time(1:N),x2(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end),x2(N:end),'r','linewidth',2);title('Handle Position(m)');
    subplot(3,3,4);set(gca,'fontsize',14);plot(time(1:N),v2(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end),v2(N:end),'r','linewidth',2);title('Handle Velocity(m/s)');
    subplot(3,3,7);set(gca,'fontsize',14);plot(time(1:N),Lce(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end),Lce(N:end),'r','linewidth',2);title('Lce/Lceopt');xlabel('Time(s)');
    subplot(3,3,6);set(gca,'fontsize',14);plot(time(1:N),a(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end),a(N:end),'r','linewidth',2);title('Muscle Activation');
    subplot(3,3,5);set(gca,'fontsize',14);plot(time(1:N),u(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end),u(N:end),'r','linewidth',2);title('Control');
    subplot(3,3,8);set(gca,'fontsize',14);plot(time(1:N), model.MuscleForce(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end), model.MuscleForce(N:end),'r','linewidth',2);title('Muscle Force (N)');xlabel('Time(s)');
    subplot(3,3,9);set(gca,'fontsize',14);plot(time(1:N), model.CableForce(1:N),'linewidth',2);
    hold on;set(gca,'fontsize',14);plot(time(N:end), model.CableForce(N:end),'r','linewidth',2);title('Cable Force (N)');xlabel('Time(s)');

    end
    
%==============================================================
function [s] = ipoptmessage(info)

 	if info == 0; s = 'solved';return;end;
 	if info == 1; s = 'solved to acceptable level'; return; end;
 	if info == 2; s = 'infeasible problem detected'; return; end;
 	if info == 3; s = 'search direction becomes too small'; return; end;
 	if info == 4; s = 'diverging iterates'; return; end;
 	if info == 5; s = 'user requested stop'; return; end;
     
 	if info == -1; s = 'maximum number of iterations exceeded'; return; end;
 	if info == -2; s = 'restoration phase failed'; return; end;
 	if info == -3; s = 'error in step computation'; return; end;
 	if info == -10; s = 'not enough degrees of freedom'; return; end;
 	if info == -11; s = 'invalid problem definition'; return; end;
 	if info == -12; s = 'invalid option'; return; end;
 	if info == -13; s = 'invalid number detected'; return; end;

 	if info == -100; s = 'unrecoverable exception'; return; end;
 	if info == -101; s = 'non-IPOPT exception thrown'; return; end;
 	if info == -102; s = 'insufficient memory'; return; end;
 	if info == -199; s = 'internal error'; return; end;
  
end