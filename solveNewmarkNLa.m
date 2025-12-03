function varargout=solveNewmarkNLa(beta,gamma,ts,fDoF,pDoF,M,Ct,Kt,Ft,u0,v0,a0,up,vp,ap,varargin)
% Newmark's method for a nonlinear problem.
%
% INPUT:
% beta - Newmark beta parameter.
% gamma - Newmark gamma parameter.
% ts - #TimeSteps+1-by-1 matrix of integration time points.
% fDoF - list of IDs of free degrees of freedom.
% pDoF - list of IDs of prescribed (fixed) degrees of freedom.
% M - (#Nodes*#Eq)-by-(#Nodes*#Eq) mass matrix or a function handle
%   returning such matrix. The input of the function is t_, u_, v_ -
%   time and #Nodes*#Eq-by-1 displacement and velocity vectors respectively.   
% Ct - (#Nodes*#Eq)-by-(#Nodes*#Eq) tangent damping matrix or a function
%   handle returning such matrix. The input of the function is t_, u_, v_ -
%   time and #Nodes*#Eq-by-1 displacement and velocity vectors 
%   respectively.  
% Kt - (#Nodes*#Eq)-by-(#Nodes*#Eq) tangent stiffness matrix or a function
%   handle returning such matrix. The input of the function is t_, u_, v_ -
%   time and #Nodes*#Eq-by-1 displacement and velocity vectors 
%   respectively.  
% Ft - (#Nodes*#Eq)-by-1 total force vector (Fint-Fext) or a function
%   handle returning such vector. The input is of the function t_, u_, v_ -
%   time and #Nodes*#Eq-by-1 displacement and velocity vectors 
%   respectively.  
% u0 - (#Nodes*#Eq)-by-1 vector of initial displacements.
% v0 - (#Nodes*#Eq)-by-1 vector of initial velocities or an empty array.
% a0 - (#Nodes*#Eq)-by-1 vector of initial accelerations or an empty array.
% up - (#Nodes*#Eq)-by-1 vector of prescribed displacements or a function
%   handle returning such vector, or an empty array. The input of the
%   function handle is t_ (time).
% vp - (#Nodes*#Eq)-by-1 vector of prescribed velocities or a function
%   handle returning such vector, or an empty array. The input of the
%   function handle is t_ (time).
% ap - (#Nodes*#Eq)-by-1 vector of prescribed accelerations or a function
%   handle returning such vector, or an empty array. The input of the
%   function handle is t_ (time).
% Options in the format ...,'Name',value,...
% 'iteration', outputFunction - function called at each iteration (of
%   Newton's iterative method) and supplied with ([i j],ts(i),u,v,a,fDoF,
%   pDoF,Mval,Cval,Kval,Fval,u0,v0,a0,upval,vpval,apval) as input. The
%   output must be one or several matrices. The outputs are collected and
%   returned at the end of the procedure.
% 'increment', outputFunction - function called at each time increment and
%   supplied with (i,ts(i),u,v,a,fDoF,pDoF,Mval,Cval,Kval,Fval,u0,v0,a0,upval,vpval,apval)
%   as input. The output must be one or several matrices. The outputs are
%   collected and returned at the end of the procedure. 
% 'maxiter', maxIter - maximum number of Newton's iterations at each
%   timestep. If this number is reached, Newton method aborts and an error
%   message is thrown. Default value - 10. 
% 'reltol', relTol - relative tolerance. Iterations stop if the norm
%   of the residual is within this tolerance relatively to the total force
%   vector.
% 'abstol', absTol - absolute tolerance. Iterations stop if the norm of the
%   residual is less or equal to absTol.
%
% OUTPUT:
% If no 'iteration' or 'increment' option specified, then:
%   u - (#Nodes*#Eq)-by-1 vector of displacements at the final time point.    
%   v - (#Nodes*#Eq)-by-1 vector of velocities at the final time point.    
%   a - (#Nodes*#Eq)-by-1 vector of accelerations at the final time point.
% If 'iteration' or 'increment' option specified, then:
%   argsOut - a cell array containing the outputs of outputFunction
%     collected at each time step or Newton iteration. Each element
%     argsOut{k} contains NI-by-n1(k)-by-...-nq(k) matrix, where
%     n1(k)-by-...-nq(k) is the size of k-th output of outputFunction and
%     NI is the total number of iterations or increments.
%
% COMMENTS:
% (I) By default, the Dirichlet (essential) boundary conditions are given
% by prescibed displacements. In order to allow for different boundary
% conditions definition, pDoF can be replaced by {pDoFu, pDoFv, pDoFa} -
% lists of nodes at which displacements, velocities and accelerations are
% prescribed respectively. The corresponding entries of up, vp, and ap are
% used directly and via finite difference approximations, which normally
% match the approximation used in Newmark method.
% (II) The mass matrix must be non-singular (e.g. a lumped diagonal matrix)
% (III) Empty array [] can be passed as a value for input parameters C, up,
% vp, ap, v0 and a0. In this case the arguments are replaced by zero
% arrays of appropriate size.
% (IV) If both 'iteration' and 'increment' options are specified, then
% outputs of the former are not collected, yet the corresponding
% outputFunction is called at the end of each iteration and can be used
% e.g. for plotting.
% (V) Output collection at each iteration is inefficient as the target
% array changes its size every time.
%
% Examples:
% u=solveNewmarkNLa(beta,gamma,ts,fDoF,pDoF,M,Ct,Kt,Ft,u0,v0,a0,up,vp,ap);
% u=solveNewmarkNLa(beta,gamma,ts,fDoF,pDoF,M,C,K,@(t,u,v)(K*u-Fext),...
%   u0,v0,a0,up,vp,ap); %linear problem
% outputFunction1=@(i,t,u,v,a,fDoF,pDoF,M,Ct,Kt,Ft,u0,v0,a0,up,vp,ap){t,u,v,a};
% outputFunction2=@(i,t,u,v,a,fDoF,pDoF,M,Ct,Kt,Ft,u0,v0,a0,up,vp,ap)...
%   plot(t,u(1),'*');
% argsOut=solveNewmarkLa(beta,gamma,ts,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap,...
%   'increment',outputFunction2,'iteration',outputFunction1);

maxiter=10; %number of maximum Newton iteration at each
reltol=1e-3; %relative tolerance for the residual
abstol=1e-8; %absolute tolerance for the residual
%% Parsing input
numGDoF=size(u0,1);

% Make sure all M,C,K,F are function handles of type @(t,u,v)
if ~isa(M,'function_handle')
    M=@(t,u,v)M;
elseif nargin(M)~=3
    error('Too many input parameters in the function handle provided');
end
if isempty(Ct)
    Ct=zeros(numGDoF,numGDoF);
end
if ~isa(Ct,'function_handle')
    Ct=@(t,u,v)Ct;
elseif nargin(Ct)~=3
    error('Too many input parameters in the function handle provided');
end
if ~isa(Kt,'function_handle')
    Kt=@(t,u,v)Kt;
elseif nargin(Kt)~=3
    error('Too many input parameters in the function handle provided');
end
if ~isa(Ft,'function_handle')
    Ft=@(t,u,v)Ft;
elseif nargin(Ft)~=3
    error('Too many input parameters in the function handle provided');
end
% Make sure up, vp, ap are function handles of type @(t)
if ~isa(up,'function_handle')
    up=@(t)up;
end
if ~isa(vp,'function_handle')
    vp=@(t)vp;
end
if ~isa(ap,'function_handle')
    ap=@(t)ap;
end

incOutput=0;
iterOutput=0;
if isempty(varargin)
    %nop
elseif mod(numel(varargin),2)==0
    for i=1:(numel(varargin)/2)
        switch varargin{(i-1)*2+1}
            case 'increment'
                incOutput=1;
                incOutputFun=varargin{i*2};
            case 'iteration'
                iterOutput=1;
                iterOutputFun=varargin{i*2};
            case 'maxiter'
                maxiter=varargin{i*2};
            case 'reltol'
                reltol=varargin{i*2};
            case 'abstol'
                abstol=varargin{i*2};
        end
    end
end

% Boundary conditions
gDoF=[1:numGDoF]';
if iscell(pDoF) % If Dirichlet BC are mixed (e.g. displacemenet, velocity etc combined)
    pDoFu=pDoF{1}; pDoFv=pDoF{2}; pDoFa=pDoF{3};
    fDoFu=setdiff(gDoF,pDoFu);    fDoFv=setdiff(gDoF,pDoFv);    fDoFa=setdiff(gDoF,pDoFa); %#ok<NASGU>
    pDoF=union([pDoFu;pDoFv;pDoFa],[]);
else % Otherwise we just maintain the generality
    pDoFu=pDoF; pDoFv=[]; pDoFa=[];
    fDoFu=fDoF; fDoFv=gDoF; fDoFa=gDoF; %#ok<NASGU>
end

%% Initial conditions
% Make initial conditions comply with the boundary conditions.
if isempty(a0)
    a0=zeros(size(u0));
end
if isempty(v0)
    a0=zeros(size(u0));
end
ap_=ap(ts(1));
vp_=vp(ts(1));
up_=up(ts(1));
u0(pDoFu)=up_(pDoFu);
v0(pDoFv)=vp_(pDoFv);
a0(pDoFa)=ap_(pDoFa);

%%
fig=gcf; % outputFunction may change current figure. We revert to current figure in the end
numTimesteps=numel(ts);

% Initialise u, v, a
a=a0; v=v0; u=u0;

% output at the initial time
i=1; % time incrementation counter
jj=0; % total counter of Newton iterations
if incOutput&&~iterOutput
    % evaluate initial values of M, C, K, F
    M_=M(ts(1),u0,v0);Ct_=Ct(ts(1),u0,v0);Kt_=Kt(ts(1),u0,v0);Ft_=Ft(ts(1),u0,v0);
    out_temp=incOutputFun(i,ts(i),u,v,a,fDoF,pDoF,M_,Ct_,Kt_,Ft_,u0,v0,a0,up,vp,ap);
    nOutput=numel(out_temp);
    output=cell(nOutput,1);
    for k=1:nOutput
        output{k}=zeros([numTimesteps size(out_temp{k})]);
        output{k}(i,:)=out_temp{k}(:);
    end
end
if iterOutput
    % evaluate initial values of M, C, K, F
    M_=M(ts(1),u0,v0);Ct_=Ct(ts(1),u0,v0);Kt_=Kt(ts(1),u0,v0);Ft_=Ft(ts(1),u0,v0);
    out_temp=iterOutputFun([i,0],ts(i),u,v,a,...
        fDoF,pDoF,M_,Ct_,Kt_,Ft_,u0,v0,a0,up_,vp_,ap_);
    nOutput=numel(out_temp);
    output=cell(nOutput,1);
    for k=1:nOutput
        output{k}=zeros([numTimesteps size(out_temp{k})]);
        output{k}(jj+1,:)=out_temp{k}(:); % changes size every iteration
    end
end


% time incrementation
for i=2:numTimesteps
    dt=ts(i)-ts(i-1);
    ap_=ap(ts(i));    vp_=vp(ts(i));    up_=up(ts(i));

    u_tilde=u + dt*v + 1/2*dt^2*(1-2*beta)*a; %u predictor
    v_tilde=v + dt*(1-gamma)*a; %v predictor
    v_tilde(pDoFv)=vp_(pDoFv);
    u_tilde(pDoFu)=up_(pDoFu);

    u=u_tilde; v=v_tilde; % initial guess corresponds to predictor
    a=zeros(size(u)); % adjusted to essential boundary conditions
    % define a at prescribed DoF (essential boundary conditions)
    if gamma ~=0
        a(pDoFv)=1/dt/gamma*(vp_(pDoFv)-v_tilde(pDoFv));
    else
        a(pDoFv)=1/dt/0.5*(vp_(pDoFv)-v_tilde(pDoFv)); % as if gamma = 1/2
    end
    if beta ~=0
        a(pDoFu)=1/dt^2/beta*(up_(pDoFu)-u_tilde(pDoFu));
    else
        a(pDoFu)=1/dt^2/0.25*(up_(pDoFu)-u_tilde(pDoFu));% as if beta = 1/4
    end
    a(pDoFa)=ap_(pDoFa);
   
    j=0;%iteration counter at current time step
    jj=jj+1;
     
    da=zeros(size(a)); %init acceleration correction
    % evaluate the residual at the predited state
    Ct_=Ct(ts(i),u,v);
    Kt_=Kt(ts(i),u,v);
    M_=M(ts(i),u,v);
    Ft_=Ft(ts(i),u,v);
    r=M_*a+Ft_;
    if iterOutput % output at the first (predictor) iteration
        out_temp=iterOutputFun([i,j],ts(i),u,v,a,...
            fDoF,pDoF,M_,Ct_,Kt_,Ft_,u0,v0,a0,up_,vp_,ap_);
        for k=1:nOutput
            output{k}(jj+1,:)=out_temp{k}(:); % changes size every iteration
        end
    end
    
    while (j<maxiter)&&...
        ( (norm(r(fDoF))>abstol)||(norm(r(fDoF))>reltol*norm(Ft_)) )
        j=j+1; jj=jj+1;
        
        % solve for acceletaion correction
        A = M_ + Ct_*dt*gamma + Kt_*dt^2*beta;
        da(fDoF)=-A(fDoF,fDoF)\(r(fDoF)+A(fDoF,pDoF)*a(pDoF));

        % update accelerations, velocities and displacements
        a(fDoF)=a(fDoF)+da(fDoF);
        % update displacement and velocities
        u(fDoF)=u_tilde(fDoF)+1/2*dt^2*2*beta*a(fDoF);
        v(fDoF)=v_tilde(fDoF)+dt*gamma*a(fDoF);
        
        % evaluate new residual at the corrected state
        M_=M(ts(i),u,v);
        Ft_=Ft(ts(i),u,v);
        Ct_=Ct(ts(i),u,v);
        Kt_=Kt(ts(i),u,v);
        r=M_*a+Ft_;

        if iterOutput
            out_temp=iterOutputFun([i,j],ts(i),u,v,a,...
                fDoF,pDoF,M_,Ct_,Kt_,Ft_,u0,v0,a0,up_,vp_,ap_);
            for k=1:nOutput
                output{k}(jj+1,:)=out_temp{k}(:); % changes size every iteration
            end
        end
    end
    
    if (norm(r(fDoF))>abstol)||(norm(r(fDoF))>reltol*norm(Ft_))
        error('Newton''s iterative method failed to converge');
    end
    
    if incOutput
        M_=M(ts(1),u0,v0);Ct_=Ct(ts(1),u0,v0);Kt_=Kt(ts(1),u0,v0);Ft_=Ft(ts(1),u0,v0);
        out_temp=incOutputFun(i,ts(i),u,v,a,fDoF,pDoF,M_,Ct_,Kt_,Ft_,u0,v0,a0,up,vp,ap);
        if ~iterOutput
            for k=1:nOutput
                output{k}(i,:)=out_temp{k}(:);
            end
        end
    end
end


if ~(incOutput||iterOutput)
    varargout={u,v,a};
else
    varargout={output};
end

figure(fig); % revert to current figure, as it was in the beginning.
end %function