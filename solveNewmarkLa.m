function varargout=solveNewmarkLa(beta,gamma,ts,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap,varargin)
% Newmark's method for a linear problem.
%
% INPUT:
% beta - Newmark beta parameter.
% gamma - Newmark gamma parameter.
% ts - #TimeSteps+1-by-1 matrix of integration time points.
% fDoF - list of IDs of free degrees of freedom.
% pDoF - list of IDs of prescribed (fixed) degrees of freedom.
% M - (#Nodes*#Eq)-by-(#Nodes*#Eq) mass matrix or a function handle
%   returning such matrix. The input of the function is t_ (time). 
% C - (#Nodes*#Eq)-by-(#Nodes*#Eq) damping matrix or a function handle
%   returning such matrix. The input of the function is t_ (time). 
% K - (#Nodes*#Eq)-by-(#Nodes*#Eq) stiffness matrix or a function handle
%   returning such matrix. The input of the function is t_ (time). 
% F - (#Nodes*#Eq)-by-1 external force vector or a function handle
%   returning such vector. The input of the function is t_ (time). 
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
% 'increment', outputFunction - function called at each time increment and
%   supplied with (i,ts(i),u,v,a,fDoF,pDoF,Mval,Cval,Kval,Fval,u0,v0,a0,upval,vpval,apval) as
%   input. The output must be one or several matrices. The outputs are
%   collected and returned at the end of the procedure. 
%
% OUTPUT:
% If no 'increment' option specified, then:
%   u - (#Nodes*#Eq)-by-1 vector of displacements at the final time point.    
%   v - (#Nodes*#Eq)-by-1 vector of velocities at the final time point.    
%   a - (#Nodes*#Eq)-by-1 vector of accelerations at the final time point.
% If 'increment' option specified, then:
%   argsOut - a cell array containing the outputs of outputFunction
%     collected at each time step or Newton iteration. Each element
%     argsOut{k} contains NI-by-n1(k)-by-...-nq(k) matrix, where
%     n1(k)-by-...-nq(k) is the size of k-th output of outputFunction and
%     NI is the total number of increments.
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
%
% Examples:
% u=solveNewmarkLa(beta,gamma,ts,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap);
% outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
% argsOut=solveNewmarkLa(beta,gamma,ts,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap,...
%   'increment',outputFunction);

%% Parsing input
numGDoF=size(u0,1);

% Make sure all M,C,K,F,up,vp,ap are function handles of type @(t).
if ~isa(M,'function_handle')
    M=@(t)M;
elseif nargin(M)~=1
    error('Too many input parameters in the function handle provided');
end
if isempty(C)
    if ~isa(K,'function_handle')
        C=zeros(size(K));
    else
        C=zeros(size(K(ts(1))));
    end
end
if ~isa(C,'function_handle')
    C=@(t)C;
elseif nargin(C)~=1
    error('Too many input parameters in the function handle provided');
end
if ~isa(K,'function_handle')
    K=@(t)K;
elseif nargin(K)~=1
    error('Too many input parameters in the function handle provided');
end
if ~isa(F,'function_handle')
    F=@(t)F;
elseif nargin(F)~=1
    error('Too many input parameters in the function handle provided');
end
if isempty(up)
    up=zeros(numGDoF,1);
end
if ~isa(up,'function_handle')
    up=@(t)up;
end
if isempty(vp)
    vp=zeros(numGDoF,1);
end
if ~isa(vp,'function_handle')
    vp=@(t)vp;
end
if isempty(ap)
    ap=zeros(numGDoF,1);
end
if ~isa(ap,'function_handle')
    ap=@(t)ap;
end

incOutput=0;
if isempty(varargin)
elseif mod(numel(varargin),2)==0
    for i=1:(numel(varargin)/2)
        switch varargin{(i-1)*2+1}
            case 'increment'
                incOutput=1;
                incOutputFun=varargin{i*2};
            case 'iteration'
                %ignored
            case 'maxiter'
                %ignored
        end
    end
end

% Boundary conditions
gDoF=[1:numGDoF]';
if iscell(pDoF) % If Dirichlet BC are mixed (e.g. displacemenet, velocity etc combined)
    pDoFu=pDoF{1}; pDoFv=pDoF{2}; pDoFa=pDoF{3};
    fDoFu=setdiff(gDoF,pDoFu);
    fDoFv=setdiff(gDoF,pDoFv);
    fDoFa=setdiff(gDoF,pDoFa);
    pDoF=union([pDoFu;pDoFv;pDoFa],[]);
else % Otherwise we just maintain the generality
    pDoFu=pDoF; pDoFv=[]; pDoFa=[];
    fDoFu=fDoF; fDoFv=gDoF; fDoFa=gDoF;
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
i=1;
if incOutput
    % evaluate initial values of M, C, K, F
    M_=M(ts(1));C_=C(ts(1));K_=K(ts(1));F_=F(ts(1));
    out_temp=incOutputFun(i,ts(i),u,v,a,fDoF,pDoF,M_,C_,K_,F_,u0,v0,a0,up,vp,ap);
    nOutput=numel(out_temp);
    output=cell(nOutput,1);
    for k=1:nOutput
        output{k}=zeros([numTimesteps size(out_temp{k})]);
        output{k}(i,:)=out_temp{k}(:);
    end
end

% time incrementation
for i=2:numTimesteps
    dt=ts(i)-ts(i-1);
    u_tilde=u + dt*v + 1/2*dt^2*(1-2*beta)*a; %u predictor
    v_tilde=v + dt*(1-gamma)*a; %v predictor
    
    % evaluate M, C, K, F
    M_=M(ts(i));
    C_=C(ts(i));
    K_=K(ts(i));
    F_=F(ts(i));
    
    % define a at prescribed DoF (boundary conditions)
    ap_=ap(ts(i));
    vp_=vp(ts(i));
    up_=up(ts(i));

    a(pDoFa)=ap_(pDoFa);
    v(pDoFv)=vp_(pDoFv);
    u(pDoFu)=vp_(pDoFu);
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
    
    % compute new accelerations
    rhs=F_-C_*v_tilde-K_*u_tilde; % right-hand side
    A=M_ + dt*gamma*C_ + dt^2*beta*K_; % operator on the left-hand side
    a(fDoF)=A(fDoF,fDoF) \ (rhs(fDoF) - A(fDoF,pDoF)*a(pDoF));

    % update displacement and velocities
    u(fDoF)=u_tilde(fDoF)+1/2*dt^2*2*beta*a(fDoF);
    v(fDoF)=v_tilde(fDoF)+dt*gamma*a(fDoF);

    if incOutput
        out_temp=incOutputFun(i,ts(i),u,v,a,fDoF,pDoF,M_,C_,K_,F_,u0,v0,a0,up,vp,ap);
        for k=1:nOutput
            output{k}(i,:)=out_temp{k}(:);
        end
    end
    
end


if ~(incOutput)
    varargout={u,v,a};
else
    varargout={output};
end

figure(fig);
end %function