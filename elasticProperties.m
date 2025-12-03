function prop=elasticProperties(varargin)
% Computes various representations for properties of a linearly elastic
% material defined in terms of any two parameters.
%
% INPUT:
% parName1, parVal1, parName2, parVal2 - where
%   parName1, parName2 - strings, parameter names (see below).
%   parVal1, parVal2 - real scalars, parameter values.
% fieldName (optional) - a single string specifying a requested property of the
%   structure
% OUTPUT:
% prop - elasticProperties structure, if no fieldName specified. Otherwise,
%   the value of the requested field.
%
% Examples:
%   prop = elasticProperties('lambda',8e7,'mu',8e7);
%   prop = elasticProperties('youngsModulus',2e8,'poissonsRatio',0.25);
%   prop.CFull % 9-by-9 matrix
%   elasticProperties('lambda',8e7,'mu',8e7,'CPlaneStressFull') % 4-by-4
%       % matrix
%
% elasticProperties structure fields:
%  'lambda', ('la' can be used for input) - first Lame parameter,
%  'mu' - shear modulus, second Lame parameter,
%  'youngsModulus', ('E'  can be used for input) - Young's modulus,
%  'poissonsRatio', ('nu'  can be used for input)- Poisson's ratio,
%  'CPlaneStressFull' (output only) - full form plane stress 4-by-4
%    elasticity tensor, 
%  'CPlaneStrainFull' (output only) - full form plane strain 4-by-4
%    elasticity tensor,  
%  'CPlaneStressEng' (output only) - short form plane stress 3-by-3
%    elasticity tensor,  
%  'CPlaneStrainEng' (output only) - short form plane strain 3-by-3
%    elasticity tensor.  
%  'CFull' (output only) - full form 3D 9-by-9 elasticity tensor
%  'CShort' (output only) - short form 3D 6-by-6 elasticity tensor


la=[];
mu=[];
E=[];
nu=[];
if (numel(varargin)~=4)&&(numel(varargin)~=5)
    error('Two name-value pairs expected as input');
end
for i=1:2
    name=varargin{(i-1)*2+1};
    value=varargin{(i-1)*2+2};
    switch name
        case {'lambda','la'}
            la=value;
        case {'mu','G'}
            mu=value;
        case {'youngsModulus','E'}
            E=value;
        case {'poissonsRatio','nu'}
            nu=value;
        otherwise
            error(['Unidentified elastic property: ' name])
    end
end

if ~isempty(la)&&~isempty(E)
    mu=(E-3*la+sqrt(E^2+9*la^2+2*E*la))/4;
elseif ~isempty(la)&&~isempty(nu)
    mu=la*(1-2*nu)/2/nu;
elseif ~isempty(mu)&&~isempty(E)
    la=mu*(E-2*mu)/(3*mu-E);
elseif ~isempty(mu)&&~isempty(nu)
    la=2*mu*nu/(1-2*nu);
elseif ~isempty(E)&&~isempty(nu)
    la=E*nu/(1+nu)/(1-2*nu);
    mu=E/2/(1+nu);
elseif ~isempty(mu)&&~isempty(la)
    %do nothing
else
    error('something went wrong');
end %if

E=mu*(3*la+2*mu)/(la+mu);
nu=la/2/(la+mu);

planeStrainCFull=...
    [la+2*mu,0,0,la;...
    0,mu,mu,0;...
    0,mu,mu,0;...
    la,0,0,la+2*mu];

planeStrainCEng=...
    [la+2*mu,la,0;...
    la,la+2*mu,0;...
    0,0,mu];

planeStressCFull=...
    [la+2*mu-la^2/(2*mu+la),0,0,la-la^2/(2*mu+la);...
    0,mu,mu,0;...
    0,mu,mu,0;...
    la-la^2/(2*mu+la),0,0,la+2*mu-la^2/(2*mu+la)];

planeStressCEng=...
    [la+2*mu-la^2/2/mu,la-la^2/2/mu,0;...
    la-la^2/2/mu,la-la^2/2/mu+2*mu,0;...
    0,0,mu];

CShort=...
    [la+2*mu, la, la, 0, 0, 0;...
    la, la+2*mu, la, 0, 0, 0;...
    la, la, la+2*mu, 0, 0, 0;...
    0,0,0,mu,0,0;
    0,0,0,0,mu,0;
    0,0,0,0,0,mu];

CFull=[...  
    la+2*mu,0,0,    0,la,0,         0,0,la;
    0,mu,0,         mu,0,0,         0,0,0;
    0,0,mu,         0,0,0,          mu,0,0;

    0,mu,0,         mu,0,0,         0,0,0;
    la,0,0,         0,la+2*mu,0,    0,0,la;
    0,0,0,          0,0,mu,         0,mu,0;

    0,0,mu,         0,0,0,          mu,0,0;
    0,0,0,          0,0,mu,         0,mu,0;
    la,0,0          0,la,0,         0,0,la+2*mu ...
    ];

prop=struct(...
    'lambda',la,...
    'mu',mu,...
    'youngsModulus',E,...
    'poissonsRatio',nu,...
    'CPlaneStrainFull',planeStrainCFull,...
    'CPlaneStrainEng',planeStrainCEng,...
    'CPlaneStressFull',planeStressCFull,...
    'CPlaneStressEng',planeStressCEng,...
    'CShort',CShort,....
    'CFull',CFull ...
    );

if numel(varargin)==5
    prop=prop.(varargin{5});
end
end %function