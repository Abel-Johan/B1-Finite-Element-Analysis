function B=tensorRotate(A,varargin)
% Rotate a 2D second-order tensor A by angle alpha
%
% INPUT:
%   A - k-by-4 matrix, k-th row represents a 2d tensor as a11, a21, a12, a22
%       k-by-3 matrix, k-th row represents a 2d strain tensor as e11, e22, g12,
%         if form='engineering', 
%       or as e11, e22, e12, if form='tensorial'
%   alpha - a scalar or k-by-m matrix or 1-by-m matrix of rotation angles.
%   form (when A is k-by-3) - a string 'engineering' or
%       'tensorial', which specifies how shear components are represented.
%
% OUTPUT:
%   B - k-by-l*m matrix, so that each k-by-l block corresponds to a k-by-l
%       matrix A rotated respectively by m-th column of angles.
%
% Examples:
%   tensorRotate([1 0 0 0],pi/4); % ans=[.5 .5 .5 .5]
%   tensorRotate([.5 .5 1],-pi/4,'engineering'); % ans=[1 0 0]

alpha=varargin{1};
if numel(varargin)==2
    form=varargin{2};
end

%a11 a21 a12 a22
if size(A,2)==4
    if isscalar(alpha)
        B=zeros(size(A));
        c2=cos(alpha)^2;
        s2=sin(alpha)^2;
        cs=cos(alpha)*sin(alpha);
        B(:,1)=A(:,1)*c2 - (A(:,2)+A(:,3))*cs + A(:,4)*s2;
        B(:,2)=A(:,2)*c2 + (A(:,1)-A(:,4))*cs - A(:,3)*s2;
        B(:,3)=A(:,3)*c2 + (A(:,1)-A(:,4))*cs - A(:,2)*s2;
        B(:,4)=A(:,4)*c2 + (A(:,3)+A(:,2))*cs + A(:,1)*s2;
    elseif (size(alpha,1)==size(A,1)) && (size(alpha,2)==1)
        B=zeros(size(A));
        c2=cos(alpha).^2;
        s2=sin(alpha).^2;
        cs=cos(alpha).*sin(alpha);
        B(:,1)=A(:,1).*c2 - (A(:,2)+A(:,3)).*cs + A(:,4).*s2;
        B(:,2)=A(:,2).*c2 + (A(:,1)-A(:,4)).*cs - A(:,3).*s2;
        B(:,3)=A(:,3).*c2 + (A(:,1)-A(:,4)).*cs - A(:,2).*s2;
        B(:,4)=A(:,4).*c2 + (A(:,3)+A(:,2)).*cs + A(:,1).*s2;
    elseif (size(alpha,1)==size(A,1)) && (size(alpha,2)>1)
        c2=cos(alpha).^2;
        s2=sin(alpha).^2;
        cs=cos(alpha).*sin(alpha);
        B=zeros(size(A,1),size(A,2)*size(alpha,2));
        for j=1:size(alpha,2)
            B(:,(j-1)*4+1)=...
                A(:,1).*c2(:,j) - (A(:,2)+A(:,3)).*cs(:,j) + A(:,4).*s2(:,j);
            B(:,(j-1)*4+2)=...
                A(:,2).*c2(:,j) + (A(:,1)-A(:,4)).*cs(:,j) - A(:,3).*s2(:,j);
            B(:,(j-1)*4+3)=...
                A(:,3).*c2(:,j) + (A(:,1)-A(:,4)).*cs(:,j) - A(:,2).*s2(:,j);
            B(:,(j-1)*4+4)=...
                A(:,4).*c2(:,j) + (A(:,3)+A(:,2)).*cs(:,j) + A(:,1).*s2(:,j);
        end
    else
        error ('Dimension mismatch')
    end
    
elseif size(A,2)==3
    if numel(varargin)~=2
        error('Additional parameter ''engineering'' or ''tensorial'' must be specified');
    end
    if strcmpi(form,'engineering')
        % e11, e22, g12 --> a11, a21, a12, a22
        A2=zeros(size(A,1),4);
        A2(:,1)=A(:,1);
        A2(:,2)=A(:,3)/2;
        A2(:,3)=A(:,3)/2;
        A2(:,4)=A(:,2);
        B2=rotTensor(A2,alpha);
        B=zeros(size(A,1),size(A,2)*size(alpha,2));
        for j=1:size(alpha,2)
            % a11, a21, a12, a22 --> e11, e22, g12
            B(:,(j-1)*size(A,2)+1)=B2(:,(j-1)*size(A2,2)+1);
            B(:,(j-1)*size(A,2)+2)=B2(:,(j-1)*size(A2,2)+4);
            B(:,(j-1)*size(A,2)+3)=2*B2(:,(j-1)*size(A2,2)+2);%=2*B2(:,3)
        end
    elseif strcmpi(form,'tensorial')
        % e11, e22, e12 --> a11, a21, a12, a22
        A2=zeros(size(A,1),4);
        A2(:,1)=A(:,1);
        A2(:,2)=A(:,3);
        A2(:,3)=A(:,3);
        A2(:,4)=A(:,2);
        B2=rotTensor(A2,alpha);
        B=zeros(size(A,1),size(A,2)*size(alpha,2));
        for j=1:size(alpha,2)
            % a11, a21, a12, a22 --> e11, e22, e12
            B(:,(j-1)*size(A,2)+1)=B2(:,(j-1)*size(A2,2)+1);
            B(:,(j-1)*size(A,2)+2)=B2(:,(j-1)*size(A2,2)+4);
            B(:,(j-1)*size(A,2)+3)=B2(:,(j-1)*size(A2,2)+2);%=B2(:,3)
        end
    else
        error('Parameter ''form'' must be ''engineering'' or ''tensorial'' ');
    end
else
    error('Wrong input format');
end