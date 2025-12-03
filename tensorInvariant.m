function invar=tensorInvariant(A,str)
% Compute an invariant of a second-order tensor.
% INPUT:
%   A - k-by-4 matrix, k-th row represents a 2d tensor as A11, A21,
%       A12, A22;
%       or k-by-9 matrix, k-th row represents a 3d tensor as A11, A21,
%       A31, A12, A22, A32, A13, A23, A33;
%       or k-by-3 matrix, k-th row represents a symmetric 2d tensor as A11,
%       A22, A12;
%       or k-by-6 matrix, k-th row represents a symmetric 3d tensor as A11,
%       A22, A33, A23, A13, A12.
%   str - a string identifier for the requested quantity: 'I1', 'I2', 'I3',
%       'J1', 'J2', 'J3', 'mises', 'tresca', 'eigen'.
%
% OUTPUT:
% invar - k-by-1 matrix, k-th element is the requested invariant for the
%   k-th input tensor;
%   or k-by-2 or k-by-3 matrix, k-th row contains eigenvalues for the k-th
%   input tensor.
%
% Remark: J2 is computed as -I2 of the deviatoric part. 
% Von Mises stress for 2D is consistent with 3D plane stress setting. 
% On 'eigen', eigenvalues are returned.
%
% Examples:
% I1=tensorInvariant(A,str,'I1')

k=size(A,1);
if size(A,2)==3
    B=zeros(k,4);
    B(:,1)=A(:,1); %A11
    B(:,2)=A(:,3); %A21
    B(:,3)=A(:,3); %A12
    B(:,4)=A(:,2); %A22
    is3d=0;
elseif size(A,2)==4
    B=A;
    is3d=0;
elseif size(A,2)==9
    B=A;
    is3d=1;
elseif size(A,2)==6
    B=zeros(k,9);
    B(:,1)=A(:,1); %A11
    B(:,2)=A(:,6); %A21
    B(:,3)=A(:,5); %A31
    B(:,4)=A(:,6); %A12
    B(:,5)=A(:,2); %A22
    B(:,6)=A(:,4); %A32
    B(:,7)=A(:,5); %A13
    B(:,8)=A(:,4); %A23
    B(:,9)=A(:,3); %A33
    is3d=1;
else
    error('Dimension mismatch');
end
 
if strcmp(str,'I1')
    if ~is3d
        invar=B(:,1)+B(:,4);
    else
        invar=B(:,1)+B(:,5)+B(:,9);
    end
elseif strcmp(str,'I2')
    if ~is3d
        invar=B(:,1).*B(:,4)-B(:,2).*B(:,3);
    else
        invar=B(:,1).*B(:,5) + B(:,1).*B(:,9) + B(:,9).*B(:,5) ...
            - (B(:,8).*B(:,6) + B(:,7).*B(:,3) + B(:,4).*B(:,2));
    end
elseif strcmp(str,'I3')
    if ~is3d
        invar=B(:,1).*B(:,4)-B(:,2).*B(:,3);
        warning('I3 requested for a 2d tensor. I2 (determinant) is returned instead');
    else
        invar=B(:,1).*B(:,5).*B(:,9) + B(:,4).*B(:,8).*B(:,3) ...
            + B(:,2).*B(:,6).*B(:,7) ...
            - B(:,3).*B(:,5).*B(:,7)  - B(:,2).*B(:,4).*B(:,9) ...
            - B(:,1).*B(:,8).*B(:,6);
    end
elseif strcmp(str,'J1')
        invar=0;
elseif strcmp(str,'J2') %equals to -I2 of deviatoric part (convention)
    if ~is3d
        invar= 1/4 * (B(:,1)-B(:,4)).^2 + B(:,3).*B(:,4);
    else
        invar= 1/6* ( (B(:,1)-B(:,5)).^2 + (B(:,1)-B(:,9)).^2 + (B(:,5)-B(:,9)).^2 )...
            + B(:,8).*B(:,6) + B(:,4).*B(:,2) + B(:,3).*B(:,7);
    end
elseif strcmp(str,'J3')
    if ~is3d
        invar= 1/4 * (B(:,1)-B(:,4)).^2 + B(:,3).*B(:,4);
        warning('J3 requested for a 2d tensor. J2=-I2(dev(A)) is returned instead');
    else
        invar = 2/27*(B(:,1).^3+B(:,5).^3+B(:,9).^3)...
            + B(:,4).*B(:,8).*B(:,3) + B(:,7).*B(:,6).*B(:,2)...
            + 1/3*B(:,1) .* ( B(:,4).*B(:,2) + B(:,7).*B(:,3) )...
            + 1/3*B(:,5) .* ( B(:,4).*B(:,2) + B(:,8).*B(:,6) )...
            + 1/3*B(:,9) .* ( B(:,7).*B(:,3) + B(:,8).*B(:,6) )...
            -2/3*(B(:,3).*B(:,5).*B(:,7)+B(:,2).*B(:,4).*B(:,9)+B(:,1).*B(:,8).*B(:,6))...
            +4/9*B(:,1).*B(:,5).*B(:,9)...
            -1/3*( B(:,1).^2.*B(:,5) + B(:,5).^2.*B(:,1) + ...
                   B(:,5).^2.*B(:,9) + B(:,9).^2.*B(:,5) + ...
                   B(:,9).^2.*B(:,1) + B(:,1).^2.*B(:,9) );
    end
elseif strcmpi(str,'tresca')
    if ~is3d
        invar=sqrt( 1/4 * (B(:,1)-B(:,4)).^2 + B(:,2).*B(:,3) );
    else
        % should be reviewed in the future
        evs=zeros(k,3);
        for i=1:k
            evs(i,:)=eig(reshape(B(i,:),3,3));
        end
        invar=max(abs(evs(:,[1 2 3])-evs(:,[2 3 1])),[],2);
    end
elseif any(strcmpi(str,{'mises','vonmises','von mises'}))
    if ~is3d
        invar=sqrt( 1/2 * ((B(:,1)-B(:,4)).^2 + B(:,1).^2 + B(:,4).^2) ...
            + 3*B(:,2).*B(:,3) );
    else
        invar=sqrt( 1/2 * ((B(:,1)-B(:,5)).^2 + (B(:,9)-B(:,5)).^2 + (B(:,1)-B(:,9)).^2) ...
            + 3*B(:,2).*B(:,4) + 3*B(:,3).*B(:,7) + 3*B(:,6).*B(:,8));
    end
elseif any(strcmpi(str,{'eigen','eig'}))
    if ~is3d
        invar=zeros(k,2);
        for i=1:k
            invar(i,:)=eig(reshape(B(i,:),2,2));
        end
    else
        invar=zeros(k,3);
        for i=1:k
            invar(i,:)=eig(reshape(B(i,:),3,3));
        end
    end
else
    error('Unknown parameter value %s',str);
end
    