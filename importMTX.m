function A=importMTX(filename,varargin)
% Import data contained in a file as a sparse matrix. Abaqus .mtx files is
% an example.
% INPUT:
%   filename - string, file name.
%   numEq (optional) - number of degrees of freedom defined at each node.
% OUTPUT:
%   A - sparse matrix imported from the file provided.
%
% Disclaimer: The author is not associated with owners of Abaqus
% trademark in any way.

    sparseFormat=dlmread(filename);

    if ~isempty(varargin)
        numEq=varargin{1};
    else
        numEq=max(sparseFormat(:,2));
    end
    
    %convert to our numbering convention
    IDs1=numEq*(sparseFormat(:,1)-1)+sparseFormat(:,2); 
    IDs2=numEq*(sparseFormat(:,3)-1)+sparseFormat(:,4);
    
    %determine the size of the matrix
    numGDoF=max([IDs1; IDs2]);

    %reflect the matrix to the other triangle 
    IDs1_=IDs1;
    IDs2_=IDs2;
    Vals_=sparseFormat(:,5);
    
    %remove the diagonal
    IDs1_(IDs1==IDs2)=[];
    IDs2_(IDs1==IDs2)=[];
    Vals_(IDs1==IDs2)=[];
    
    A=sparse([IDs1; IDs2_],[IDs2; IDs1_],[sparseFormat(:,5); Vals_],numGDoF,numGDoF);
end