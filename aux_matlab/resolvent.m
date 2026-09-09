function [S,Uout,Vout] = resolvent(L0,omega,nEig,Wq,invWf,B,C,filters,Wf,returnPhysicalSpace)
% [gain,U,V] = resolvent(L0,omega,nEig,Wq,invWf,B,C,filters,Wf,returnPhysicalSpace)
% Performs classical resolvent analysis using "eigs" function of the
% linear operator H=(L0-i*omega)^-1 using input and output matrices
% given by B and C. The adjoint is obtained with respect to the norm W.
% Inputs :
%     L0     : Corresponding time domain operator (with b.c. already
%               imposed)
%     omega  : Freauency (in rads) of the analysis
%     noEigs : number of eigen values desired.
% Optional inputs;
%     Wq       : Matrix containing input space  metric. (default = 1)
%     invWf    : Matrix with the the inverse of the output space metrics
%                 (defaults to 1)
%     B,C    : Forcing and observation matrices (default = 1)
%     filters: Structure containing a functions that applies a spatial 
%              filter and its transpose to a vector. It will be to be used 
%              before and after the resolvent operator to filter out junk 
%              modes. Typically obtained from "mesh.filters".
%              (Default=no filter)
%
%     Wf       : Matrix with the output space metrics. Used only for output
%                vetor normalization (default = Wq, if B is a square matrix)
%     returnPhysicalSpace : if true, forces and resposes are converted into
%                           the physical space, i.e., forces are multiplied 
%                           by B, and its respective responses, before 
%                           multiplication by C, are returned.
%     
% Outputs :
%     S : singular values (sqrt. of gains)
%     U : optimal outputs/responses
%     V : optimal inputs/forcing


% Check inputs
if ~exist('Wq','var')    ; Wq=1; end

if ~exist('invWf','var') ; invWf=1; end

if ~exist('Wf','var')    
    if size(B,1)==size(B,2) && size(C,1)==size(C,2) && size(B,1)==size(C,1)
        Wf=Wq; 
    else
        Wf=1 ;
    end
end

if ~exist('returnPhysicalSpace','var') ; returnPhysicalSpace=true; end

if ~exist('B','var')    
    B=1; 
elseif size(B,1)==1 || size(B,2)==1 
    B = spdiags(B,0,length(B),length(B)); 
end

if ~exist('C','var')    
    C=1; 
elseif size(C,1)==1 || size(C,2)==1 
    C = spdiags(C,0,length(C),length(C)); 
end


% Correct integration weights shape: if a vector was given, convert to diagonal matrix
if (size(invWf,1)==1 && size(invWf,2)==1)
    disp('Using a scalar invWf.')
elseif  prod(size(invWf) == size(B,2))==0
    error('invWf is a matrix, but its size is not compatible with the number of columns of B.');    
end

if (size(Wf,1)==1 && size(Wf,2)==1)
    disp('Using a scalar Wf.')
elseif  prod(size(Wf) == size(B,2))==0
    error('Wf is a matrix, but its size is not compatible with the number of columns of B.');    
end

if (size(Wq,1)==1 && size(Wq,2)==1)
    disp('Using a scalar invWq.')
elseif  prod(size(Wq) == size(C,1))==0
    error('Wq is a matrix, but its size is not compatible with the number of rows of C.');    
end


% Check matrix dimensions
size_L0 = size(L0);
size_B  = size(B);
size_C  = size(C);

if size_L0(1) ~= size_L0(2)
    error('L0 matrix is not square.')
end

nDOFs   = size_L0(1);
nForces = size_B(2);


% B has to be either a matrix with the same number of columns as L0, or a scalar 
if size_B(1) ~= size_L0(2) && prod(size_B) ~= 1
    error('B matrix has incompatible number of rows with L0.')
end
% C has to be either a matrix with the same number of rows as L0, or a scalar 
if size_C(2) ~= size_L0(1) && prod(size_C) ~= 1
    error('C matrix has incompatible number of columns with L0.')
end

%%  
disp(['    computing resolvent analysis for omega=' num2str(omega)]);


% LU-decomposition of L0-sigma*I
tic;

fprintf('  Starting LU-decomposition of L-sigma*I: ');
LsI      =  L0-1i*omega*speye(nDOFs);
[R,R_T]  =  GetInverseFunction(LsI);
time = toc;

disp(['    elapsed time - LU-decomposition of L-sigma*I: ' datestr(time/24/3600, 'HH:MM:SS')]);

% Create spatial filter function from filters in the mesh object
if exist('filters','var')
    filter    = filters.filter;
    filter_ct = filters.filter_ct;
    FILTER    = @(x) reshape( filter   ( reshape(x,[],5)),[],1) ; 
    FILTER_ct = @(x) reshape( filter_ct( reshape(x,[],5)),[],1) ; 
else
    FILTER       = @(x)x;
    FILTER_ct    = @(x)x;
end

% Function handle for resolvent operator H and Hermitian H*H
H           = @(v) C *FILTER(   R(  FILTER(   B *v)));
Htr         = @(v) B'*FILTER_ct(R_T(FILTER_ct(C'*v)));

HtrH        = @(v) invWf*Htr(Wq*H(v));

% 'eigs' parameters
opts.tol    = eps;
opts.disp   = 2;
opts.issym  = false;
opts.isreal = false;

% Input via eigendecomposition of H*H
tic

fprintf('  Starting SVD via ''eigs'' accelerated with LU decomposition : ');

[V, OMEGA] 	= eigs(HtrH, nForces, nEig, 'lm', opts);
S           = sqrt(diag(real(OMEGA))); % singular values of H given by sqrt. of eigenvalues of H*H, should be real (imag. part is O(eps))

time = toc;
disp(['    elapsed time : ' datestr(time/24/3600, 'HH:MM:SS')]);

if returnPhysicalSpace
    % Normalization and outputs
    Vout = B*V; % optimal inputs/forcing in physical space (after applying B)    
    Uout = zeros(size(Vout));
    for i = 1:nEig
        % Re scale modes for norm 1 force
        v = V(:,i);
        scale = sqrt(v'*(Wf*v));
        Vout(:,i) = Vout(:,i)/scale;

        % reponse from the forcing in filter space (before multiplication 
        % by C)
        Uout(:,i) =  FILTER(   R(  FILTER(   B *v))) /S(i)/scale; 
    end
else
    Vout = V;
    Uout = zeros(size_C(1),nEig);
    for i = 1:nEig
        % Re scale modes for norm 1 force
        v = V(:,i);
        scale = sqrt(v'*(Wf*v));
        Vout(:,i) = Vout(:,i)/scale;

        % reponse from the forcing in output space (after multiplication 
        % by C)
        Uout(:,i) =  H(V(:,i))/S(i)/scale; 
    end
end

