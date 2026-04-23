function [S,U,V] = resolvent(L0,omega,nEig,W,invW,B,C,filters)
% [gain,U_out,V_in] = ResAnalysis(LHS,RHS,F,omega,noEigs,B,C)
% Performs classical resolvent analysis using "eigs" function of the
% linear operator H=(L0-i*omega)^-1 using input and output matrices
% given by B and C. The adjoint is obtained with respect to the norm W.
% Inputs :
%     L0     : Corresponding time domain operator (with b.c. already
%               imposed)
%     omega  : Freauency (in rads) of the analysis
%     noEigs : number of eigen values desired.
% Optional inputs;
%     W,invW : Matrix containing the metric and its inverse. (default = 1)
%     B,C    : Forcing and observation matrices (default = 1)
%     filters: Structure containing a functions that applies a spatial 
%              filter and its transpose to a vector. It will be to be used 
%              before and after the resolvent operator to filter out junk 
%              modes. Typically obtained from "mesh.filters".
%              (Default=no filter)
% Outputs :
%     S : singular values (sqrt. of gains)
%     U : optimal outputs/responses
%     V : optimal inputs/forcing

% reduce to generalized to standard EVP
nDOF    = size(L0,1);

disp(['    computing resolvent analysis for omega=' num2str(omega)]);

%check inputs
if ~exist('W','var')    ;    W=1; end
if ~exist('invW','var') ; invW=1; end
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



% Correct integration weights shape (if a vector was given, convert to
% diagonal matrix)
if (size(W,1)==1 && size(W,2)==1)
    disp('Using a scalar metric.')
elseif (size(W,1)==1 || size(W,2)==1)
    invW    = 1./W;
    W       = diag(sparse(  W ));
    invW    = diag(sparse(invW));
elseif ~exist('invW') || prod(size(W) == size(invW))==0
    error('W is a matrix and its inverse was not provided or does not has the same size as W.');
    
end

% LU-decomposition of L0-sigma*I
tic
fprintf('  Starting LU-decomposition of L-sigma*I: ');
LsI                 = L0-1i*omega*speye(nDOF);
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

HtrH        = @(v) invW*Htr(W*H(v));

% 'eigs' parameters
opts.tol    = eps;
opts.disp   = 2;
opts.issym  = false;
opts.isreal = false;

% Input via eigendecomposition of H*H
tic

fprintf('  Starting SVD via ''eigs'' accelerated with LU decomposition : ');

[V, OMEGA] 	= eigs(HtrH, nDOF, nEig, 'lm', opts);
S           = sqrt(diag(real(OMEGA))); % singular values of H given by sqrt. of eigenvalues of H*H, should be real (imag. part is O(eps))

time = toc;
disp(['    elapsed time : ' datestr(time/24/3600, 'HH:MM:SS')]);

% Normalization and outputs
U = zeros(size(V));
for i = 1:nEig
    % normalization in inner product norm
    V(:,i)  = V(:,i)/sqrt(V(:,i)'*W*V(:,i));
    U(:,i) =  H(V(:,i))/S(i); %C*(qq*(UU\(LL\(pp*(rr\(B*(V(:,i))))))))/S(i);
end
