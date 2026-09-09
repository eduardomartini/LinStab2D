function  mesh = AddDW2mesh(mesh)
%CREATEDW Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('FDorder')
        if isfield(mesh,'FDorder')
            useStandard_DF = mesh.FDorder<0;
            useSBP_FD      = mesh.FDorder>0;
            FDorder = abs(mesh.FDorder);
        else
            error('CreateDW : FD order not defined as an argument nor specified in the mesh object');
        end
    end

    if isfield(mesh,'x_periodicity')
        periodic=mesh.x_periodicity;
    else
        periodic=false;        
    end
    if periodic
            %extend the domain FDorder points on both sides, keeping mask and
            %all
            [NY,NX] = size(mesh.X);
            mesh_bak = mesh;

            dx = mesh.X(2,2)-mesh.X(1,1);
            mesh.X = [mesh.X(:,end-FDorder+1:end)-NX*dx,mesh.X,mesh.X(:,1:FDorder)+NX*dx];
            mesh.Y = [mesh.Y(:,end-FDorder+1:end)      ,mesh.Y,mesh.Y(:,1:FDorder)      ];

            mesh.ngp = numel(mesh.X);
            mesh.usedInd = 1:mesh.ngp;
    end
    
    N=mesh.ngp;
    [NY,NX] = size(mesh.X);

    Dx              = sparse(N,N);
    D2x             = sparse(N,N);
    Dy              = sparse(N,N);
    D2y             = sparse(N,N);
    W               = nan (NY,NX); 
    W(mesh.usedInd) = 1;
    
    if mesh.y_symmetry
        Dy_symm  = sparse(N,N);
        D2y_symm = sparse(N,N);
        Dy_asymm = sparse(N,N);
        D2y_asymm= sparse(N,N);

        W_symm  = zeros (NY,NX); 
        W_asymm = zeros (NY,NX); 
        W_symm (mesh.usedInd)=1;
        W_asymm(mesh.usedInd)=1;
    end
    
    X = mesh.X(mesh.usedInd);
    Y = mesh.Y(mesh.usedInd);
   
    xList=sort(unique(X(:)));
    yList=sort(unique(Y(:)));
    
    dx = mesh.X(2,2)-mesh.X(1,1);
    dy = mesh.Y(2,2)-mesh.Y(1,1);
    
    ySym = dy/2;
    
    % Create X derivative Matrices
    for y = yList'
        % find all used points with a given x coordinate. 
        pos = find(Y==y);
        
        %guarantee that points are in increasing order
        [~,order] = sort(X(pos));
        pos=pos(order);
        
        %find holes 
        xl=X(pos);
       
        dx_line=(xl(2:end)-xl(1:end-1))/dx;
        domList = [0;find(dx_line>1.5);length(xl)];
        
        %create a FD scheme in each domain
        for idom = 1:length(domList)-1
           
            dom = domList(idom)+1:domList(idom+1);
            NX  = length(dom);
            
            dom = domList(idom)+1:domList(idom+1);
            if NX<2
                error('CreateDW : Cannot construct derivative with less than three points!')
            elseif NX<FDorder+1
                FDorder_curr = 2;
            else
                FDorder_curr=FDorder;
            end
            
             if useStandard_DF
                [Dx_1D ,~] = getNonCompactFDmatrix(NX,dx,1,FDorder_curr);
                [D2x_1D,~] = getNonCompactFDmatrix(NX,dx,2,FDorder_curr);
            end
            
            if useSBP_FD                
                [   Dx_1D     ,D2x_1D     , ...
                    ~         ,     ~     , ...
                    ~         ,     ~      ] = Dmats_SBP(NX,dx,FDorder_curr);
            end
        
            dom_pos = pos(dom);
            Dx      (dom_pos,dom_pos) = Dx_1D        ;
            D2x     (dom_pos,dom_pos) = D2x_1D       ;
          
            
            %trapeizodal rule
            
            meshpos=mesh.usedInd(dom_pos);
            W     (meshpos)          = W     (meshpos)*dx;
            W     (meshpos([1,end])) = W     (meshpos([1,end]))/2;
            
        end
    end
    
    % Create Y derivative Matrices

    for x = xList'
        % find all used points with a given x coordinate. 
        pos = find(X==x);
        
        %guarantee that points are in increasing order
        [~,order] = sort(Y(pos));
        pos=pos(order);
        
        %find holes 
        yl=Y(pos);
        dy_line=(yl(2:end)-yl(1:end-1))/dy;
        domList = [0;find(dy_line>1.5);length(yl)];
        
        %create a FD scheme in each domain
        for idom = 1:length(domList)-1
           
            dom = domList(idom)+1:domList(idom+1);
            Ny  = length(dom);

            if Ny<2
                error('CreateDW : Cannot construct derivative with less than three points!')
            elseif Ny<=FDorder+1
                FDorder_curr = 2;
            else
                FDorder_curr=FDorder;
            end
                
            if useStandard_DF
                [Dy_1D ,~] = getNonCompactFDmatrix(Ny,dy,1,FDorder_curr);
                [D2y_1D,~] = getNonCompactFDmatrix(Ny,dy,2,FDorder_curr);

                [Dy_1D_ext ,~] = getNonCompactFDmatrix(Ny*2,dy,1,FDorder_curr);
                [D2y_1D_ext,~] = getNonCompactFDmatrix(Ny*2,dy,2,FDorder_curr);
                
                flip2= @(A)flip(A,2);
                Dy_1D_sym   = Dy_1D_ext (end/2+1:end,end/2+1:end)   + flip2(Dy_1D_ext (end/2+1:end,1:end/2));
                D2y_1D_sym  = D2y_1D_ext(end/2+1:end,end/2+1:end)   + flip2(D2y_1D_ext(end/2+1:end,1:end/2));

                Dy_1D_asym  = Dy_1D_ext(end/2+1:end,end/2+1:end)    - flip2(Dy_1D_ext (end/2+1:end,1:end/2));
                D2y_1D_asym = D2y_1D_ext(end/2+1:end,end/2+1:end)   - flip2(D2y_1D_ext(end/2+1:end,1:end/2));
            end
            
            if useSBP_FD
                [   Dy_1D      , D2y_1D     , ...
                    Dy_1D_sym  , Dy_1D_asym , ...
                    D2y_1D_sym , D2y_1D_asym ] = Dmats_SBP(Ny,dy,FDorder_curr);
            end
            
            dom_pos = pos(dom);
            
            Dy      (dom_pos,dom_pos) = Dy_1D        ;
            D2y     (dom_pos,dom_pos) = D2y_1D       ;

            %trapeizodal rule
            meshpos=mesh.usedInd(dom_pos);
            W(meshpos)               = W(meshpos)*dy;
            W(meshpos([1,end]))      = W(meshpos([1,end]))/2;
            
            if (abs(yl(1)-ySym)<.1*dy) && mesh.y_symmetry
                Dy_symm  (dom_pos,dom_pos) = Dy_1D_sym    ;
                D2y_symm (dom_pos,dom_pos) = D2y_1D_sym   ;
                Dy_asymm (dom_pos,dom_pos) = Dy_1D_asym   ;
                D2y_asymm(dom_pos,dom_pos) = D2y_1D_asym  ;
                
                W_symm (meshpos)           = W_symm (meshpos)*dx;
                W_symm (meshpos(end))      = W_symm (meshpos(end))/2;
                W_asymm(meshpos)           = W_asymm(meshpos)*dx;
                W_asymm(meshpos( 1 ))      = W_asymm(meshpos( 1 )).*3/4;
                W_asymm(meshpos(end))      = W_asymm(meshpos(end))  ./2;
            end
        end
    end
    
    

    if periodic
        %impose periodicity and remove peridic extension
        xlim=[mesh_bak.X(1,1),mesh_bak.X(end,end)];    

        dom = find( mesh.X>=xlim(1) & mesh.X<=xlim(2) );
        
        Nx=NX-FDorder*2;
        %Create matrix to impose periodicity
        B = sparse(Nx,Nx+FDorder*2);
        for i=1:Nx
            B(i,i+FDorder)=1;
        end
        for i=1:FDorder
            B(end+(i-1)-(FDorder-1),i)=1;
            B(i,end+(i-1)-(FDorder-1))=1;

        end 
        BB = sparse(kron(B,speye((Ny))));
        Dx  = Dx (dom,:)*BB';
        D2x = D2x(dom,:)*BB';
%         Dx  = Dx(dom,dom);
%         D2x = D2x(dom,dom);
        Dy          = Dy (dom,dom);
        D2y         = D2y(dom,dom);
        
        W           = W(dom);
        
        [NY,NX]=size(mesh_bak.X);

        W       = reshape(W         ,NY,NX);

        if mesh.y_symmetry
            
            Dy_symm     = Dy_symm(dom,dom);
            D2y_symm    = D2y_symm(dom,dom);
            Dy_asymm    = Dy_asymm(dom,dom);
            D2y_asymm   = D2y_asymm(dom,dom);
        
            W_symm      = W_symm(dom);
            W_asymm     = W_asymm(dom);
            
            W_symm  = reshape(W_symm    ,NY,NX);
            W_asymm = reshape(W_asymm   ,NY,NX);
        end
        mesh.X       = reshape(mesh.X(dom),NY,NX);
        mesh.Y       = reshape(mesh.Y(dom),NY,NX);
        mesh.ngp     = length(dom);
        mesh.usedInd = (1:mesh.ngp)';
    end

    %finalize with mixed second derivatives
    mesh.Dxy = Dx*Dy ; 
    mesh.Dyx = Dy*Dx ; 
    
    mesh.Dx  = Dx ; 
    mesh.D2x = D2x; 

    mesh.Dy  = Dy ;
    mesh.D2y = D2y; 
    
    mesh.W   = W ; 

    
    if mesh.y_symmetry
        mesh.Dy_symm  = Dy_symm ;
        mesh.D2y_symm = D2y_symm; 

        mesh.Dy_asymm  = Dy_asymm ;
        mesh.D2y_asymm = D2y_asymm; 

        mesh.Dxy_asymm = Dy_asymm*Dx; 
        mesh.Dyx_asymm = Dx*Dy_asymm; 

        mesh.Dxy_symm = Dy_symm*Dx; 
        mesh.Dyx_symm = Dx*Dy_symm; 

        mesh.W_symm    = W_symm ; 
        mesh.W_asymm   = W_asymm ; 
    end
end