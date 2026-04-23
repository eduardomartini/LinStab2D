function [defMesh] = DeformMesh(mesh,X2,Y2,imposeXcoordinatesPeriodicity)
    % newMesh = DeformMesh(mesh,X2,Y2)
    % Deforms the 'mesh' object to the coordinates given by X2 and Y2.
    % Derivative matrices and integration weights are updated accordingly. 

    if ~exist('imposeXcoordinatesPeriodicity') ; imposeXcoordinatesPeriodicity = false;    end
    
    fprintf( '--- Deforming mesh ...  '); tic();
    
    y_symmetry = mesh.y_symmetry;
    
        
    
    defMesh                 = mesh;
    
    % If periodic conditions along x are imposed, we need to remove the linear trend in x before computing the Jacobians. 
    % The corresponding costant is then added a posteriori to the Jacobians. This avoids the discontinuity of x on the boundaries.
    if imposeXcoordinatesPeriodicity && ~isempty(X2)
        x1 = min(X2(:));
        x2 = max(X2(:));
        
        lin_x = linspace(x1,x2,size(X2,2));
        lin_X = repmat(lin_x,size(X2,1),1);
        
        X2 = X2-lin_X;
        
        extra_Dx_X = (x2-x1)/( max(mesh.X(:))-min(mesh.X(:)) );
    else
        lin_X = 0;
        extra_Dx_X = 0;
    end

    mDx  = mesh.Dx;
    mD2x = mesh.D2x;
    mDy  = mesh.Dy;
    mD2y = mesh.D2y;
    
    %% Construct new mesh and prepare diff matrices to be transformed

    ngp = mesh.ngp;
    
    Dx          = mesh.Dx;
    D2x         = mesh.D2x;
    Dy          = mesh.Dy;
    D2y         = mesh.D2y;
    Dxy         = mesh.Dxy;
    Dyx         = mesh.Dyx;
    
    if y_symmetry
        Dy_symm     = mesh.Dy_symm;
        Dy_asymm    = mesh.Dy_asymm;
        D2y_symm    = mesh.D2y_symm;
        D2y_asymm   = mesh.D2y_asymm;
        Dxy_symm    = mesh.Dxy_symm;
        Dxy_asymm   = mesh.Dxy_asymm;
        Dyx_symm    = mesh.Dyx_symm;
        Dyx_asymm   = mesh.Dyx_asymm;
    end

    
    if isempty(Y2)
        defMesh.X               = X2+lin_X;

        X2 = X2(mesh.usedInd);
    
        x2p  = Dx*X2(:)+extra_Dx_X;
        x2pp = D2x*X2(:);
        
        Dx_new  = spdiags(x2p.^-1           ,0,ngp,ngp)*Dx;
        D2x_new = spdiags(x2p.^-2           ,0,ngp,ngp)*D2x - ...
                  spdiags(x2pp.*x2p.^-3     ,0,ngp,ngp)*Dx ;
              
              
        det = abs(x2p); % transformation determinant. Used to update int. weights.

                
        Dxy_new = Dx_new*Dy ;
        Dyx_new = Dy*Dx_new ;
        Dy_new  = Dy;
        D2y_new = D2y;
        
        if y_symmetry
            
            Dxy_symm_new  = Dx_new*Dy_symm ;
            Dyx_symm_new  = Dy_symm*Dx_new ;
            Dy_symm_new   = Dy_symm;
            D2y_symm_new  = D2y_symm;
            
            Dxy_asymm_new = Dx_new*Dy_asymm ;
            Dyx_asymm_new = Dy_asymm*Dx_new ;
            Dy_asymm_new  = Dy_asymm;
            D2y_asymm_new = D2y_asymm;

        
        end
        
    elseif isempty(X2)
        defMesh.Y               = Y2;

        Y2 = Y2(mesh.usedInd);

        y2p  = Dy*Y2(:);
        y2pp = D2y*Y2(:);
        
        Dy_new  = spdiags(y2p.^-1       ,0,ngp,ngp)*Dy;
        D2y_new = spdiags(y2p.^-2       ,0,ngp,ngp)*D2y -  ...
                  spdiags(y2pp.*y2p.^-3 ,0,ngp,ngp)*Dy ;        
        
        Dx_new  = Dx;
        D2x_new = D2x;
          
        Dxy_new = Dx*Dy_new ;
        Dyx_new = Dy_new*Dx ;
              
        det = abs(y2p); % transformation determinant. Used to update int. weights.

        if y_symmetry
            Dy_symm_new  = spdiags(y2p.^-1,0,ngp,ngp)*Dy_symm;
            D2y_symm_new = spdiags(y2p.^-2,0,ngp,ngp)*D2y_symm - ...
                           spdiags(y2pp.*y2p.^-3,0,ngp,ngp)*Dy_symm;        

            Dy_asymm_new  = spdiags(y2p.^-1,0,ngp,ngp)*Dy_asymm;
            D2y_asymm_new = spdiags(y2p.^-2,0,ngp,ngp)*D2y_asymm - ...
                            spdiags(y2pp.*y2p.^-3,0,ngp,ngp)*Dy_asymm;        
            
            
            Dxy_symm_new = Dx*Dy_symm_new ;
            Dyx_symm_new = Dy_symm_new*Dx ;

            Dxy_asymm_new = Dx*Dy_asymm_new ;
            Dyx_asymm_new = Dy_asymm_new*Dx ;
            
        end
    else
        defMesh.X               = X2+lin_X ;
        defMesh.Y               = Y2        ;
        
        X2 = X2(mesh.usedInd);
        Y2 = Y2(mesh.usedInd);

        %% Compute transformation derivatives
        %     Jacobian
        J = zeros(2,2,ngp);
        J(1,1,:) = mDx*X2(:) + extra_Dx_X;
        J(1,2,:) = mDx*Y2(:);

        J(2,1,:) = mDy*X2(:);
        J(2,2,:) = mDy*Y2(:);             
        % 
        dJ = zeros(2,2,ngp,2);
        dJ(1,1,:,1) = mD2x*X2(:);
        dJ(1,2,:,1) = mD2x*Y2(:);

        dJ(2,1,:,1) = mDx*(mDy*X2(:));
        dJ(2,2,:,1) = mDx*(mDy*Y2(:)); 

        dJ(1,1,:,2) = mDy*(mDx*X2(:));
        dJ(1,2,:,2) = mDy*(mDx*Y2(:));

        dJ(2,1,:,2) = mD2y*X2(:);
        dJ(2,2,:,2) = mD2y*Y2(:); 

        % Inverse transform derivative
        Ji  = zeros(2,2,ngp);
        det = J(1,1,:).*J(2,2,:)- J(2,1,:).*J(1,2,:);
        Ji(1,1,:)  =  J(2,2,:)./det;
        Ji(2,1,:)  = -J(2,1,:)./det;
        Ji(1,2,:)  = -J(1,2,:)./det;
        Ji(2,2,:)  =  J(1,1,:)./det;

        %% Compute new derivative matrices
        %% First Derivatives
        % f(xi)=f(xi(sigj))
        % df/dxi = df/dsigj dsigj/dxi  =  df/dsigj (dxj/dsigi)^-1  = Ji_ij * df/dsigj

        Dx_new          = spdiags(squeeze(Ji(1,1,:)),0,ngp,ngp)*Dx ...
                        + spdiags(squeeze(Ji(1,2,:)),0,ngp,ngp)*Dy;
        Dy_new          = spdiags(squeeze(Ji(2,1,:)),0,ngp,ngp)*Dx ...
                        + spdiags(squeeze(Ji(2,2,:)),0,ngp,ngp)*Dy;

        if y_symmetry
            Dy_symm_new     = spdiags(squeeze(Ji(2,1,:)),0,ngp,ngp)*Dx ...
                            + spdiags(squeeze(Ji(2,2,:)),0,ngp,ngp)*Dy_symm;
            Dy_asymm_new    = spdiags(squeeze(Ji(2,1,:)),0,ngp,ngp)*Dx ...
                            + spdiags(squeeze(Ji(2,2,:)),0,ngp,ngp)*Dy_asymm;
        end

        %% Second Derivatives
        % ddf/dxidxj = Ji_ik  ddf/dsigkdsigl Ji_jl  +         d(Ji_ik)/dxj   * df/dsigk 
        %            = Ji_ik  ddf/dsigkdsigl Ji_jl  + Ji_jl * d(Ji_ik)/dsigl * df/dsigk 
        %            = Ji_ik  ddf/dsigkdsigl Ji_jl  + Q_ijk                  * df/dsigk 

        JiT = Ji;
        JiT(1,2,:) = Ji(2,1,:);
        JiT(2,1,:) = Ji(1,2,:);
        for i=1:4
            A=zeros(2,2); A(i)=1;
            a(:,:,:,i)  =  multiprod(JiT,multiprod(A, Ji));
        end
        axx=squeeze(a(1,1,:,:));
        axy=squeeze(a(1,2,:,:));
        ayx=squeeze(a(2,1,:,:));
        ayy=squeeze(a(2,2,:,:));

        % d(Ji)/dsigl = - Ji*dJ/dsigl*Ji 

        dJi(:,:,:,1)=-multiprod(multiprod(Ji,dJ(:,:,:,1)),Ji);
        dJi(:,:,:,2)=-multiprod(multiprod(Ji,dJ(:,:,:,2)),Ji);

        % Q_ijk = Ji_jl * d(Ji_ik)/dsigl
        Q = zeros(2,2,2,ngp);
        for i=1:2; for j=1:2; for k=1:2
                    Q(i,j,k,:) = Ji(j,1,:).*dJi(i,k,:,1) + Ji(j,2,:).*dJi(i,k,:,2);
        end; end; end

        Dxy = Dx*Dy;
        Dyx = Dy*Dx;

        diags = @(x)spdiags(x,0,ngp,ngp);

        D2x_new =  diags(axx(:,1))*D2x + diags(axy(:,1))*Dxy + ...
                diags(ayx(:,1))*Dyx + diags(ayy(:,1))*D2y + ...
                diags(squeeze(Q(1,1,1,:)))*Dx + diags(squeeze(Q(1,1,2,:)))*Dy ;
        Dxy_new =  diags(axx(:,2))*D2x + diags(axy(:,2))*Dxy + ...
                diags(ayx(:,2))*Dyx + diags(ayy(:,2))*D2y + ...
                diags(squeeze(Q(1,2,1,:)))*Dx + diags(squeeze(Q(1,2,2,:)))*Dy ;
        Dyx_new =  diags(axx(:,3))*D2x + diags(axy(:,3))*Dxy + ...
                diags(ayx(:,3))*Dyx + diags(ayy(:,3))*D2y + ...
                diags(squeeze(Q(2,1,1,:)))*Dx + diags(squeeze(Q(2,1,2,:)))*Dy ;
        D2y_new =  diags(axx(:,4))*D2x + diags(axy(:,4))*Dxy + ...
                diags(ayx(:,4))*Dyx + diags(ayy(:,4))*D2y + ...
                diags(squeeze(Q(2,2,1,:)))*Dx + diags(squeeze(Q(2,2,2,:)))*Dy ;

        if y_symmetry

            Dxy_symm_new =  diags(axx(:,2))*D2x            + diags(axy(:,2))*Dxy_symm + ...
                         diags(ayx(:,2))*Dyx_symm       + diags(ayy(:,2))*D2y_symm + ...
                         diags(squeeze(Q(1,2,1,:)))*Dx  + diags(squeeze(Q(1,2,2,:)))*Dy_symm ;
            Dyx_symm_new =  diags(axx(:,3))*D2x            + diags(axy(:,3))*Dxy_symm + ...
                         diags(ayx(:,3))*Dyx_symm       + diags(ayy(:,3))*D2y_symm + ...
                         diags(squeeze(Q(2,1,1,:)))*Dx  + diags(squeeze(Q(2,1,2,:)))*Dy_symm ;
            D2y_symm_new =  diags(axx(:,4))*D2x            + diags(axy(:,4))*Dxy_symm + ...
                         diags(ayx(:,4))*Dyx_symm       + diags(ayy(:,4))*D2y_symm + ...
                         diags(squeeze(Q(2,2,1,:)))*Dx  + diags(squeeze(Q(2,2,2,:)))*Dy_symm ;

            Dxy_asymm_new=  diags(axx(:,2))*D2x            + diags(axy(:,2))*Dxy_asymm + ...
                         diags(ayx(:,2))*Dyx_asymm       + diags(ayy(:,2))*D2y_asymm + ...
                         diags(squeeze(Q(1,2,1,:)))*Dx  + diags(squeeze(Q(1,2,2,:)))*Dy_asymm ;
            Dyx_asymm_new=  diags(axx(:,3))*D2x            + diags(axy(:,3))*Dxy_asymm + ...
                         diags(ayx(:,3))*Dyx_asymm       + diags(ayy(:,3))*D2y_asymm + ...
                         diags(squeeze(Q(2,1,1,:)))*Dx  + diags(squeeze(Q(2,1,2,:)))*Dy_asymm ;
            D2y_asymm_new=  diags(axx(:,4))*D2x            + diags(axy(:,4))*Dxy_asymm + ...
                         diags(ayx(:,4))*Dyx_asymm       + diags(ayy(:,4))*D2y_asymm + ...
                         diags(squeeze(Q(2,2,1,:)))*Dx  + diags(squeeze(Q(2,2,2,:)))*Dy_asymm ;
        end       
    end            
    % update intregration weights
    Wnew        = zeros(size(mesh.W));
    Wnew(mesh.usedInd)       = mesh.W(mesh.usedInd).*abs(det(:));
    
    if y_symmetry
        Wnew_symm   = zeros(size(mesh.W));
        Wnew_asymm  = zeros(size(mesh.W));
        Wnew_symm(mesh.usedInd)  = mesh.W_symm(mesh.usedInd).*abs(det(:));
        Wnew_asymm(mesh.usedInd) = mesh.W_asymm(mesh.usedInd).*abs(det(:));
    end
            
    defMesh.Dx        = Dx_new ;
    defMesh.Dy        = Dy_new ;

    defMesh.D2x       = D2x_new ;
    defMesh.Dxy       = Dxy_new ;
    defMesh.Dyx       = Dyx_new ;
    defMesh.D2y       = D2y_new ;
    
    defMesh.W         = Wnew;
    
    if y_symmetry
        defMesh.Dy_symm   = Dy_symm_new ;
        defMesh.Dy_asymm  = Dy_asymm_new ;

        defMesh.Dxy_symm  = Dxy_symm_new ;
        defMesh.Dyx_symm  = Dyx_symm_new ;
        defMesh.D2y_symm  = D2y_symm_new ;

        defMesh.Dxy_asymm = Dxy_asymm_new ;
        defMesh.Dyx_asymm = Dyx_asymm_new ;
        defMesh.D2y_asymm = D2y_asymm_new ;

        defMesh.W_symm    = Wnew_symm;
        defMesh.W_asymm   = Wnew_asymm;
    end
    
fprintf( ' Done in %.0f seconds.\n',toc); 
