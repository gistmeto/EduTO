function [U,K] = fun_fem_solve(fem,opt)

% Initializing FEM Variable
    K=sparse(fem.ndof,fem.ndof);
    F=zeros(fem.ndof,1);

% Build K Matirx 
    kelist=kron(fem.Emin+(fem.E0-fem.Emin)*opt.erho.^opt.penal,ones(144,1))'.*fem.K_S;    
    K=sparse(fem.is,fem.js,kelist); K=(K+K')/2;

% Apply boundary conditions 
    K(fem.bcdof,:)=0;
    K(:,fem.bcdof)=0;
    for i=1:length(fem.bcdof)
        K(fem.bcdof(i),fem.bcdof(i))=1;
    end
    F(fem.bcdof)=fem.bcval;

% Apply traction force
    for i = 1:size(fem.fdof,1)
        F(fem.fdof(i)) = F(fem.fdof(i)) + fem.fval(i);
    end

% Solve 
    U = K\F;
end

