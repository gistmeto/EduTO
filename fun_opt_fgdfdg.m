function [f,g,dfdx,dgdx] = fun_opt_fgdfdg(fem,opt)

% Evaluation of objective (f) and constraint (g) functions
    f = fem.U'*fem.K*fem.U; 
    g = (fem.Ve*opt.erho)/opt.VT-opt.volfrac;    

% Derivatives    
    DerhoDnrho = opt.Ten';
    DnrhoDfdv = (sech(opt.bt*opt.fdv)).^2*opt.bt/(2*tanh(opt.bt)); % TODO
    
% Evaluation of dfdx
    dkelist = kron((fem.E0-fem.Emin)*opt.penal*opt.erho.^(opt.penal-1),ones(144,1))'.*fem.K_S;
    ulist1=reshape(kron(fem.U(fem.edof)',ones(1,12)),fem.ne*144,1)';
    ulist2=reshape(kron(fem.U(fem.edof)',ones(12,1)),fem.ne*144,1)';
    DfDerho=-sum(reshape((ulist1.*dkelist.*ulist2),144,fem.ne))';
    DfDdv=(DerhoDnrho*DfDerho).*DnrhoDfdv;
    dfdx = DfDdv(opt.dof_dd);
    
% Evaluation of dgdx
    DgDerho=fem.Ve'/opt.VT;    
    DgDdv=(DerhoDnrho*DgDerho).*DnrhoDfdv;
    dgdx= DgDdv(opt.dof_dd);
    dgdx=dgdx';
        
end
