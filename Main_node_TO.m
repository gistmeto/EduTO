clear; clc; close all; tic;
modelname = 'Cyllinder2_Beam';

%% Pre-Processing
    [inputs] = fun_pre_inputsload(modelname);      % Load input parameters
    [msh]    = fun_pre_mshload(modelname);         % Load mesh information
    [fem]    = fun_pre_feminit(inputs,msh);        % Initialization for FEM
    [opt]    = fun_pre_optinit(inputs,fem);        % Initialization for Optimization
    fprintf('Elapsed time for Pre-Processing:%s\n',(datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')))
    
%% Main processing - Topology optimization 
    while opt.bt < inputs.bt_fn
        % Three-field projection and element density calculation
        opt.fdv=opt.chol_Kft'\(opt.chol_Kft\(opt.Tft*opt.nv));  
        opt.nrho = max(min(tanh(opt.bt*opt.fdv)/(2*tanh(opt.bt))+0.5,1),-1); 
        opt.erho=opt.Ten*opt.nrho;                          

        % FEM, Sensitivity analysis, and objective error
        [fem.U,fem.K] = fun_fem_solve(fem,opt);                  
        [f,g,dfdx,dgdx]= fun_opt_fgdfdg(fem,opt);                
        opt.fhis(opt.iter)=f; opt.ghis(opt.iter)=g;             
        if (opt.iter>1) opt.deltaf =...                        
           abs((opt.fhis(opt.iter)-opt.fhis(opt.iter-1))/opt.fhis(opt.iter-1));end    

        % Display and plotting
        disp(sprintf(['Iter:%d, f:%.4f, Volume:%.4f, deltaf:%.5f, beta:%.2f'] ...
        , opt.iter, f, g+opt.volfrac, opt.deltaf,opt.bt));
        fun_opt_plot(fem,opt);        

        % Update design varialbe using MMA optimizer
        [dvnew,~,~,~,~,~,~,~,~,opt.MMA.low,opt.MMA.upp] = ...
            mmasub(1,length(opt.dv),opt.iter,opt.dv,opt.dvmin,opt.dvmax,opt.dvold,opt.dvolder,...
            f,dfdx,g,dgdx,opt.MMA.low,opt.MMA.upp,opt.MMA.a0,opt.MMA.a,opt.MMA.c,opt.MMA.d);        
        opt.iter = opt.iter+1;       
        opt.dvolder = opt.dvold; opt.dvold = opt.dv; opt.dv = dvnew;      
        opt.nv(opt.dof_dd,1) = opt.dv;
        
        % Continuation
        if (exist('cont_sw','var')==0) && (opt.deltaf < inputs.conv)
            cont_sw=1; cont_iter=0;
            disp('Continuation for Heaviside Proejction')
        elseif exist('cont_sw','var')
            cont_iter=cont_iter+1;
            if (mod(cont_iter,inputs.bt_ns)==1)
                opt.bt=opt.bt*inputs.bt_ic;
            end
        end    
    end
    fprintf('Elapsed time for Main Processing:%s\n',(datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')))

%% Post-processing
    fun_post(fem,opt,inputs); 
    fprintf('Elapsed time for Post-Processing:%s\n',(datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS')))