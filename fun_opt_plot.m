function [] = fun_opt_plot(fem,opt)

% Obejctive function history
    figure(1); clf(1); h=figure(1);
    set(h, 'Position', [5, 40, 400, 350]);
    plot(opt.fhis,'-k','linewidth',2);grid on;
    xlabel('Iteration');
    ylabel('Objective');
    
% Constraint function history
    figure(2);clf(2); h=figure(2);
    set(h, 'Position', [410, 40, 400, 350]);
    plot(opt.ghis+opt.volfrac,'-b','linewidth',2);hold on;
    xlabel('Iteration');
    ylabel('Volume');
    plot([1 opt.iter],[opt.volfrac opt.volfrac],':b','linewidth',2);grid on;ylim([0 1]);
    
% Plot density distribution
    figure(3); clf(3); h=figure(3);
    set(h, 'Position', [815, 40, 600, 350]);
    scatter3(fem.X(:,1),fem.X(:,2),fem.X(:,3),3,-opt.nrho)
    view([-0.5,-3,3.5]); axis equal; grid on;colormap('gray');drawnow

% Display computation time 
    if rem(opt.iter,10)==0; toc; end    
    
end
