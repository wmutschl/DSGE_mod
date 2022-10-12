@#define orderApp = 2

var a va siga junk;
varobs a va siga; % reporting variables
varexo eps_va eps_siga;
parameters RHOA RHOSIGA Ass SIGAss;

RHOA = 0.98;
RHOSIGA = 0.99;
Ass = 1;
SIGAss = 0.02;

model;
0 = -log(exp(a)/Ass) + exp(siga)*log(exp(va));
0 = -log(exp(va)) + RHOA*exp(siga(-1))/exp(siga)*log(exp(va(-1))) + eps_va;
0 = -log(exp(siga)/SIGAss) + RHOSIGA*log(exp(siga(-1))/SIGAss) + eps_siga;
junk=0.9*junk(+1); % junk equation because 2nd and 3rd order approximation are not implemented for purely backward models
end;

steady_state_model;
a = log(Ass);
va = log(1);
siga = log(SIGAss);
end;

shocks;
var eps_va = 1^2;
var eps_siga = 0.0265^2;
end;
steady;
stoch_simul(order=@{orderApp},irf=0,periods=0,nocorr,nodecomposition,ar=0);

%% MANUAL SIMULATION
verbatim; 
    options_.periods = 2000;
    % create joint index for state variables x and varobs variables y (we will focus only on these variables below)
    indx = M_.nstatic + (1:M_.nspred); % state variables in DR order
    indy = [];
    for j=1:length(options_.varobs) %use for loop to keep same ordering as declared in varobs
        indy = [indy find(ismember(M_.endo_names(oo_.dr.order_var),options_.varobs{j}))]; % reporting variables in DR order
    end
    indxy = [indx(:);indy(:)]; % joint index for state and observable variables
    xy_ss = oo_.dr.ys(oo_.dr.order_var); xy_ss = xy_ss(indxy); % steady-state of x and y in DR order

    % select rows in perturbation matrices for state and observable variables only
    gx = oo_.dr.ghx(indxy,:); gu = oo_.dr.ghu(indxy,:);
    if options_.order > 1
        gxx = oo_.dr.ghxx(indxy,:); gxu = oo_.dr.ghxu(indxy,:); guu = oo_.dr.ghuu(indxy,:); gss = oo_.dr.ghs2(indxy,:);
    end

    % get shock series
    randn('seed',1)
    exo = transpose(chol(M_.Sigma_e)*randn(M_.exo_nbr,options_.periods)); % draw from standard normal and multiply with standard deviations

    % initialize
    xhat = zeros(length(indx),1);
    yhat = zeros(length(indy),options_.periods);

    % do the simulations
    for t = 1:options_.periods
        u = exo(t,:)';
        xyhat = gx*xhat + gu*u;
        if options_.order > 1
            xyhat = xyhat + 1/2*gxx*kron(xhat,xhat) + gxu*kron(xhat,u) + 1/2*guu*kron(u,u) + 1/2*gss;
        end         
        xhat = xyhat(1:length(indx)); % update states (in deviation from steady-state)
        yhat(:,t) = xyhat((length(indx)+1):end); % update observables
    end    
    
    TBL = zeros(length(options_.varobs),4);
    for j = 1:length(indy)
        TBL(j,1) = mean(yhat(j,:)');
        TBL(j,2) = std(yhat(j,:)');
        TBL(j,3) = skewness(yhat(j,:)');
        TBL(j,4) = kurtosis(yhat(j,:)');
    end
    disp(array2table(TBL,'RowNames',options_.varobs,'VariableNames',{'MEAN','STD','SKEWNESS','KURTOSIS'}));
    
end; //verbatim end