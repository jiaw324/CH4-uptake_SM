
% We fit the Quadratic function  model
%
% $$ y = theta(1) x^ theta(2) (1-x)^theta(3) $$
%

%%\
% First clear some variables from possible previous runs.
clear data model options

%%
% import data from files 

data_measure=readtable('./CH4_oxi_20240206\D010-OC096.xlsx'); 


%%
% Next, create a data structure for the observations and control
% variables. Typically one could make a structure |data| that
% contains fields |xdata| and |ydata|.

data.xdata = data_measure.WFPS;   % wpfs 
data.ydata = data_measure.CH4FLUX; % CH4 flux (CH4/mgCm-2h-1) 

%%
% Here is a plot of the data set.
figure(1); clf
% figure
plot(data.xdata,data.ydata,'s','MarkerSize',12);
xlim([0 1]); xlabel('WFPS'); ylabel('CH_4 flux [\mumole-CH_4/(g h)]');  %for CH4_oxi_20240206
% xlim([0 1]); xlabel('WFPS'); ylabel('CH_4 flux [\mug CH_4-C/(m^2 h)]'); %for CH4_oxi_20240206
xticks([0 0.2 0.4 0.6 0.8 1])
xticklabels({'0','20%','40%','60%','80%','100%'})
set(gca,'FontSize',18)
%
% For the MCMC run we need the sum of squares function. For the
% plots we shall also need a function that returns the model.
% Both the model and the sum of squares functions are
% easy to write as one line anonymous functions using the @
% construct. 

modelfun = @(x,theta) theta(1)*x.^theta(2).*(1-x).^theta(3); %objective function 
ssfun    = @(theta,data) sum((data.ydata-modelfun(data.xdata,theta)).^2);


%%
% In this case the initial values for the parameters are easy to guess
% by looking at the plotted data. As we alredy have the sum-of-squares
% function, we might as well try to minimize it using |fminsearch|.
% given initial theta
theta0 = [50,0.6,1.2];  % the values are important for sampling cannot be fitted by squaratic well
 [tmin,ssmin] = fminsearch(ssfun,theta0,[],data);  % find local minimizer of unconstrained nonlinear multivariable function using Nelder-Mead algorithm, but may produce initials resulting in NaN  
% tmin = theta0;
 % [tmin,ssmin] = fminunc(ssfun,theta0,[],data);  % find local minimizer of unconstrained nonlinear multivariable function using 'quasi-newton' (default) or specifiying other algorithms,   
%[] = fmincon(...); %find local minimizer of constrained nonlinear multivariable function 
n = length(data.xdata);
p = 2;
mse = ssmin/(n-p); % estimate for the error variance
%%
% The Jacobian matrix of the model function is easy to calculate so we use
% it to produce estimate of the covariance of theta. This can be
% used as the initial proposal covariance for the MCMC samples by
% option |options.qcov| below.

J = [tmin(1)*(1-data.xdata).^tmin(3)*tmin(2).*data.xdata.^(tmin(2)-1),...
    tmin(1)*data.xdata.^tmin(2)*tmin(3).*(1-data.xdata).^(tmin(3)-1),...
    data.xdata.^tmin(2).*(1-data.xdata).^tmin(3)];
% J = [theta0(1)*(1-data.xdata).^theta0(3)*theta0(2).*data.xdata.^(theta0(2)-1),...
%     theta0(1)*data.xdata.^theta0(2)*theta0(3).*(1-data.xdata).^(theta0(3)-1),...
%     data.xdata.^theta0(2).*(1-data.xdata).^theta0(3)];

J(J==Inf) = 0; % make sure tcov is not Inf 

% [jm,jn] = size(J);
% for i=1:jm
%     for j=1:jn
%         if(~isreal(J(i,j)))
%             J(i,j) = real(0);
%         end
%     end
% end

tcov = inv(J'*J)*mse;


%%
% We have to define three structures for inputs of the |mcmcrun|
% function: parameter, model, and options.  Parameter structure has a
% special form and it is constructed as Matlab cell array with curly
% brackets. At least the structure has, for each parameter, the name
% of the parameter and the initial value of it. Third optional
% parameter given below is the minimal accepted value. With it we set
% a positivity constraits for both of the parameters.

% tmin = theta0;  % turn on when the estimated initials doesn't work

params = {
    {'theta1', tmin(1), 0,10000}   
    {'theta2', tmin(2), 0,1000}
    {'theta3', tmin(3), 0,1000}
    };

%%
% In general, each parameter line can have up to 7 elements: 'name',
% initial_value, min_value, max_value, prior_mu, prior_sigma, and
% targetflag

%%
% The |model| structure holds information about the model. Minimally
% we need to set |ssfun| for the sum of squares function and the
% initial estimate of the error variance |sigma2|.

model.ssfun  = ssfun;
model.sigma2 = mse; % (initial) error variance from residuals of the lsq fit

%%
% If we want to sample the error variance sigma2 as an extra model
% parameter, we need to set a prior distribution for it. A convenient
% choice is the conjugate inverse chi-squared distribution, which
% allows Gibbs sampling step for sigma2 after every
% Metropolis-Hastings update for the other parameters. This is
% acchieved by |options.updatesigma=1|, below. The default prior is
% uninformative, but we can set the prior parameters with the
% following options. Option |model.N| for the number of observatios is
% needed for 'updatesigma', if it is not given, the code tries to
% guess |N| from the data.

model.N = length(data.ydata);  % total number of observations
model.S20 = model.sigma2;      % prior mean for sigma2
model.N0  = 4;                 % prior accuracy for sigma2

%%
% The |options| structure has settings for the MCMC run. We need at
% least the number of simulations in |nsimu|. Here we also set the
% option |updatesigma| to allow automatic sampling and estimation of the
% error variance. The option |qcov| sets the initial covariance for
% the Gaussian proposal density of the MCMC sampler.

options.method = 'dram';
options.nsimu = 5000;
options.updatesigma = 1;
options.qcov = tcov; % covariance from the initial fit
% options.updatesigma = 0; % means signa2 is fixed


%%
% The actual MCMC simulation run is done using the function
% |mcmcrun|.


rng('default')  %set a seed to fix initial sampler 

[res,chain,s2chain] = mcmcrun(model,data,params,options);


i1 = fix(0.5*options.nsimu);
i2 = fix(options.nsimu);

theta = mean(chain(i1:i2,:));

% % constrain theta to make sure a+b+c>0
% n0 = 0;
% chain_sum = [0,0,0];
% for i = 1:options.nsimu
%    if sum(chain(i,:))>0
%       % i
%         n0 = n0+1;
%         chain(i,:);
%         chain_sum = chain_sum + chain(i,:);
%    end
% end
% theta = chain_sum/n0;

% redirect terminal I/O to file
diary ./CH4_oxi_20240206\D010-OC096_mcmc
% diary disp(fname3)
RR = (corrcoef(data.ydata,modelfun(data.xdata,theta))).^2
% RR = (corrcoef(data.ydata,modelfun(data.xdata,mean(chain)))).^2
%RR = (corrcoef(data.ydata(2:n),modelfun(data.xdata(2:n),mean(chain)))).^2
RMSE = sqrt(mean((modelfun(data.xdata,theta)-data.ydata).^2))
NRMSE = sqrt(mean((modelfun(data.xdata,theta)-data.ydata).^2))/mean(data.ydata)    %normalize RMSE
cc = theta(1)
alpha = theta(2)
beta = theta(3)
wfps_opt = alpha/(alpha+beta)
Fmax = cc*(alpha/(alpha+beta))^alpha*(beta/(alpha+beta))^beta
diary off
%%
% 
% During the run, a status window is showing the estimated time to
% the end of the simulation. The simulation can be ended by Cancel
% button and the chain generated so far is returned.

%%
% After the run the we have a structure |res| that contains some
% information about the run and a matrix outputs |chain| and
% |s2chain| that contain the actual MCMC chains for the parameters
% and for the observation error variance.

% %%
% % The |chain| variable is |nsimu| × |npar| matrix and it can be
% % plotted and manipulated with standard Matlab functions. MCMC toolbox
% % function |mcmcplot| can be used to make some useful chain plots and
% % also to plot 1 and 2 dimensional marginal kernel density estimates of
% % the posterior distributions.
%
figure(2); clf
mcmcplot(chain,[],res,'chainpanel');


 cd ./CH4_oxi_20240206
 saveas(gcf,'D001-OC103_mcmc_chain','jpg')
 cd ../


%%
% The |'pairs'| options makes pairwise scatterplots of the CH4_oxi_20240206s of
% the |chain|.
% 
% figure(3); clf
% mcmcplot(chain,[],res,'pairs');
% % 
% %%
% % If we take square root of the |s2chain| we get the chain for error
% % standard deviation. Here we use |'hist'| option for the histogram of
% % the chain.
% 
% figure(4); clf
% mcmcplot(sqrt(s2chain),[],[],'hist')
% title('Error std posterior')
% 
% % add prior distribution to the plot, if it was informative
% if res.N0>0
%   xl = xlim; xx = linspace(xl(1),xl(2));
%   hold on
%   plot(xx,invchi1pf(xx,res.N0,sqrt(res.S20)))  
%   hold off
%   legend('posterior','prior')
% end

%%
% A point estimate of the model parameters can be calculated as the
% mean of the |chain|. Here we plot the fitted model using the
% posterior means of the parameters.

x = linspace(0,1)';
figure(1);
hold on
% plot(x,modelfun(x,theta),'-k','LineWidth',3)
plot(x,modelfun(x,mean(chain(i1:i2,:))),'-k','LineWidth',3)
plot(x,0*x,':')
hold off
legend({'data','model'},'Location','best')
 
 
%%
% Instead of just a point estimate of the fit, we should also study
% the predictive posterior distribution of the model. The |mcmcpred|
% and |mcmcpredplot| functions can be used for this purpose. By them
% we can calculate the model fit for a randomly selected subset of the
% chain and calculate the predictive envelope of the model. The grey
% areas in the plot correspond to 50%, 90%, 95%, and 99% posterior
% regions.

figure(5); clf 

% out = mcmcpred(res,chain(i1:i2,:),s2chain(i1:i2,1),x,modelfun);
% out = mcmcpred(res,chain(i1:i2,:),[],x,modelfun);
out = mcmcpred(res,chain(i1:i2,:),s2chain,x,modelfun);
mcmcpredplot(out);

% If s2chain has been given to mcmcpred, then the plot shows 95%
% probability limits for new observations and for model parameter
% uncertainty. If s2chain is not used then the plot contains 50%,
% 90%, 95%, and 99% predictive probability limits due parameter uncertainty.

hold on

plot(data.xdata,data.ydata,'s','MarkerSize',12,'MarkerEdgeColor','b');%,'MarkerFaceColor','{}');
xlim([0 1]); xlabel('WFPS'); ylabel('CH_4 flux [\mumole-CH_4/(g h)]');  %for CH4_oxi_20240206
% xlim([0 1]); xlabel('WFPS'); ylabel('CH_4 flux [\mug CH_4-C/(m^2 h)]'); %for CH4_oxi_20240206
xticks([0 0.2 0.4 0.6 0.8 1])
xticklabels({'0','20%','40%','60%','80%','100%'})
set(gca,'FontSize',18)


 cd ./CH4_oxi_20240206
 saveas(gcf,'D001-OC103_mcmc','jpg')
 cd ../


