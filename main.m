%% Main result for real data experiment (sparse matrix recovery)
clear all; close all;
rng(2019)

% Read data: the first three columns are responses, 
% the 4-5th column are confounders, the rest columns are SNPs
dat = csvread('GenDat.csv',1,2);

n = size(dat,1);
d1 = 2;
d2 = 250;  % set the number of considered SNPs, 
            % must less than size(dat,2) - d1 - 3

%% fit varying coefficient index model for each choice of response
Y = dat(:,1:3);
X = dat(:,4:5);
% random choose d2 SNPs
SelectSNP = 6 + sort(randsample(size(dat,2) - d1 - 3, d2));
Z = dat(:, SelectSNP);
c = 119/215;  % ratio of the first population

% define two score function and tunning parameters
s11 = @(x) c*normpdf(x,0,1).*x + (1-c)*normpdf(x,50,1).*(x-50);
s12 = @(x) c*normpdf(x,0,1) + (1-c)*normpdf(x,50,1);
s21 = @(x) tpdf(x,13)*14.*x./(13+x.^2) + tpdf(x-50,13)*14.*(x-50)./(13+(x-50).^2);
s22 = @(x) tpdf(x,13) + tpdf(x-50,13);
s1 = @(x) s11(x)./s12(x);
s2 = @(x) s21(x)./s22(x);

gam = 5*sqrt(log(d2)/n); % CLIMN tunning parameter
lam = sqrt(log(d1*d2)/n); % soft-thresholding tunning parameter
tau = (n/log(d1*d2))^(1/6); % hard-truncation tunning parameter

% Step 1: calculate truncated score variable for confounder x
SX = [s1(X(:,1)) s2(X(:,2))];
SX(abs(SX)>tau) = 0; % we only truncate x

% Step 2: estimate sparse precision matrix for covariance z
OmegaHat = CLIMECovZ(Z, n, d2, gam);
OmegaHat(abs(OmegaHat)<1e-7) = 0;  % computational error

save OmegaHat.mat OmegaHat % save precision matrix estimator for convenience
%load OmegaHat.mat

% Step 3: estimate signals with closed form
BetaHat = zeros(d1,d2,3);
for RR = 1:3
    betahat = SX'*(Y(:,RR).*Z)*OmegaHat/n;
    BetaHat(:,:,RR) = Soft_Thres(betahat,lam);
end

% Plot signal magnitude
for RR = 1:3
    for Conf = 1:2
        figure(1000 + RR*100 + Conf)
        plot(BetaHat(Conf,:,RR),'-r','LineWidth',2.5)
        xlabel('selected SNPs', 'Fontsize', 16)
        ylabel('signal magnitude','Fontsize', 16)
        filename = ['Figures/RealRes' num2str(RR) 'Con' num2str(Conf) '.png'];
        print('-dpng', filename)
    end
end


%% Run Nonparametric Bootstrapping for getting confidence interval
% set number of replications
NN = 100;
% run nonparametric bootstrapping
NonbootBetaHat = zeros(d1,d2,NN,3);
for i = 1:NN
    i
    SamID = randsample(n,n,true);
    bootY = Y(SamID,:);
    bootSX = SX(SamID,:);
    bootZ = Z(SamID,:);
    bootOmegaHat = CLIMECovZ(bootZ, n, d2, gam);
    bootOmegaHat(abs(bootOmegaHat)<1e-7) = 0;
    for RR = 1:3
        betahat = bootSX'*(bootY(:,RR).*bootZ)*bootOmegaHat/n;
        NonbootBetaHat(:,:,i,RR) = Soft_Thres(betahat,lam);
    end
end

% get 95% confidence interval
NonbootUpCI = zeros(d1,d2,3);
NonbootLoCI = zeros(d1,d2,3);
for k1 = 1:d1
    for k2 = 1:d2
        for RR = 1:3
            NonbootLoCI(k1,k2,RR) = quantile(NonbootBetaHat(k1,k2,:,RR),0.025);
            NonbootUpCI(k1,k2,RR) = quantile(NonbootBetaHat(k1,k2,:,RR),0.975);
        end
    end
end

% plot with confidence interval
for RR = 1:3
    for Conf = 1:2
        figure(2000 + RR*100 + Conf)
        plot(1:d2, BetaHat(Conf,:,RR),'-r',...
            1:d2, NonbootUpCI(Conf,:,RR),'-b',...
            1:d2, NonbootLoCI(Conf,:,RR),'-g','LineWidth',1.5)
        xlabel('selected SNPs', 'Fontsize', 16)
        ylabel('signal magnitude','Fontsize', 16)
        filename = ['Figures/NonCIRealRes' num2str(RR) 'Con' num2str(Conf) '.png'];
        print('-dpng', filename)
    end
end

%% Run parametric bootstrapping for getting confidence interval
% set number of replications
NNN = 10000;
% run parametric bootstrapping
ParbootBetaHat = zeros(d1,d2,NNN,3);
for i = 1:NNN
    bootX1 = [normrnd(0,1,119,1); normrnd(50,1,96,1)]; % generate confounder 1
    bootX2 = [trnd(13,60,1); 50 + trnd(13,59,1);trnd(13,48,1);50 + trnd(13,48,1)]; % generate confounder 2
    bootSX = [s1(bootX1) s2(bootX2)];
    bootSX(abs(bootSX)>tau) = 0;
    for RR = 1:3
        betahat = bootSX'*(Y(:,RR).*Z)*OmegaHat/n;
        ParbootBetaHat(:,:,i,RR) = Soft_Thres(betahat,lam);
    end
end

% get 95% confidence interval
ParbootUpCI = zeros(d1,d2,3);
ParbootLoCI = zeros(d1,d2,3);
for k1 = 1:d1
    for k2 = 1:d2
        for RR = 1:3
            ParbootLoCI(k1,k2,RR) = quantile(ParbootBetaHat(k1,k2,:,RR),0.025);
            ParbootUpCI(k1,k2,RR) = quantile(ParbootBetaHat(k1,k2,:,RR),0.975);
        end
    end
end

% plot with confidence interval
for RR = 1:3
    for Conf = 1:2
        figure(3000 + RR*100 + Conf)
        plot(1:d2, BetaHat(Conf,:,RR),'-r',...
            1:d2, ParbootUpCI(Conf,:,RR),'-b',...
            1:d2, ParbootLoCI(Conf,:,RR),'-g','LineWidth',1.5)
        xlabel('selected SNPs', 'Fontsize', 16)
        ylabel('signal magnitude','Fontsize', 16)
        filename = ['Figures/ParCIRealRes' num2str(RR) 'Con' num2str(Conf) '.png'];
        print('-dpng', filename)
    end
end











