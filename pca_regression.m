function [coeff_loadings,score,latent,explained,rsq_PCR] = pca_regression(Xdata,Ydata,nPC)
%pca_regression Function runs principal component analysis and conducts PC
%regression with nPC components, including plots
%   Arguments:
%   Xdata = observations x variables input data matrix
%   Ydata = response variable array
%   nPC = number of PCs to regress
%
%   Outputs:
%   coeff_loadings = Principal component coefficients are recipe for counting any given PC
%       each column of coeff contains coefficients for one principal component.
%       columns are in order of descending component variance, latent.
%   score = how each individual observation is composed of the PCs. matrix of PCs x observations.
%   latent = principal component variances, the eigenvalues of the covariance matrix, returned as a column vector
%       an eigenvalue is the total amount of variance in the variables in the dataset explained by the common factor.
%       PCs are typically retained for further analysis when their eigenvalues (of the transformation matrix) are larger than 1.
%   explained = the contribution of each PC to variability in data
%   rsq_PCR = r2 for PCR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Normalize data
% normalize manually % (data-nanmean(data))./nanstd(data)
nXdata=normalize(Xdata);

% % remove rows with any NaNs
% Xdata(find(sum(isnan(nXdata)')>0),:)=[];
% Ydata(find(sum(isnan(nXdata)')>0),:)=[];
% nXdata(find(sum(isnan(nXdata)')>0),:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conduct PCA
[coeff_loadings,score,latent,~,explained]=pca(nXdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot PCA
% biplot of 2 PCs
close all
figure('position',[100 50 1200 900],'paperpositionmode','auto');
subplot(2,2,1)
for n=1:size(Xdata,2)
    labels{n,:}=string(n);
end
biplot(coeff_loadings(:,1:2),'scores',score(:,1:2),'varlabels',string(labels));
box on
title('BIPLOT')

% scree plot
subplot(2,2,2)
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
title('SCREE PLOT')
grid on

% profile plot
subplot(2,2,3)
plot(corr(score(:,1:size(Xdata,2)),nXdata(:,1:size(Xdata,2)),'rows','complete')','o-')
xlabel('Variable')
ylabel('Correlation')
ylim([-1 1])
title('PROFILE PLOT')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Principal component regression
% https://www.mathworks.com/help/stats/partial-least-squares-regression-and-principal-components-regression.html?prodcode=ST&language=en
if nargin>1
    if nargin<3
        % set nPC as PCs with eigenvalues >=1
        nPC=sum(latent>=1);
    end
    % PCR is a linear regression of the response variable on n principal components
    betaPCR = regress(Ydata-mean(Ydata), score(:,1:nPC));
    
    % to make the PCR results easier to interpret in terms of the original data, transform to regression coefficients for the original, uncentered variables
    betaPCR = coeff_loadings(:,1:nPC)*betaPCR;
    betaPCR = [nanmean(Ydata) - nanmean(nXdata)*betaPCR; betaPCR];
    yfitPCR = [ones(size(nXdata,1),1) nXdata]*betaPCR;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate principal component regression r-square
    % total sum of squares
    TSS = sum((Ydata-nanmean(Ydata)).^2);
    % residual sum of squares
    RSS_PCR = nansum((Ydata-yfitPCR).^2);
    rsq_PCR = 1 - RSS_PCR/TSS
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot principal component regression
    % plot fitted vs. observed response for the PCR fit
    subplot(2,2,4)
    plot(Ydata,yfitPCR,'bo');
    xlabel('Observed Response');
    ylabel('Fitted Response');
    title([num2str(nPC),' PCR; rsq = ',num2str(rsq_PCR)]);
    grid on
    set(findall(gcf,'-property','FontName'),'FontName','Gill Sans MT','FontSize',13)
end
end