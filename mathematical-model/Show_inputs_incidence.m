clear all; 

load epid_data;
load popn_data;


mat = permute(squeeze(data_aligned(:,:,:,end)),[1,3,2]);                   % Dims: 1.Week 2.City 3.Age gp

tis = {'0-4 yrs','5-11 yrs','12-17 yrs','18-49 yrs','50-64 yrs','>65 yrs'};

% Get the per-capita rates
den = permute(repmat(agenums,[1,1,52]),[3,2,1]);
dat = mat./den*1e3;

ff=figure; fs = 14; lw = 1.5;
for ii = 1:6
    subplot(2,3,ii); hold on;
    
    inds = 1:20;
    plot(inds, dat(inds,:,ii),'linewidth',lw);
    set(gca,'ColorOrderIndex',1);
    
    inds = 20:30;
    plot(inds, dat(inds,:,ii),'--','linewidth',lw);
    
    if ismember(ii,[1,4])
        ylabel({'Weekly FluSurv-NET incidence','per 1000 age-specific population'});
    end
    if ismember(ii,4:6)
        xlabel('Calendar week');
    end
    
    set(gca,'fontsize',fs);
    title(tis{ii});
end
subplot(2,3,2); legend('Oakland','Contra Costa','Location','NorthWest');
set(ff,'Position',[680   402   899   576]);


% --- Get the total burden, helping to figure out excess hospitalisations
mat2 = squeeze(sum(dat(1:20,:,:),1));
pha = (1 - mat2(1,:)./mat2(2,:))