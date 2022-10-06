clear all; load Model_setup3;

ff=figure; 

subplot(1,2,1);
hold on; fs = 14;

h = bar(p.vacc*100);
set(gca,'XTickLabel',{'0-4yo','5-11yo','12-17yo','18-49yo','50-64yo','>65yo'}, 'fontsize', fs);

yl = ylim; yl(2) = 80; ylim(yl);
xlabel('Age group');
ylabel('Vaccination coverage, 2017/18 season');

xlim([0.2 6.8]);
w = h.BarWidth;
y = p.vacc(2)*100;
rectangle('Position',[2-w/2,0,w,y+11],'LineWidth',2,'LineStyle','--');
xtickangle(45);
title('Coverage');

subplot(1,2,2); hold on;
[num, txt, raw] = xlsread('Vaccine_data.xlsx','Vaccine_Data');
[r,c] = find(strcmp(raw,'VE'));
mat   = cell2mat(raw(r(2)+[1:6],c(2)+[0:2]));
mat2  = mat(:,[2,1,3])';

plot(1:6, mat2(2,:), '.', 'markersize', 24);
hilo = diff(mat2,1);
errorbar(1:6,mat2(2,:),hilo(1,:),hilo(2,:),'linestyle','None','linewidth',2);
xlim([0.6 6.5]);
set(gca,'XTickLabel',{'0-4yo','5-11yo','12-17yo','18-49yo','50-64yo','>65yo'}, 'fontsize', fs);
xtickangle(45);
xlabel('Age group');
ylabel('Vaccination efficacy (A+B)');
title('Efficacy');

set(ff, 'Position', [680   612   899   366]);