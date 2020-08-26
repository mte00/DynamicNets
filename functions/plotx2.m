function plotx2(t,y)
set(gcf,'DefaultAxesColorOrder',[0.8 0.1 0.1;1 0 0;1 0 0;0 0 1]);
cu=y(:,2);
cl=y(:,3);

h=t;
h=h;
hh=fill([h(1); h(1:end); flipud([h(1:end); h(end)])],[cu(1); cl(1:end); flipud([cu(1:end); cl(size(cl,1))])],'b');
set(hh,'edgecolor',[0 0.75 0.75]);
set(hh,'facecolor',[0 0.75 0.75]);

hold on;

 plot(h,y(:,1),'b','LineWidth',2);

% hold on;
% zz=zeros(size(y,1),1);
% plot(h,zz,'b-');
