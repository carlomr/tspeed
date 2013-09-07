function [a,b] = load_and_plot(fname,numrecvs, start)

if nargin==1
    start = 0;
end

for i=1:numrecvs
    fname = sprintf('%s_rcv_%i.out',fname, i-1);

    a{i}=load(fname0);
%b=load(fname1);
a(1:start-1,:) = zeros(start-1,2);
b(1:start-1,:) = zeros(start-1,2);
figure
plot(linspace(0,1,length(a(1:end,1))), a(1:end,1),'k--');
hold on;
plot(linspace(0,1,length(b(1:end,1))), b(1:end,1),'r-');
set(gca, 'Ytick', []);
title('u_x');;
legend('r_1','r_2', 'location', 'best')
set(findall(gcf,'type','text'),'fontSize',14);
set(gca, 'fontsize',14)
%legend boxoff

figure
hold on;
plot(linspace(0,1,length(a(1:end,1))), a(1:end,2),'k--');
plot(linspace(0,1,length(b(1:end,1))), b(1:end,2),'r');
xlabel('t');
set(gca, 'Ytick', []);
title('u_y');;
legend('r_1','r_2', 'location', 'best')
set(findall(gcf,'type','text'),'fontSize',14);
set(gca, 'fontsize',14)
%legend boxoff

