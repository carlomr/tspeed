function  load_and_plot(fn,startrecvs,numrecvs, start)

if nargin==3
    start = 0;
end

ofsx = 0;
ofsy = 0;
for i=startrecvs:numrecvs+startrecvs-1
    fname = sprintf('%s_rcv_%i.out',fn, i-1);

    a{i}=load(fname);
    a{i}(1:start-1,:) = zeros(start-1,2);
    myofsx = max(a{i}(:,1));
    myofsy = max(a{i}(:,2));
    if myofsy>ofsy
	ofsy = myofsy;
    end
    if myofsx>ofsx
	ofsx = myofsx;
    end

end
t = linspace(0,1,length(a{end}(:,1)));

figure;
subplot(2,1,1);
hold on;

title('u_x');
xlabel('t');
for i=startrecvs:numrecvs+startrecvs-1
    j = i-startrecvs+1;
    plot(t(:),a{i}(:,1)*3+j*ofsx,'r','LineWidth',1.5);
end
hold off
axis tight
set(gca,'YTick',linspace(1*ofsx,numrecvs*ofsx,floor(numrecvs/2)));
set(gca,'YTickLabel', num2cell(1:2:numrecvs));
set(findall(gcf,'type','text'),'fontSize',14);
set(gca, 'fontsize',14)
subplot(2,1,2);
hold on;
title('u_y');
xlabel('t');
for i=startrecvs:numrecvs+startrecvs-1
    j = i-startrecvs+1;
    plot(t(:),a{i}(:,2)*3+j*ofsy,'k','LineWidth',1.5);
end
set(gca,'YTick',linspace(1*ofsy,numrecvs*ofsy,floor(numrecvs/2)));
set(gca,'YTickLabel', num2cell(1:2:numrecvs));
set(findall(gcf,'type','text'),'fontSize',14);
set(gca, 'fontsize',14)
axis tight;
