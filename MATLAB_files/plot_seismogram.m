function plot_seismogram(r, dt)

nr = length(r(1,1,:));
ns = length(r(1,:,1));
t = 2*dt:dt:dt*(ns+1);
ofs = max(max(max(r)))/2;

figure;
subplot(2,1,1);
hold on;
for i=1:nr
    plot(t(:),r(1,:,i)+i*ofs,'r','LineWidth',2);
end
hold off
axis tight
set(gca,'YTick',linspace(1*ofs,nr*ofs,nr));
set(gca,'YTickLabel', num2cell(1:nr));
subplot(2,1,2);
hold on;
for i=1:nr
    plot(t(:),r(2,:,i)+i*ofs,'k','LineWidth',2);
end
axis tight
set(gca,'YTick',linspace(1*ofs,nr*ofs,nr));
set(gca,'YTickLabel', num2cell(1:nr));
