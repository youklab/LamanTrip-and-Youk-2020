figure(21);
cur_pos = get(gcf,'position');
cur_w = 450;
cur_h = 320;
set(gcf,'units','pixels','position',[cur_pos(1),cur_pos(2),cur_w,cur_h]);
clf
hold on

font_name = 'Helvetica';
font_size_tick = 18;
font_size_label = 18;

set(get(gca,'XAxis'),'FontSize', font_size_tick, 'FontName', font_name);
set(get(gca,'YAxis'),'FontSize', font_size_tick, 'FontName', font_name);

set(gca, 'Layer', 'top');
xlabel('time (hours)', 'FontName', font_name, 'FontSize',font_size_label);
ylabel('colonies per mL', 'FontName', font_name, 'FontSize',font_size_label);
set(gca,'YScale','log');
set(gca, 'YMinorTick','On');
set(gca, 'TickLength', [0.025, 0.025]);
set(gca, 'LineWidth', 1.5);

mu=0.25;
Kmin=30000;
eps = 0.96;
Tmin = 37.9;
Tmax = 40.2;
Tscal = (Tmax-Tmin)/mu;

colors = { [1.0 0.05 0.0] [0.1 0.7 0.1] [0 0 0.85] };
draw_step = 1;

scale_cultures = 2;
tmax = 600;

T=40.8;
A0 = [11111 33333 100000];
curves = {[] [] []};
for (j = 1:3)
    curves{j} = [A0(j)];
    Mt = 0;
    At = A0(j);
    pa_t=0;
    pd = (T-Tmin)/Tscal;
    if (pd < 0)
        pd = 10^-10;
    elseif(pd>1)
        pd=1;
    end
    K = Kmin;
    Nt = A0;
    t=0;
    while (At > 0 && t < 50)
        pa_t = mu * Mt / (K + Mt);
        newborn = binornd(At,pa_t,1,1);
        At = At + newborn - binornd(At,pd,1,1);
        Mt = Mt + At;
        Nt = Nt + newborn;
        curves{j} = [curves{j} At];
        t=t+1;
    end
end

clf
hold on
step_1 = 3;
step_2 = 3;
step_3 = 3;
t1 = [0,4,8,12:10:(length(curves{1})-1)];
t2 = [0,4,8,12:10:(length(curves{2})-1)];
t3 = [0,4,8,12:10:(length(curves{2})-1)];

data1 = curves{1}(t1+1);
data2 = curves{2}(t2+1);
data3 = curves{3}(t3+1);

clf
cur_pos = get(gcf,'position');
cur_w = 450;
cur_h = 320;
set(gcf,'units','pixels','position',[cur_pos(1),cur_pos(2),cur_w,cur_h]);
hold on

sz = 80;
col_expon = [0.8 0.4 0.2];
col_power = [0.1 0.4 0.7];

xdat = 0:max(t1*1.05);

plot(t1, data1,  '.','MarkerSize',28,'MarkerFaceColor',[0.3 0.3 0.7]);%, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
plot(t2, data2,  '.','MarkerSize',28,'MarkerFaceColor',[0.7 0.3 0.3]);%, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);
plot(t3, data3,  '.','MarkerSize',28,'MarkerFaceColor',[0.3 0.7 0.3]);%, 'filled', 'MarkerFaceColor', [0.3 0.3 0.3]);

b = [ones(length(t2(1:3)),1), t2(1:3)'] \ log(data2(1:3)');
ydat = exp(b(1)) * exp(b(2) * xdat);
h4= plot(xdat, ydat, ':', 'Color', col_expon,  'LineWidth',2);

b = [ones(length(t1(1:3)),1), t1(1:3)'] \ log(data1(1:3)');
ydat = exp(b(1)) * exp(b(2) * xdat);
h4= plot(xdat, ydat, ':', 'Color', col_expon,  'LineWidth',2);

b = [ones(length(t3(1:3)),1), t3(1:3)'] \ log(data3(1:3)');
ydat = exp(b(1)) * exp(b(2) * xdat);
h4= plot(xdat, ydat, ':', 'Color', col_expon,  'LineWidth',2);

set(gca,'YLim',[10^1,10^5]);
set(gca,'XLim',[0,25]);
xlabel('time at 41.0^oC (hours)');
ylabel('survivors per mL');

font_name = 'Helvetica';
font_size = 18;
font_size_tick = 18;
font_size_label = 18;
set(get(gca,'XAxis'),'FontSize', font_size, 'FontName', font_name);
set(get(gca,'YAxis'),'FontSize', font_size, 'FontName', font_name);

set(gca,'YScale','log');
set(gca, 'YMinorTick','Off');
set(gca, 'TickLength', [0.025, 0.025]);
set(gca, 'LineWidth',1);

set(gca, 'YTick', [10^-6 10^-4 10^-2 10^0 10^2 10^4 10^6]);
set(gca, 'YTickLabel',['10^{-6}'; '10^{-4}'; '10^{-2}';  ' 10^{0}'; ' 10^{2}'; ' 10^{4}'; ' 10^{6}']);

set(gca, 'XTick', [0 25 50]);
set(gca, 'XTickLabel',[' 0 '; ' 25'; ' 50']);
