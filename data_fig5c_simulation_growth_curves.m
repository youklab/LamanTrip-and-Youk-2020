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
ylabel('cells / mL', 'FontName', font_name, 'FontSize',font_size_label);
set(gca,'YScale','log');
set(gca, 'YMinorTick','On');
set(gca, 'TickLength', [0.025, 0.025]);
set(gca, 'LineWidth', 1.5);

set(gca, 'YTick', [10^2 10^3 10^4 10^5 10^6 10^7]);
set(gca, 'YTickLabel',['10^{2}'; '10^{3}'; '10^{4}'; '10^{5}'; '10^{6}'; '10^{7}']);
set(gca,'YLim',[10^2 10^7]);


mu=0.25;
Kmin=30000;
eps = 0.96;
Tmin = 37.9;
Tmax = 40.2;
Tscal = (Tmax-Tmin)/mu;

colors = { [1.0 0.05 0.0] [0.1 0.7 0.1] [0 0 0.85] };
draw_step = 100;

scale_cultures = 5;
tmax = 600;

T=39.1;
A0 = 385;
n_cultures = 8;
for (j = 1:3)
    for (i=1:n_cultures)
        Mt = 0;
        At = A0;
        if (j == 1)
            At = round(A0 * (1+2*(i/n_cultures)));
        elseif (j==3)
            At = round(A0 * (1+2*(i/n_cultures)));
        end
        pa_t=0;
        pd = (T-Tmin)/Tscal;
        if (pd < 0)
            pd = 10^-10;
        elseif(pd>1)
            pd=1;
        end
        K = Kmin;
        Nt = A0;
        Nt_old = Nt;
        t=0;
        told=t;
        while (At > 0 && Nt < 5*10^6 || At==0 && t<tmax)
            pa_t = mu * Mt / (K + Mt);
            newborn = binornd(At,pa_t,1,1);
            At = At + newborn - binornd(At,pd,1,1);
            Mt = Mt + At;
            Nt = Nt + newborn;
            t = t + 1;
            if (mod(t,draw_step)==0 || Nt_old/Nt < 0.2)
                plot([told t], [Nt_old Nt] ,'o--', 'Color', colors{j}, 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 0], 'MarkerSize',4);
                told = t;
                Nt_old = Nt;
                drawnow
            end
        end
        plot([told t], [Nt_old Nt] ,'o--', 'Color', colors{j}, 'LineWidth', 1.5, 'MarkerEdgeColor', [0 0 0], 'MarkerSize',4);
        told = t;
        Nt_old = Nt;
        drawnow
    end
    A0 = A0*scale_cultures;
end