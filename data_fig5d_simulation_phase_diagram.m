figure(20)
cur_pos = get(gcf,'position');
cur_w = 450;
cur_h = 350;
set(gcf,'units','pixels','position',[cur_pos(1),cur_pos(2),cur_w,cur_h]);
clf
hold on

set(gca, 'Layer', 'top');
set(gca, 'YTick', [10^1 10^3 10^5 10^7]);
set(gca, 'YTickLabel',['10^{1}'; '10^{3}'; '10^{5}'; '10^{7}']);
set(gca, 'XTick', [37 38 39 40 41]);
set(gca, 'XTickLabel',['37'; '38'; '39'; '40'; '41']);

set(get(gca,'XAxis'),'FontSize', font_size_tick, 'FontName', font_name);
set(get(gca,'YAxis'),'FontSize', font_size_tick, 'FontName', font_name);

xlabel('temperature (^oC)', 'FontName', font_name, 'FontSize',font_size_label);
ylabel('cells per mL', 'FontName', font_name, 'FontSize',font_size_label);
set(gca,'YScale','log');
set(gca, 'YMinorTick','Off');
set(gca, 'TickLength', [0.025, 0.025]);
set(gca, 'LineWidth', 2);

set(gca,'XLim',[37 41]);
set(gca,'YLim',[5 10^7]);

c_rg_shade = [0.55 0.9 0.55];
c_g_pt = [0 0 0.65];
c_ng_pt = [0.9 0.05 0.0];

sz = 70;
face_alpha = 1;

mu=0.25;
Kmin=30000;
eps = 0.96;
Tmin = 37.9;
Tmax = 40.2;
Tscal = (Tmax-Tmin)/mu;

data_g = [];
data_ng = [];

scale_cultures = 10;

temps = 37.0:0.5:41;
sim_data = [];
for (w=1:length(temps))
    A0 = round(10*unifrnd(0.95,1.05));
    n_cultures = 8;
    for (j = 1:7)
        q = temps(w)+normrnd(0,0.025);
        T=q;
        counter = 0;
        if (A0 < 10^6 || T>40.1)
            for (i=1:n_cultures)
                Mt = 0;
                At = A0;
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
                while (At > 0 && pa_t < pd)
                    pa_t = mu * Mt / (K + Mt);
                    newborn = binornd(At,pa_t,1,1);
                    At = At + newborn - binornd(At,pd,1,1);
                    Mt = Mt + At;
                end
                
                if (pa_t > pd)
                    counter = counter +1;
                end
                sim_data = [sim_data; temps(w) A0 (pa_t>pd)*1];
            end
            if (counter == 0)
                data_ng = [data_ng; q A0];
                scatter(q,A0, sz,'^', 'filled', 'MarkerFaceColor',c_ng_pt, 'MarkerFaceAlpha',face_alpha);
                drawnow
            elseif (counter == n_cultures)
                data_g = [data_g; q A0];
                scatter(q,A0, sz,'^', 'filled', 'MarkerFaceColor',c_g_pt, 'MarkerFaceAlpha',face_alpha);
                drawnow
            end
        end
        A0 = round(A0*scale_cultures*unifrnd(0.95,1.05));
    end
end
temps = 37.25:0.5:41;
sim_data = [];
for (w=1:length(temps))
    A0 = round(30*unifrnd(0.95,1.05));
    n_cultures = 8;
    for (j = 1:6)
        q = temps(w)+normrnd(0,0.025);
        T=q;
        counter = 0;
        if (A0 < 10^6 || T>40.1)
            for (i=1:n_cultures)
                Mt = 0;
                At = A0;
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
                while (At > 0 && pa_t < pd)
                    pa_t = mu * Mt / (K + Mt);
                    newborn = binornd(At,pa_t,1,1);
                    At = At + newborn - binornd(At,pd,1,1);
                    Mt = Mt + At;
                end
                
                if (pa_t > pd)
                    counter = counter +1;
                end
                sim_data = [sim_data; temps(w) A0 (pa_t>pd)*1];
            end
            if (counter == 0)
                data_ng = [data_ng; q A0];
                scatter(q,A0, sz,'^', 'filled', 'MarkerFaceColor',c_ng_pt, 'MarkerFaceAlpha',face_alpha);
                drawnow
            elseif (counter == n_cultures)
                data_g = [data_g; q A0];
                scatter(q,A0, sz,'^', 'filled', 'MarkerFaceColor',c_g_pt, 'MarkerFaceAlpha',face_alpha);
                drawnow
            end
        end
        A0 = round(A0*scale_cultures*unifrnd(0.95,1.05));
    end
end

c_g_shade = [0.7 0.9 1.0];
c_g_pt = [0 0 0.65];

c_ng_shade = [1.0 0.62 0.55];
c_ng_pt = [0.9 0.05 0.0];

c_rg_shade = [0.55 0.9 0.55];
c_rg_pt = [0.1 0.5 0.1];

boundary_g = [];
temp = 37.8;
numb = 0;
curdim = 1;
while (curdim(1)>0)
    cur_data = data_g(data_g(:,1)>temp,:);
    curdim=size(cur_data);
    if (curdim(1) > 0)
        [temp ind] = min(cur_data(:,1));
        boundary_g = [boundary_g; cur_data(ind,:)];
        boundary_g(boundary_g(:,2)>cur_data(ind,2),:)=[];
        temp = cur_data(ind,1);
    end
end
boundary_rg = [
    min(boundary_g(:,1)) 1;
    boundary_g];

boundary_stat = [
    35 10^8;
    35  max(boundary_g(:,2));
    max(boundary_g(:,1)) max(boundary_g(:,2))
    ];

boundary_g = [
    35 1;
    min(boundary_g(:,1)) 1;
    boundary_g;
    max(boundary_g(:,1)) 10^8;
    35 10^8
    ];

boundary_ng = [];
temp = 40.35;
numb = 0;
curdim = 1;
while (curdim(1)>0)
    cur_data = data_ng(data_ng(:,1)<temp,:);
    curdim=size(cur_data);
    if (curdim(1) > 0)
        [temp ind] = max(cur_data(:,1));
        boundary_ng = [boundary_ng; cur_data(ind,:)];
        boundary_ng(boundary_ng(:,2)<cur_data(ind,2),:)=[];
        temp = cur_data(ind,1);
    end
end
boundary_rg = [
    boundary_rg;
    boundary_ng;
    min(boundary_ng(:,1)) 1];

boundary_stat = [
    boundary_stat;
    max(boundary_ng(:,1)) max(boundary_ng(:,2));
    max(boundary_ng(:,1)) 10^8;
    ];

boundary_ng = [
    boundary_ng;
    min(boundary_ng(:,1)) 1;
    42.0 1;
    42.0 10^8;
    max(boundary_ng(:,1)) 10^8;
    ];

fill(gca, boundary_rg(:,1),boundary_rg(:,2), c_rg_shade, 'LineStyle', '-','EdgeColor',c_rg_pt);
fill(gca, boundary_g(:,1),boundary_g(:,2), c_g_shade, 'LineStyle', '--','EdgeColor',c_g_pt, 'LineWidth',1.5);

fill(gca, boundary_stat(:,1),boundary_stat(:,2), [0.9 0.9 0.9], 'LineStyle', '--','EdgeColor',[0.4 0.4 0.4], 'LineWidth',1.5);

fill(gca, boundary_ng(:,1),boundary_ng(:,2), c_ng_shade, 'LineStyle', '--','EdgeColor',c_ng_pt, 'LineWidth',1.5);
scatter(data_g(:,1),data_g(:,2), sz,'^', 'filled', 'MarkerFaceColor',c_g_pt, 'MarkerFaceAlpha',face_alpha);
scatter(data_ng(:,1),data_ng(:,2), sz,'^', 'filled', 'MarkerFaceColor',c_ng_pt, 'MarkerFaceAlpha',face_alpha);

T=[37.5:0.1:40.0 40.0:0.01:41];
pd = (T-Tmin)./Tscal;
pd(pd<0)=0;
pd(pd>1)=1;
phase_bound = (Kmin.*(1-eps/2)*eps.*pd.^2./(mu-eps.*pd));
phase_bound(pd==0) = 0;
phase_bound(pd==1) = 10^7;
plot(T, phase_bound, '-', 'Color',c_rg_shade/2,'LineWidth',5);
