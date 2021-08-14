mu=0.25;
Kmin=30000;
eps = 0.96;
Tmin = 37.9;
Tmax = 40.2;
Tscal = (Tmax-Tmin)/mu;

tmax = 600;

figure(22);
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
xlabel('# of cells / mL', 'FontName', font_name, 'FontSize',font_size_label);
ylabel('population half-life (hours)', 'FontName', font_name, 'FontSize',font_size_label);
set(gca,'YScale','log');
set(gca, 'YMinorTick','On');
set(gca, 'TickLength', [0.025, 0.025]);
set(gca, 'LineWidth', 1.5);

set(gca,'XScale','log');
set(gca,'YScale','log');

Tcrit = Tmax;
nr_reps = 3;
for (T=[Tcrit,39.1,40,40.8,42,43])
    nr_densities = 5;
    if (T==39.1)
        A0 = [7500 9500];
        for (i = 1:3)
            A0 = [round(A0(1)*0.5) A0];
        end
        cur_col= [0.85 0.75 0.0  ];
    end
    if (T==40)
        A0 = [150000];
        for (i = 1:5)
            A0 = [round(A0(1)*0.33) A0];
        end
        cur_col= [0.8  0.6  0.1  ];
    end
    if (T==Tcrit)
        A0 = [660000];
        for (i = 1:2)
            A0 = [A0 round(A0(end)*0.33)];
        end
        cur_col= [0.75 0.45 0.2  ];
    end
    if (T==40.8)
        A0 = [5000000];
        for (i = 1:nr_densities)
            A0 = [A0 round(A0(end)*0.2)];
        end
        cur_col= [0.7  0.3  0.3  ];
    end
    if (T==42)
        A0 = [5000000];
        for (i = 1:nr_densities)
            A0 = [A0 round(A0(end)*0.3)];
        end
        cur_col= [0.65 0.15 0.4  ];
    end
    
    if (T==43)
        A0 = [5000000];
        for (i = 1:3)
            A0 = [A0 round(A0(end)*0.4)];
        end
        cur_col= [0.6  0.0  0.5  ];
    end
    curves = {[] [] []};
    
    death_rate = zeros(nr_densities,nr_reps);
    density = zeros(nr_densities,1);
    
    for (j = 1:length(A0))
        for (k=1:nr_reps)
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
            t = 0;
            while (At > 0 && t < 10)
                pa_t = mu * Mt / (K + Mt);
                newborn = binornd(At,pa_t,1,1);
                At = At + newborn - binornd(At,pd,1,1);
                Mt = Mt + At;
                Nt = Nt + newborn;
                curves{j} = [curves{j} At];
                t=t+1;
            end
            P=polyfit(0:t,log(curves{j}),1);
            death_rate(j,k) = -log(2)/P(1);
            
            density(j) = A0(j);
        end
    end
    if (T==39.1)
        density(6) = density(5)*(9500/7500);
        death_rate(6,1) = 500;
        death_rate(6,2) = 500;
        death_rate(6,3) = 500;
    end
    if (T==40)
        density(7) = density(end)*3;
        death_rate(7,1) = 500;
        death_rate(7,2) = 500;
        death_rate(7,3) = 500;
    end
    errorbar(density,mean(death_rate,2),std(death_rate')/sqrt(3),std(death_rate')/sqrt(3),    'o--','LineWidth',1.5,'MarkerSize',6,'Color',cur_col);
end

set(gca, 'XTick', [10^3 10^5 10^7]);
set(gca, 'YTick', [1 10 100 500]);
set(gca, 'XTickLabel',['10^3'; '10^5'; '10^7']);
set(gca, 'YTickLabel',['10^0'; '10^1'; '10^2';' 500']);

set(gca,'XLim',[300 10^7]);
set(gca,'YLim',[10^0 5*10^2]);
set(gca, 'YMinorTick','Off');
set(gca, 'XMinorTick','Off');
