% 'Equation for Describing Solute Transport in Field soils with Preferential
% Flow Paths'
% Kim et al., 2005

% The authors propose a simple equation that can predict the breakthrough
% of solutes without excessive data requirements under steady flow
% conditions.

% GENERALIZED PREFERENTIAL FLOW MODEL (GPFM): the soil is conceptually
% divided into a distribution zone and a conveyance zone. The distribution
% zone acts as a linear reservoir resulting in an exponential loss of
% solutes through preferential flow paths from this zone. In the conveyance
% zone, the transport of solutes is described with the
% convective-dispersive equation.

% INPUT DATA: apparent water content of the distribution zone, and solute
% velocities and dispersivities in the conveyance zone.

% The thickness of the distribution zone depends on land use and tillage
% practices.

close all
clear all
clc

BTC = xlsread('Data_Spectro.xlsx',1); % 1) t, 2) C, 3) V per fraction (m), 4) Flow rate (ml/s)

t = BTC(:,1); % (h)
C_tracer = BTC(:,2);
x = 50; % (cm) Height of the chamber
q = 0.3; % (cm/s) Steady state flow rate

    d = 26; % (cm) depth of distribution zone
    rho_b = 1.5; % (g/cm3) bulk density of the soil. 1.7 g/cm3 for coarse sand; 1.65 g/cm3 for fine sand.
    rho = 2.65; % (g/cm3) mean density of quarz
    n = 1 - rho_b/rho; % Porosity of the media.
    theta_s = n; % (cm3/cm3) We take the saturated water content as being equal to the porosity of the media.
    Kd = 0.19;  % Adsorption-desorption Distribution Coefficient; parameter for
        % understanding the mobility of a compound in the environment
        % (partitioning), and its distribution between water, sludge, soil and
        % sediments. The higher Kd the compound is less likely to move in soil.
w = d*(rho_b*Kd + theta_s); % (m) apparent water content in the distribution zone;
                            % For nonadsorbed chemicals, w is simply
                            % defined as the volume water in the
                            % distribution zone per unit area. For adsorbed
                            % chemicals, w =d(rho*Kd + theta_s).

%v = q/(beta*(rho*Kd + theta));
v = 0.625%v = 0.625; % (cm/s)
theta = theta_s; % (cm3/cm3) is used as the equilibrium volumetric moisture content in the finger (Selker et al., 1992)
beta = q/(v*theta);
eta = q/w;
C0 = 0.85; 
Dt = 240; 

    % a) Newton-Raphson's algorithm 
syms P

f = @(P) 0.5.*C0.*erfc((x-v.*t(find(t>Dt):end))./(sqrt(4.*P.*t(find(t>Dt):end)))) - exp(v.*x.*(1- sqrt(1-4.*P.*eta./v.^2))./(2.*P)-eta.*t(find(t>Dt):end)).*erfc((x-v.*t(find(t>Dt):end).*sqrt(1-4.*P*eta/v.^2))./(sqrt(4.*P.*t(find(t>Dt):end))))   -    0.5.*C0.*erfc((x-v.*(t(find(t>Dt):end)-Dt))./(sqrt(4.*P.*(t(find(t>Dt):end)-Dt)))) + exp(v.*x.*(1- sqrt(1-4.*P.*eta./v.^2))./(2.*P)-eta.*(t(find(t>Dt):end)-Dt)).*erfc((x-v.*(t(find(t>Dt):end)-Dt).*sqrt(1-4.*P.*eta./v.^2))./(sqrt(4.*P.*(t(find(t>Dt):end)-Dt))))  - C_tracer(find(t>Dt):end);
    
f_prime = diff(f,P,1); 

P_new = -(0.1^2-1)*v^2/(4*eta)  ; 
f_prime_eval = eval(subs(f_prime,P,P_new));

%nanmean(eval(subs(f,P,P_new))./f_prime_eval)
tolerance = 0.1;
estimate = 1;

for i=1:20
    if estimate> tolerance
        if f_prime_eval ~= 0
            try
                check = sqrt(1-4.*P_new.*eta./v.^2);
                f_prime_eval = eval(subs(f_prime,P,P_new));    
                P_old = P_new;
                P_new = P_old - nanmean(eval(subs(f,P,P_old))./f_prime_eval);
                estimate =  abs(P_new - P_old);
                j = i+1;
            catch
                P_new = P_old;
                estimate = abs(P_new - P_old);
                j=j+1;
            end
        end
    end
end

P = P_new;
    C_C0_fit = 0.5.*erfc((x-v.*t)./(sqrt(4.*P.*t))) - exp(v.*x.*(1- sqrt(1-4.*P.*eta./v.^2))./(2.*P)-eta.*t).*erfc((x-v.*t.*sqrt(1-4.*P.*eta./v.^2))./(sqrt(4.*P.*t)))./P;
    C_C0_delay = zeros(length(t),1);

  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
            C_C0_delay(i) =0.5*erfc((x-v*(t(i)-Dt))/(sqrt(4*P*(t(i)-Dt)))) - exp(v*x*(1- sqrt(1-4*P*eta/v^2))/(2*P)-eta*(t(i)-Dt))*erfc((x-v*(t(i)-Dt)*sqrt(1-4*P*eta/v^2))/(sqrt(4*P*(t(i)-Dt))))/P;
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;

% Goodness of fit:
t_shift = 0; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

sim_cal_GPFM = F_shifted(:,2);
obs_cal = C_C0_tracer;

NSE_cal = 1 - nansum((F_shifted(:,2) - C_C0_tracer).^2)./nansum((C_C0_tracer - nanmean(C_C0_tracer)).^2);
A = [F_shifted(:,2),C_C0_tracer];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_cal = corr(1,2);
RMSE_cal = sqrt(nanmean((C_C0_tracer - F_shifted(:,2)).^2));

% ----------------------------------------------------------------------- %
% Results:
% ----------------------------------------------------------------------- %

% figure(1)
% pfit = line(t, C_C0_tracer);
%     set(pfit, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor', [.7 .7 .7])
% hold on
%     hfit = line(t, C_C0_fit);
%     set(hfit, 'LineWidth', 2, 'Color', [0.2 0.3 0.5]);        
% set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
% legend('Data','LSQ')
% hylabel = ylabel('C/C_0');   
% hxlabel = xlabel('Time (s)');
%     set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
% %print -dpng GPFM_CUM_FIT.png
% %close
 
figure(2)
pfit = line(t, C_C0_tracer);
    set(pfit, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor', [.7 .7 .7])
hold on
    hfit = line(t, F);
    set(hfit, 'LineWidth', 2, 'Color', [0.2 0.3 0.5]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
hylabel = ylabel('C/C_0');   
hxlabel = xlabel('Time (s)');
    set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
hold on


% ----------------------------------------------------------------------- %
%%%%%% CDE %%%%%%%%

Dt = 210;  %s
C0 = 0.85; %mg/cm^3
 
% Newton-Raphson's algorithm 
syms p
 
f = @(p) 0.5.*C0.*(erfc((x - v.*t(find(t>Dt):end))./sqrt(4.*p.*t(find(t>Dt):end))) + exp(v.*x./p).*erfc((x+v.*t(find(t>Dt):end))./(sqrt(4.*p.*t(find(t>Dt):end))))) - 0.5.*C0.*(erfc((x - v.*(t(find(t>Dt):end)-Dt))./sqrt(4.*p.*(t(find(t>Dt):end) - Dt))) + exp(v.*x./p).*erfc((x+v.*(t(find(t>Dt):end) - Dt))./(sqrt(4.*p.*(t(find(t>Dt):end) - Dt))))) - C_tracer(find(t>Dt):end);  
f_prime = diff(f,p,1); 
 
fun = 0.5.*C0.*(erfc((x - v.*t(find(t>Dt):end))./sqrt(4.*p.*t(find(t>Dt):end))) + exp(v.*x./p).*erfc((x+v.*t(find(t>Dt):end))./(sqrt(4.*p.*t(find(t>Dt):end))))) - 0.5.*C0.*(erfc((x - v.*(t(find(t>Dt):end)-Dt))./sqrt(4.*p.*(t(find(t>Dt):end) - Dt))) + exp(v.*x./p).*erfc((x+v.*(t(find(t>Dt):end) - Dt))./(sqrt(4.*p.*(t(find(t>Dt):end) - Dt))))) - C_tracer(find(t>Dt):end);  

P_new =  7; 
f_prime_eval = eval(subs(f_prime,p,P_new));
P_new = P_new - nanmean(eval(subs(fun,p,P_new))./f_prime_eval);
tolerance = 0.1;
estimate = 1;
 
for i=1:20
    if estimate> tolerance
        if f_prime_eval ~= 0
        P_old = P_new;
        P_new = P_old - nanmean(eval(subs(fun,p,P_old))./f_prime_eval);
        estimate =  abs(P_new - P_old);
        f_prime_eval = eval(subs(f_prime,p,P_new));
        j = i+1;
        end
    end
end
 
P_CDE = P_new;
    
    %P = 16; C0 = 1; v= 0.625; 
    C_C0_fit = 0.5.*(erfc((x - v.*t)./sqrt(4.*P_CDE.*t))+exp(v.*x./P_CDE).*erfc((x+v.*t)./(sqrt(4.*P_CDE.*t))));
    C_C0_delay = zeros(length(t),1);
  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
    C_C0_delay(i) = 0.5*(erfc((x - v*(t(i)- Dt))/sqrt(4*P_CDE*(t(i)-Dt))) + exp(v*x/P_CDE)*erfc((x+v*(t(i)- Dt))/(sqrt(4*P_CDE.*(t(i) - Dt)))));
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;
 
t_shift = 0; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

sim_cal_CDE = F_shifted(:,2);
obs_cal_CDE = C_C0_tracer;
 
NSE_CDE = 1 - nansum((F_shifted(:,2) - C_C0_tracer).^2)./nansum((C_C0_tracer - nanmean(C_C0_tracer)).^2);
A = [F_shifted(:,2),C_C0_tracer];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_CDE = corr(1,2);
RMSE_CDE = sqrt(nanmean((C_C0_tracer - F_shifted(:,2)).^2));

    hfit = line(F_shifted(:,1), F_shifted(:,2));
    set(hfit,'LineStyle',':', 'LineWidth', 2, 'Color', [0.3010, 0.7450, 0.9330]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
hylabel = ylabel('C/C_0');   
hxlabel = xlabel('Time (s)');
    set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
    legend('Measurements','GPFM', 'CDE')
hold off
% ----------------------------------------------------------------------- %
% TESTING:
% ----------------------------------------------------------------------- %
%%%%% GPFM %%%%%
%%% 1st peak %%%

BTC = [];
BTC = xlsread('Data_Spectro.xlsx',4); % 1) t, 2) C, 3) V per fraction (m), 4) Flow rate (ml/s)

% % parameters for 10Jun19a (first peak fitting)
t = BTC(1:end,1); % (h)
C_tracer = BTC(:,2);
x = 50; % (cm)
v = 1.1; % (cm/s) velocity of the solute
Dt = 11;  %s

%inflection = find(BTC(:,1)==72);
%  P = P_new;
    C_C0_fit = 0.5.*erfc((x-v.*t)./(sqrt(4.*P.*t))) - exp(v.*x.*(1- sqrt(1-4.*P.*eta./v.^2))./(2.*P)-eta.*t).*erfc((x-v.*t.*sqrt(1-4.*P.*eta./v.^2))./(sqrt(4.*P.*t)))./P;
    C_C0_delay = zeros(length(t),1);
  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
            C_C0_delay(i) =0.5*erfc((x-v*(t(i)-Dt))/(sqrt(4*P*(t(i)-Dt)))) - exp(v*x*(1- sqrt(1-4*P*eta/v^2))/(2*P)-eta*(t(i)-Dt))*erfc((x-v*(t(i)-Dt)*sqrt(1-4*P*eta/v^2))/(sqrt(4*P*(t(i)-Dt))))/P;
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;

t_shift = 10; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

end_peak_1 = 27;

sim_test1_pk1_GPFM = F_shifted(1:end_peak_1,2);
obs_test1_pk1_GPFM = C_C0_tracer(1:end_peak_1);

NSE_test = 1 - nansum((F_shifted(1:end_peak_1,2) - C_C0_tracer(1:end_peak_1)).^2)./nansum((C_C0_tracer(1:end_peak_1) - nanmean(C_C0_tracer(1:end_peak_1))).^2);
A = [F_shifted(1:end_peak_1,2),C_C0_tracer(1:end_peak_1)];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_test = corr(1,2);
RMSE_test = sqrt(nanmean((C_C0_tracer(1:end_peak_1) - F_shifted(1:end_peak_1,2)).^2));
  
% figure(3)
% pfit = line(F_shifted(:,1), C_C0_tracer);
%     set(pfit, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor', [.7 .7 .7])
% hold on
%     hfit = line(F_shifted(:,1), C_C0_fit(1:length(F_shifted)));
%     set(hfit, 'LineWidth', 2, 'Color', [0.2 0.3 0.5]);        
% set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
% legend('Data','LSQ')
% hylabel = ylabel('C/C_0');   
% hxlabel = xlabel('Time (s)');
%     set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
% %print -dpng GPFM_CUM_TEST.png
% %close
 
figure(4)
pfit = line(F_shifted(:,1), C_C0_tracer);
    set(pfit, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor', [.7 .7 .7])
hold on
    hfit = line(F_shifted(1:end_peak_1,1), F_shifted(1:end_peak_1,2));
    set(hfit, 'LineWidth', 2, 'Color', [0.2 0.3 0.5]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
    hold on
% ----------------------------------------------------------------------- %
%%%%%%%% CDE %%%%%%%%%%%%
%%%% 1st peak %%%%%%%%%%%
v = 1.1; % (cm/s) velocity of the solute
Dt = 11; %s
% C0 =1.35; %mg/L
P_CDE = 5.67;
 C_C0_fit = 0.5.*(erfc((x - v.*t)./sqrt(4.*P_CDE.*t))+exp(v.*x./P_CDE).*erfc((x+v.*t)./(sqrt(4.*P_CDE.*t))));
    C_C0_delay = zeros(length(t),1);
  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
    C_C0_delay(i) = 0.5*(erfc((x - v*(t(i)- Dt))/sqrt(4*P_CDE*(t(i)-Dt))) + exp(v*x/P_CDE)*erfc((x+v*(t(i)- Dt))/(sqrt(4*P_CDE.*(t(i) - Dt)))));
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;

t_shift = 0; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

end_peak_1 = 27; %28

sim_test1_pk1_CDE = F_shifted(1:end_peak_1,2);
obs_test1_pk1_CDE = C_C0_tracer(1:end_peak_1);

NSE_test_CDE = 1 - nansum((F_shifted(1:end_peak_1,2) - C_C0_tracer(1:end_peak_1)).^2)./nansum((C_C0_tracer(1:end_peak_1) - nanmean(C_C0_tracer(1:end_peak_1))).^2);
A = [F_shifted(1:end_peak_1,2),C_C0_tracer(1:end_peak_1)];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_test_CDE = corr(1,2);
RMSE_test_CDE = sqrt(nanmean((C_C0_tracer(1:end_peak_1) - F_shifted(1:end_peak_1,2)).^2));

hfit = line(F_shifted(1:end_peak_1,1), F_shifted(1:end_peak_1,2));
    set(hfit,'LineStyle',':',  'LineWidth', 2, 'Color', [0.3010, 0.7450, 0.9330]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
    hold on
%-------------------------------------------------------------------------%
%%%%% GPFM %%%%%%%
%%%% 2nd peak %%%%

BTC = [];
BTC = xlsread('Data_Spectro.xlsx',5); % 1) t, 2) C, 3) V per fraction (m), 4) Flow rate (ml/s)
%BTC = xlsread('BTC_DATA_WETLAND.xlsx');

% parameters for 10Jun19a (second peak fitting)
t = [];
t = BTC(1:end,1); % (h)
C_tracer = [];
C_tracer = BTC(:,2);
x = 50; % (cm)
v = 0.625; % (cm/s) velocity of the solute
Dt = 14;  %s

%inflection = find(BTC(:,1)==72);
%  P = P_new;
    C_C0_fit = 0.5.*erfc((x-v.*t)./(sqrt(4.*P.*t))) - exp(v.*x.*(1- sqrt(1-4.*P.*eta./v.^2))./(2.*P)-eta.*t).*erfc((x-v.*t.*sqrt(1-4.*P.*eta./v.^2))./(sqrt(4.*P.*t)))./P;
    C_C0_delay = zeros(length(t),1);
  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
            C_C0_delay(i) =0.5*erfc((x-v*(t(i)-Dt))/(sqrt(4*P*(t(i)-Dt)))) - exp(v*x*(1- sqrt(1-4*P*eta/v^2))/(2*P)-eta*(t(i)-Dt))*erfc((x-v*(t(i)-Dt)*sqrt(1-4*P*eta/v^2))/(sqrt(4*P*(t(i)-Dt))))/P;
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;

t = t + 36;

t_shift = -7; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

begin_peak_2 = 13; %13

sim_test1_pk2_GPFM = F_shifted(begin_peak_2:end,2);
obs_test1_pk2_GPFM = C_C0_tracer(begin_peak_2:end);

NSE_test2 = 1 - nansum((F_shifted(begin_peak_2:end,2) - C_C0_tracer(begin_peak_2:end)).^2)./nansum((C_C0_tracer(begin_peak_2:end) - nanmean(C_C0_tracer(begin_peak_2:end))).^2);
A = [F_shifted(begin_peak_2:end,2),C_C0_tracer(begin_peak_2:end)];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_test2 = corr(1,2);
RMSE_test2 = sqrt(nanmean((C_C0_tracer(begin_peak_2:end) - F_shifted(begin_peak_2:end,2)).^2));
  

%pfit = line(F_shifted(:,1), C_C0_tracer);
%     set(pfit, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor', [.7 .7 .7])
 hold on
    hfit = line(F_shifted(begin_peak_2:end,1), F_shifted(begin_peak_2:end,2));
    set(hfit, 'LineWidth', 2, 'Color', [0.2 0.3 0.5]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
    set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
hold on

% ----------------------------------------------------------------------- %
%%%%%%% CDE %%%%%%%%%%%
%%%%%%%% 2nd peak %%%%%
t = BTC(1:end,1); % (h)
v = 0.625; % (cm/s) velocity of the solute
Dt = 14; %s
%C0 =1.35; %mg/L

 C_C0_fit = 0.5.*(erfc((x - v.*t)./sqrt(4.*P_CDE.*t))+exp(v.*x./P_CDE).*erfc((x+v.*t)./(sqrt(4.*P_CDE.*t))));
    C_C0_delay = zeros(length(t),1);
  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
    C_C0_delay(i) = 0.5*(erfc((x - v*(t(i)- Dt))/sqrt(4*P_CDE*(t(i)-Dt))) + exp(v*x/P_CDE)*erfc((x+v*(t(i)- Dt))/(sqrt(4*P_CDE.*(t(i) - Dt)))));
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;

t = t + 36;

t_shift = -14;%-20; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

begin_peak_2 = 11;%9;

sim_test1_pk2_CDE = F_shifted(begin_peak_2:end,2);
obs_test1_pk2_CDE = C_C0_tracer(begin_peak_2:end);

NSE_test2_CDE = 1 - nansum((F_shifted(begin_peak_2:end,2) - C_C0_tracer(begin_peak_2:end)).^2)./nansum((C_C0_tracer(begin_peak_2:end) - nanmean(C_C0_tracer(begin_peak_2:end))).^2);
A = [F_shifted(begin_peak_2:end,2),C_C0_tracer(begin_peak_2:end)];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_test2_CDE = corr(1,2);
RMSE_test2_CDE = sqrt(nanmean((C_C0_tracer(begin_peak_2:end) - F_shifted(begin_peak_2:end,2)).^2));

hfit = line(F_shifted(begin_peak_2:end,1), F_shifted(begin_peak_2:end,2));
    set(hfit,'LineStyle',':', 'LineWidth', 2, 'Color', [0.3010, 0.7450, 0.9330]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
hylabel = ylabel('C/C_0');   
hxlabel = xlabel('Time (s)');
    set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
    legend('Measurements','GPFM','CDE')
    
    
% ----------------------------------------------------------------------- %
% TESTING:
% ----------------------------------------------------------------------- %
%%%%% GPFM %%%%%
%% no PF
BTC = [];
BTC = xlsread('Data_Spectro.xlsx',3); % 1) t, 2) C, 3) V per fraction (m), 4) Flow rate (ml/s)

% % parameters for 10Jun19a (first peak fitting)
t = BTC(1:end,1); % (h)
C_tracer = BTC(:,2);
x = 50; % (cm)
v = 0.625 ; % (cm/s) velocity of the solute
C0 = 0.85;
Dt = 55 %s

%inflection = find(BTC(:,1)==72);
%  P = P_new;
    C_C0_fit = 0.5.*erfc((x-v.*t)./(sqrt(4.*P.*t))) - exp(v.*x.*(1- sqrt(1-4.*P.*eta./v.^2))./(2.*P)-eta.*t).*erfc((x-v.*t.*sqrt(1-4.*P.*eta./v.^2))./(sqrt(4.*P.*t)))./P;
    C_C0_delay = zeros(length(t),1);
  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
            C_C0_delay(i) =0.5*erfc((x-v*(t(i)-Dt))/(sqrt(4*P*(t(i)-Dt)))) - exp(v*x*(1- sqrt(1-4*P*eta/v^2))/(2*P)-eta*(t(i)-Dt))*erfc((x-v*(t(i)-Dt)*sqrt(1-4*P*eta/v^2))/(sqrt(4*P*(t(i)-Dt))))/P;
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;

t_shift = 0; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

sim_test2_GPFM = F_shifted(:,2);
obs_test2_GPFM = C_C0_tracer(:);

NSE_test_no_PF = 1 - nansum((F_shifted(:,2) - C_C0_tracer(:)).^2)./nansum((C_C0_tracer(:) - nanmean(C_C0_tracer(:))).^2);
A = [F_shifted(:,2),C_C0_tracer(:)];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_test_no_PF = corr(1,2);
RMSE_test_no_PF = sqrt(nanmean((C_C0_tracer(1:end_peak_1) - F_shifted(1:end_peak_1,2)).^2));
  
% figure(3)
% pfit = line(F_shifted(:,1), C_C0_tracer);
%     set(pfit, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor', [.7 .7 .7])
% hold on
%     hfit = line(F_shifted(:,1), C_C0_fit(1:length(F_shifted)));
%     set(hfit, 'LineWidth', 2, 'Color', [0.2 0.3 0.5]);        
% set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
% legend('Data','LSQ')
% hylabel = ylabel('C/C_0');   
% hxlabel = xlabel('Time (s)');
%     set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
% %print -dpng GPFM_CUM_TEST.png
% %close
 
figure(5)
pfit = line(F_shifted(:,1), C_C0_tracer);
    set(pfit, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor', [.7 .7 .7])
hold on
    hfit = line(F_shifted(:,1), F_shifted(:,2));
    set(hfit, 'LineWidth', 2, 'Color', [0.2 0.3 0.5]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
    hold on
% ----------------------------------------------------------------------- %
%%%%%%%% CDE %%%%%%%%%%%%
v = 0.625; % (cm/s) velocity of the solute
Dt = 55; %s
 C0 =0.85; %mg/L

 C_C0_fit = 0.5.*(erfc((x - v.*t)./sqrt(4.*P_CDE.*t))+exp(v.*x./P_CDE).*erfc((x+v.*t)./(sqrt(4.*P_CDE.*t))));
    C_C0_delay = zeros(length(t),1);
  
    for i=1:length(t)
        if t(i)<Dt
            C_C0_delay(i) = 0;
        else
    C_C0_delay(i) = 0.5*(erfc((x - v*(t(i)- Dt))/sqrt(4*P_CDE*(t(i)-Dt))) + exp(v*x/P_CDE)*erfc((x+v*(t(i)- Dt))/(sqrt(4*P_CDE.*(t(i) - Dt)))));
        end
    end
    
F = C_C0_fit - C_C0_delay;
C_C0_tracer = C_tracer./C0;

t_shift = 0; % s
F_shifted = [t-t_shift,F];
F_shifted(any(F_shifted<0,2),:)=[];
C_C0_tracer = C_C0_tracer(1:length(F_shifted));

sim_test2_CDE = F_shifted(:,2);
obs_test2_CDE = C_C0_tracer(:);

NSE_test_CDE_no_PF = 1 - nansum((F_shifted(:,2) - C_C0_tracer(:)).^2)./nansum((C_C0_tracer(:) - nanmean(C_C0_tracer(:))).^2);
A = [F_shifted(:,2),C_C0_tracer(:)];
A(any(isnan(A),2),:)=[];
corr = corrcoef(A(:,1),A(:,2));
R_test_CDE_no_PF = corr(1,2);
RMSE_test_CDE_no_PF = sqrt(nanmean((C_C0_tracer(:) - F_shifted(:,2)).^2));

hfit = line(F_shifted(:,1), F_shifted(:,2));
    set(hfit,'LineStyle',':',  'LineWidth', 2, 'Color', [0.3010, 0.7450, 0.9330]);% [0.3010, 0.6450, 0.6330]);        
set(gca, 'Fontname', 'Helvetica','Fontsize',7, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.01 .01], 'XMinorTick', 'on', 'YMinorTick', 'on', 'Ygrid', 'on', 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3],  'LineWidth', 1);
  legend('Measurements','GPFM','CDE')
  hylabel = ylabel('C/C_0');   
hxlabel = xlabel('Time (s)');
    set([hxlabel, hylabel], 'FontSize', 10, 'FontWeight', 'bold');
