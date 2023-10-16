%% Gradient simulator (GradSim) tool for predicting protein elution
% Copyright (c) 2023, Scott H. Altern. All rights reserved.

%% This section calculates experimental Kp values from batch data
% Mass balances are employed to calculate a series of Kp values from
% sequential elution steps performed in the batch experiment. Liquid and
% solid volumes are required, as well as the load protein concentration,
% eluted sample protein concentrations, and salt concentrations.

clear; tic; close all;

% Inputs
Cp_load = 2.2;
Cp = [0.25 0 0 0 0 0.05 0.25 0.25 0.5 0.5 0.2 0.1 0.1];
Cs = [0 150 300 450 600 750 900 1050 1200 1350 1500];

% Experimental volumes
Vr_batch = 40/1000;                      % Resin volume (mL)
Vh_batch = 40/1000;                      % Holdup volume (mL)
Vl_load = 110/1000;                      % Load volume (mL)
Vl_batch = 110/1000;                     % Volume of wash added (mL)

% Store protein concentrations and salt concentrations
Cp_grad(1,:) = Cp(1,2:end-1);
Cp_FT = Cp(1,1);
Cp_strip = Cp(1,end);

% FT check
if Cp_FT > sum(Cp)
    isFT = 1;
else
    isFT = 0;
end

% Strip check
if Cp_strip > sum(Cp)
    isstrip = 1;
else
    isstrip = 0;
end

% Fix zero Cp values to avoid dividing by zero scenario
Cp_grad(Cp_grad <= 0) = 1e-4;
    
% Initialize arrays
Cs_actual = 0*Cs;
Kp = 0*Cp_grad;

% Calculate initial Kp
Kp(1) = (Cp_load*Vl_load - Cp_grad(1)*(Vh_batch + Vl_batch))/...
    (Cp_grad(1)*Vr_batch);

% Calculations for salt and Kp
for i=1:length(Cs)-1
    % Calculate salt concentrations considering hold-up
    Cs_actual(i+1) = (Cs(i+1)*Vl_batch + Cs_actual(i)*Vh_batch)/...
        (Vl_batch + Vh_batch);

    % Calculate Kp for series
    Kp(i+1) = (Cp_grad(i)*(Kp(i)*Vr_batch + Vh_batch) - ...
        Cp_grad(i+1)*(Vh_batch+Vl_batch))/(Cp_grad(i+1)*Vr_batch);
end

% Correct negative Kp values if present
Kp(Kp < 0) = 0;


%% This section predicts protein elution using batch Kp data
% Parameters for salt gradient process are defined by user in addition
% to the protein concentration and appropriate column parameters.
% Salt profile is generated using a column stage-wise mass balance 
% with T time steps and N stages. Simulated Kp profile is interpolated 
% from the salt profile using batch experiment Kp data and corresponding 
% salt steps in the batch sequential elution. Protein profile is calculated 
% using the simulated Kp values throughout the column. Elution salt conc.
% is then calculated from the location of the protein peak maximum, or
% peak first moment. Lastly, the outlet profiles are plotted with volume.

% Setting for elution salt calc.: Peak max (PM) or center of mass (CoM)
moment_setting = 'CoM'; 
 
% Set specifications for gradient process
% Sequence: Load, Wash, Gradient, Hold, Strip
cprot =   [ 1   0   0     0     0  ];    % Protein conc. (mg/ml)
csalt0 =  [ 0   0   0     1.5   0  ];    % Start salt in each sequence [M]
csalt1 =  [ 0   0   1.5   1.5   0  ];    % End salt in each sequence [M]
vol_seq = [ 5   5   40    10    10 ];    % Volumes in each sequence [CV]

% Set column parameters
N  = 10;                                 % Number of stages
Ep = 0.6;                                % Intraparticular porosity
Ee = 0.4;                                % Interstitial porosity
VL = Ee/N;                               % Mobile phase volume
VH = Ep/N;                               % Hold-up volume
Et = Ee + Ep*(1-Ee);                     % Total porosity
VR = (1-Et)/N;                           % Resin volume

% Generate stage and time arrays
vol_seq = [0 vol_seq];
vol_cu = cumsum(vol_seq);
seq = vol_seq;
stages = 1:N;
step_size = 1/N;
timesteps = 0:step_size:sum(seq);
idx = timesteps;

% Generate salt and protein inlet profiles
csalt_in = 0*timesteps; cprot_in = 0*timesteps;
for t=2:numel(timesteps)
    for s=2:numel(seq)
        if timesteps(t) > sum(seq(1:s-1)) && timesteps(t) <= sum(seq(1:s))
            csalt_in(t) = csalt0(s-1) + ((csalt1(s-1)-csalt0(s-1))/...
                seq(s))*(timesteps(t)-sum(seq(1:s-1)));
            cprot_in(t) = cprot(s-1);
            break;
        end
    end
end

% Generate salt and protein column profiles
csalt_col = zeros(length(stages),length(timesteps));
csalt_col(1,:) = csalt_in;
for n=2:numel(stages)                                                                           % Calculates values of salt concentration for all stages and times
    for t=2:numel(timesteps)                                                                       % Using mass balance equation without adsorption
        csalt_col(n,t) = ...
            (csalt_col(n-1,t-1)*VL + csalt_col(n,t-1)*VH)/(VL+VH);                           
    end
end    
csalt_out = csalt_col(end,:);

% Set Kp and Cs
Kp_exp = Kp; Cs_exp = Cs/1000;

% Add zero Kp for case where simulation salt exceeds experiment
if max(csalt1) > max(Cs_exp)
    Cs_exp = [Cs_exp max(csalt1)];
    Kp_exp = [Kp_exp 0];
end

% Generate Kp for salt profile via interpolation 
Kp_sim = pchip(Cs_exp,Kp_exp,csalt_col);

% Generate protein column profile
cprot_col = zeros(length(stages),length(timesteps));
cprot_col(1,:) = cprot_in;
for n=2:length(stages)                                                                           
    for t=2:length(timesteps)                                                                       
        cprot_col(n,t) = (cprot_col(n-1,t-1)*VL + cprot_col(n,t-1)*...
            (VH + Kp_sim(n,t-1)*VR))/(VL + VH + Kp_sim(n,t)*VR);                          
    end
end
cprot_out = cprot_col(end,:);

% Calculate normalized protein trace
[max_prot,loc] = max(cprot_out);
cprot_out_norm = cprot_out/max_prot;

% Determine peak center of mass
if strcmpi(moment_setting,'CoM')
    % Extract linear gradient region
    grad_idx = [find(idx == vol_cu(3)) find(idx == vol_cu(4))];
    cp_int = cprot_out(grad_idx(1):grad_idx(2));
    vol_int = idx(grad_idx(1):grad_idx(2));

    % Cumulative area
    cp_int = cp_int/trapz(vol_int,cp_int);
    area_cu = cumtrapz(vol_int,cp_int);

    % Center of mass location
    loc = find(area_cu >= 0.5,1,'first');
    loc = loc + find(idx == vol_cu(3)) - 1;
end

% Calculate elution salt concentration by checking to see if most
% of the protein was in FT, strip, or gradient
if isFT == 1
    elusalt_sim = 0;
elseif isstrip == 1
    elusalt_sim = 1.501;
else
    elusalt_sim = csalt_out(loc);
end       

% Convert to mM
elusalt_sim = elusalt_sim*1000;

% Print calculated elution salt concentration
fprintf('\nPredicted elution salt concentration = %4.0f mM\n\n',...
    elusalt_sim);

% Plot protein and salt profiles
figure(1)
yyaxis left
plot(idx,cprot_out_norm)
yyaxis right
plot(idx,csalt_out)

% Add plot labels
xlabel("Volume (CV)")
yyaxis left
ylim([-0.025 1.025])
ylabel("Normalized Protein Concentration")
set(gca,'xminortick','on','yminortick','on','fontsize',12)
yyaxis right
ylim([-max(csalt1)*0.025 max(csalt1)*1.025])
ylabel("Salt Concentration (M)")
set(gca,'xminortick','on','yminortick','on','fontsize',12)
grid on
box on

toc