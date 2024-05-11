
% plot parameters
cols =  ['b','r','k','y','g'];
line_type = ["-", "--", ":","-."];
line_names = ["result (mean P)", "Unc rho_m", "Unc t_s", "Unc Temp"];
test = 1;

% numerical parameters
Dr = 10;
tempITRmax = 5;
varITRmax = 3;

%% input parameters of your planet

% core composition
comp_IC = 2; %1=Fe,liq 2=Fe,sol 3=Fe3S,sol 5=FeSV,sol inner core composition
comp_mantle = 9; % 7=crust 8=FC 9=MA 10=MC mantle composition

% constants
G = 6.6743e-11;

% contraints
M_p = 0.330103e24; %+-0.000021e24
R_p = 2439.4e3; %+-0.1
CMR2 = 0.346; %+-0.014
ImIc_meas = 0.431; %+-0.025

% parameters
rho_m_ref = [2900,3300]; % mantle density
t_s_ref = [50e3,200e3]; % silicate shell thickness
T_core_ref = [2350,2550]; % core-mantle boundary temperature
T_ICB_ref = [2350,2550]; % isothermal core (X,PX)
T_CMB_ref = [1600,2000]; % (1, P6)
T_crust_ref = [0,0]; % interpolate conductivily from CMB to surface
T_surf_ref = [273.15-180,273.15+430]; % night / day

% materials
MatData = readtable("MaterialData.xlsx");%material	T_ref	p_ref	rho_0	drho/dp	drho/dT	drho/dXs	drho/dXs2	K_0	dK/dT	dK/dp	dK/dXs	dK/dXs2	a_0	da/dp	da/dT	Mw	Cp_A	Cp_B	Cp_C	Cp_D	Cp_E	k
disp(MatData.material);

% variables
R_CMB_ref = 2000e3; % 1750-2200 km core-mantle boundary
R_ICB_ref = 1000e3; % 0-1750 km inner core boundary
DR = 100e3;
X_Sulfur_ref = 0.15; % 0-1 [-] Sulfur content % in liquid outer core
DX = 0.05;

Xvars = [R_CMB_ref,R_ICB_ref,X_Sulfur_ref]; % variables (unknowns)

last_ITR = 0; % 0= keep iterating 1=sensitivity analysis on answer 2=stop
varITR = 0;
while (varITR <= varITRmax) & last_ITR <= 1
    close all;
    clc;
    figure

    Phi_der = [
        [0,0,0],
        [0,0,0],
        [0,0,0]
        ];

    if last_ITR == 0
        PhiTests = [
            [Xvars(1),Xvars(2),Xvars(3),"-",2,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1)+DR,Xvars(2),Xvars(3),"--",1,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1)-DR,Xvars(2),Xvars(3),"--",1,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2)+DR,Xvars(3),"-.",1,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2)-DR,Xvars(3),"-.",1,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2),Xvars(3)+DX,":",1,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2),Xvars(3)-DX,":",1,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)]
            ];
    elseif last_ITR == 1
        disp("last ITR");
        PhiTests = [
            [Xvars(1),Xvars(2),Xvars(3),"-",2,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2),Xvars(3),"--",1,rho_m_ref(1),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2),Xvars(3),"--",1,rho_m_ref(2),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2),Xvars(3),"-.",1,mean(rho_m_ref),t_s_ref(1),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2),Xvars(3),"-.",1,mean(rho_m_ref),t_s_ref(2),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)],
            [Xvars(1),Xvars(2),Xvars(3),":",1,mean(rho_m_ref),mean(t_s_ref),T_core_ref(1),T_ICB_ref(1),T_CMB_ref(1),T_crust_ref(1),T_surf_ref(1)],
            [Xvars(1),Xvars(2),Xvars(3),":",1,mean(rho_m_ref),mean(t_s_ref),T_core_ref(2),T_ICB_ref(2),T_CMB_ref(2),T_crust_ref(2),T_surf_ref(2)],
            [Xvars(1),Xvars(2),Xvars(3),"-",2,mean(rho_m_ref),mean(t_s_ref),mean(T_core_ref),mean(T_ICB_ref),mean(T_CMB_ref),mean(T_crust_ref),mean(T_surf_ref)]
            ];
    end
    
    PhiTestResults = [];
    delta_tar = [];
    
    for test = 1:length(PhiTests(:,1 ))
        disp([varITR,test,PhiTests(test,1),PhiTests(test,2),PhiTests(test,3)]);
        % define variables of test
        R_CMB = max(min(str2num(PhiTests(test,1)),R_p-t_s_ref(2)),0); % bounds of CMB within 0>Rp-ts
        R_ICB = max(min(str2num(PhiTests(test,2)),R_CMB),0); % bounds of ICB within 0>CMB
        X_Sulfur = max(min(str2num(PhiTests(test,3)),1),0); % bounds of XS within 0>1

        R_CMB = str2num(PhiTests(test,1));
        R_ICB = str2num(PhiTests(test,2));
        X_Sulfur = str2num(PhiTests(test,3));
        
        % define parameters of test
        rho_m = str2num(PhiTests(test,6));
        t_s = str2num(PhiTests(test,7));
        T_core = str2num(PhiTests(test,8));
        T_ICB = str2num(PhiTests(test,9));
        T_CMB = str2num(PhiTests(test,10));
        T_crust = str2num(PhiTests(test,11));
        T_surf = str2num(PhiTests(test,12));

        linewidth = str2num(PhiTests(test,5));
        
        % create planetary model
        LayerModel = [
            [0,comp_IC,"iso",T_core];         % center: R=0
            [R_ICB,4,"adia",T_ICB];   % transition radius, material, thermal behaviour (T_l,T_u)
            [R_CMB,comp_mantle,"iso",T_CMB];    % transition radius, material, thermal behaviour 
            [R_p-t_s,7,"adia",T_surf];  % transition radius, material, thermal behaviour 
            [R_p,0,"adia",T_surf]           % surface: R=Rp, mat=vacuum
            ];
        LayerBounds = cast(LayerModel(:,1),'double');
        
        tempITR = 0;
        M = 0;
        deltaM = 1;
        
        while (tempITR < tempITRmax) & ( deltaM > 0.02 )
            M_prev = M;
            %% Constructing the upward integration loop
            
            % initialise
            r0 = 0;
            M0 = 0;
            I0 = 0;
            g0 = 0;
            T0 = 0;
            if tempITR == 0
                Interior = zeros(length(1:R_p/Dr),7); % radius, mass, MOI, gravity, pressure, temperature, density
                Interior(1,6) = T_core;
            end
        
            layer = 1;
            
            % perform the upward loop
            for ii = 1:R_p/Dr
                % march up
                r  = r0 + Dr;
            
                % update current layer
                if r > LayerBounds(layer+1) | r0==0
                    if r0 ~= 0
                        layer = layer + 1;
                    end
                    layer_index = str2double(LayerModel(layer,2));
                    %disp(["upward",layer,MatData.material(layer_index)]);

                    if tempITR == 0
                        rho = MatData.rho_0(layer_index)+MatData.drho_dXs(layer_index)*X_Sulfur+MatData.drho_dXs2(layer_index)*X_Sulfur^2;
                    end
                end
            
                if tempITR > 0
                    rho = Interior(ii,7);
                else
                    Interior(ii,7) = rho;
                end
        
                % compute
                M = M0 + 4*pi*rho*r^2*Dr;
                I = I0 + (8/15)*pi* rho*((r+Dr/2)^5 - (r-Dr/2)^5);
                g = G*M/r^2;
            
                Interior(ii,1) = r;
                Interior(ii,2) = M;
                Interior(ii,3) = I;
                Interior(ii,4) = g;
            
                r0 = r;
                M0 = M;
                I0 = I;
            end
            
            %% Constructing the downward integration loop
            
            % initialise
            p0 = 0;
            r0 = R_p;
            
            % perform the downward loop
            for ii = 1:R_p/Dr
                % march down
                r  = r0 - Dr;
            
                % update layer
                if r < LayerBounds(layer) | r==R_p-Dr
                    if r0 ~= R_p 
                        layer = layer - 1;
                    end
                    layer_index = str2double(LayerModel(layer,2));
                    %disp(["downward",layer,MatData.material(layer_index)]);
            
                    % update layer properties
                    T_ref = MatData.T_ref(layer_index);
                    p_ref = MatData.p_ref(layer_index);
                    K_stiff = (MatData.K_0(layer_index) + +MatData.dK_dXs(layer_index)*X_Sulfur+MatData.dK_dXs2(layer_index)*X_Sulfur^2)*1e9;
                    alpha = MatData.a_0(layer_index);
                end
            
                rho = Interior(R_p/Dr-(ii-1),7);
                
                % compute
                p = p0 + rho*Interior(R_p/Dr-(ii-1),4)*Dr;
            
                Interior(R_p/Dr-(ii-1),5) = p;
            
                r0 = r;
                p0 = p;
            end
            
            %% Constructing the thermal loop & linear density model
            % initialise
            r0 = 0;
            
            layer = 1;
            % perform the upward loop
            for ii = 1:R_p/Dr
                % march up
                r  = r0 + Dr;
               
                % update current layer
                if r > LayerBounds(layer+1) | r0==0
                    if r0 ~= 0
                        layer = layer + 1;
                    end
                    layer_index = str2double(LayerModel(layer,2));
                    %disp(["thermal",layer,MatData.material(layer_index),LayerModel(layer,3)]);
            
                    % update layer properties
                    temp_type = LayerModel(layer,3);
            
                    T_ref = MatData.T_ref(layer_index);
                    p_ref = MatData.p_ref(layer_index);
                    rho0 = MatData.rho_0(layer_index) + MatData.drho_dXs(layer_index)*X_Sulfur+MatData.drho_dXs2(layer_index)*X_Sulfur^2;
                    K_stiff = (MatData.K_0(layer_index) +MatData.dK_dXs(layer_index)*X_Sulfur+MatData.dK_dXs2(layer_index)*X_Sulfur^2)*1e9;
                    k_th = MatData.k(layer_index)/(10^3);
                    visc = MatData.visc(layer_index);
        
                    alpha = MatData.a_0(layer_index);% + MatData.da_dp(layer_index)*(Interior(ii,5)-MatData.p_ref(layer_index))+MatData.da_dT(layer_index)*(Interior(ii,6)-MatData.T_ref(layer_index));
                    C_p = (MatData.Cp(layer_index));%+MatData.Cp_B(layer_index)*T+MatData.Cp_C(layer_index)*T^2+MatData.Cp_D(layer_index)*T^3)*1000/MatData.Mw(layer_index);
                   
                    % update layer's thermal properties
                    R_l = str2double(LayerModel(layer,1));
                    R_u = str2double(LayerModel(layer+1,1));

                    T_l = str2double(LayerModel(layer,4));
                    T_u = str2double(LayerModel(layer+1,4));

                    % effectively ignores T_crust and interpolates with
                    % neighbouring boundaries instead
                    if layer == 3
                        R_u = str2double(LayerModel(layer+2,1));
                        T_u = str2double(LayerModel(layer+2,4));
                    elseif layer == 4
                        R_l = str2double(LayerModel(layer-1,1));
                        T_l = str2double(LayerModel(layer-1,4));
                    elseif layer == 5
                        R_l = str2double(LayerModel(layer-2,1));
                        T_l = str2double(LayerModel(layer-2,4));    
                    end

                    rho_avg = Interior(int64((R_l+R_u)/(2*Dr)),7);
                    g_avg = Interior(int64((R_l+R_u)/(2*Dr)),4);

                    if visc ~= 0 % non-solids
                        %disp([num2str(rho_avg)," ",num2str(alpha)," ",num2str(g_avg)," ",num2str(abs(T_u-T_l))," ",num2str(abs(R_u-R_l))," ",num2str(k_th/(rho_avg*C_p))," ",num2str(visc)]);
                        kappa = k_th/(rho_avg*C_p);
                        Ra = (rho_avg*alpha*g_avg*abs(T_u-T_l)*abs(R_u-R_l)^3)/(kappa*visc);
                        CF = (1/((3/(10^4))*Ra+1)); % maps 0-Inf (Rauleigh) to 0-1 (conductive-convective) with halfway convective at Ra=1e4
                    else
                        Ra = 0;
                        CF = 1; % conductive
                    end
                    CF = 1;

                    R_BL = CF*(R_u-R_l)/2; % thickness of boundary layer
                    T_avg = (T_l+T_u)/2; % temp in convection zone
                    
                    %disp([num2str(Ra),num2str(CF)]);
                end
            
                % temperature parameters
                p = Interior(ii,5);
                if ii > 1
                    T = Interior(ii-1,6);
                else
                    T = T_core;
                end
            
                % temperature model
               
                if R_u == R_l
                    dTdr = 0;
                else
                    if R_l < r < R_l+R_BL
                        dTdr = (T_avg-T_l)/R_BL;
                    elseif R_u-R_BL < r < R_u
                        dTdr = (T_u-T_avg)/R_BL;
                    else 
                        dTdr = 0;
                    end
                end
                T = T + dTdr*Dr;
            
                % linear density model
                rho = rho0*(1-alpha*(T-T_ref)+(p-p_ref)/K_stiff);
            
                Interior(ii,6) = T;
                Interior(ii,7) = rho;
            
                r0 = r;
            end
            
            tempITR = tempITR + 1;
            M = Interior(R_p/Dr-1,2);
            deltaM = abs((M - M_prev) / M);
        end

        %% verification
        M = Interior(R_p/Dr-1,2);
        disp(['Measured mass is ' num2str(M_p) ' kg'])
        disp(['Numerical mass is ' num2str(M) ' kg'])
        disp(['Difference is ' num2str((M_p-M)/M_p*100) ' percent'])
        disp('----------------------------------------------')
        I = Interior(R_p/Dr-1,3);
        MOI_ratio = I / (M * R_p^2);
        disp(['Measured I/MR^2 is ' num2str(CMR2) ' -'])
        disp(['Numerical I/MR2 is ' num2str(MOI_ratio) ' -'])
        disp(['Difference is ' num2str((CMR2-MOI_ratio)/CMR2*100) ' percent'])
        disp('----------------------------------------------')
        Ic = Interior(int64(R_CMB/Dr)-1,3);
        Im = I-Ic;
        ImIc = Im/I;
        disp(['Measured Im/Ic is ' num2str(ImIc_meas) ' -'])
        disp(['Numerical Im/Ic is ' num2str(ImIc) ' -'])
        disp(['Difference is ' num2str((ImIc-ImIc_meas)/ImIc_meas*100) ' percent'])
        disp('==============================================')
        
        PhiTestResults(end+1,:) = [M/M_p,MOI_ratio/CMR2,ImIc/ImIc_meas];

        delta_tar(test,:) = [1- PhiTestResults(1,1), 1- PhiTestResults(1,2), 1- PhiTestResults(1,3)];
        
        %% plotting the 1D profiles of a planet with homogeneous denisty

        if last_ITR == 1
        
            bb = 18;
            
            % plot figure
            
            subplot(1,5,1)
            plt1 = plot(Interior(:,2)./1e23,Interior(:,1)./1e3,'linestyle',PhiTests(test,4), 'Color', cols(1), 'LineWidth', linewidth);
            hold on;
            xlabel('Integrated Mass (10^{23} kg)','FontSize',bb)
            ylabel('Radius (km)','FontSize',bb)
            set(gca,'FontSize',bb)
            
            subplot(1,5,2)
            plt2 = plot(Interior(:,4),Interior(:,1)./1e3,'linestyle',PhiTests(test,4), 'Color', cols(2), 'LineWidth', linewidth);
            hold on;
            xlabel('Gravity (m/s^2)','FontSize',bb)
            ylabel('Radius (km)','FontSize',bb)
            set(gca,'FontSize',bb)
            
            subplot(1,5,3)
            plt3 = plot((Interior(:,5))./1e9,Interior(:,1)./1e3,'linestyle',PhiTests(test,4), 'Color', cols(3), 'LineWidth', linewidth);
            hold on;
            xlabel('Pressure (GPa)','FontSize',bb)
            ylabel('Radius (km)','FontSize',bb)
            set(gca,'FontSize',bb)
            
            subplot(1,5,4)
            plt4 = plot((Interior(:,7)),Interior(:,1)./1e3,'linestyle',PhiTests(test,4), 'Color', cols(4), 'LineWidth', linewidth);
            hold on;
            xlabel('Density (kg/m3)','FontSize',bb)
            ylabel('Radius (km)','FontSize',bb)
            set(gca,'FontSize',bb)
            
            subplot(1,5,5)
            plt5 = plot((Interior(:,6)),Interior(:,1)./1e3,'linestyle',PhiTests(test,4), 'Color', cols(5), 'LineWidth', linewidth);
            hold on;
            xlabel('Temperature (K)','FontSize',bb)
            ylabel('Radius (km)','FontSize',bb)
            set(gca,'FontSize',bb)

        end
       
    end
    
    % legend
    legend('','Location', 'best','Orientation','vertical', 'LineWidth', 2)
    legend('boxoff')
    
    for j =1:length(line_type)
        plot([NaN NaN], [NaN NaN],line_type(j), 'Color', 'k', 'DisplayName', line_names(j))
    end
    hold off
    
    if last_ITR == 1
        last_ITR = 2;
    elseif (((delta_tar(end,1) && delta_tar(end,2) && delta_tar(end,3)) < 0.05) | varITR+1 == varITRmax)
        last_ITR = 1;
    end

    % compute derivatives
    Phi_der(1,1) = (PhiTestResults(2,1)-PhiTestResults(3,1)) / (2*DR); % dM/dR_CMB
    Phi_der(2,1) = (PhiTestResults(2,2)-PhiTestResults(3,2)) / (2*DR); % dI/dR_CMB
    Phi_der(3,1) = (PhiTestResults(2,3)-PhiTestResults(3,3)) / (2*DR); % dImIc/dR_CMB
    
    Phi_der(1,2) = (PhiTestResults(4,1)-PhiTestResults(5,1)) / (2*DR); % dM/dR_IC
    Phi_der(2,2) = (PhiTestResults(4,2)-PhiTestResults(5,2)) / (2*DR); % dI/dR_ICB
    Phi_der(3,2) = (PhiTestResults(4,3)-PhiTestResults(5,3)) / (2*DR); % dImIc/dR_ICB
    
    Phi_der(1,3) = (PhiTestResults(6,1)-PhiTestResults(7,1)) / (2*DX); % dM/dXS
    Phi_der(2,3) = (PhiTestResults(6,2)-PhiTestResults(7,2)) / (2*DX); % dI/dXS
    Phi_der(3,3) = (PhiTestResults(6,3)-PhiTestResults(7,3)) / (2*DX); % dImIc/dXS
    
    % adapt variables to minimize measurement difference
    dXvars = linsolve(Phi_der,transpose(delta_tar(end,:)));
    disp(PhiTestResults);
    disp(Phi_der);
    disp(["delta_tar:",delta_tar(1),"% M_p | ",delta_tar(2),"% I/MR2",delta_tar(3),"% Im/I"]);
    disp(["dXvars:",dXvars(1),"dR_CMB",dXvars(2),"dR_ICB",dXvars(3),"dR_XS"]);
    Xvars = Xvars + transpose(dXvars);
    Xvars(1) = max(min(Xvars(1),R_p-t_s_ref(2)),0); % bounds of CMB within 0>Rp-ts
    Xvars(2) = max(min(Xvars(2),Xvars(1)),0); % bounds of ICB within 0>CMB
    Xvars(3) = max(min(Xvars(3),1),0); % bounds of XS within 0>1

    varITR = varITR + 1;
   
     
end

disp('PARAMETERS');
disp(['mantle density: ', num2str(mean(rho_m_ref)), '+-', num2str((max(rho_m_ref)-min(rho_m_ref))/2), ' [kg/m3]']);
disp(['crust thickness: ', num2str(mean(t_s_ref)/1e3), '+-', num2str((max(t_s_ref)-min(t_s_ref))/(2e3)), ' [km]']);
disp(['core temperature: ', num2str(mean(T_core_ref)), '+-', num2str((max(T_core_ref)-min(T_core_ref))/2), ' [K]']);
disp(['surface temperature: ', num2str(mean(T_surf_ref)), '+-', num2str((max(T_surf_ref)-min(T_surf_ref))/2), ' [K]']);
disp('VARIABLES');
disp(['inner-core boundary: ', num2str(R_ICB/1e3),' [km]']);
disp(['core-mantle boundary: ', num2str(R_CMB/1e3),' [km]']);
disp(['sulfur fraction: ', num2str(100*X_Sulfur), ' [%]']);
disp('RESULTS');
disp(['planetary mass:', num2str(M_p), ' +error:', num2str(M_p*(delta_tar(1,1))), '+-', num2str(M_p*max(abs(delta_tar(:,1)))), ' [kg]']);
disp(['planetary mass:', num2str(M), ' +error:', num2str(100*(delta_tar(1,1))), '+-', num2str(max(abs(delta_tar(:,1)))), ' [%]']);
disp(['planetary moment of inertia ratio:', num2str(CMR2), '+error:', num2str((delta_tar(1,2))), '+-', num2str(max(abs(delta_tar(:,2)))), ' [-]']);
disp(['planetary moment of inertia ratio:', num2str(MOI_ratio), '+error:', num2str(100*(delta_tar(1,2))/CMR2), '+-', num2str(max(abs(delta_tar(:,2)))/CMR2), ' [%]']);
disp(['mantle to core moment of inertia ratio:', num2str(ImIc_meas), '+error:', num2str((delta_tar(1,3))), '+-', num2str(max(abs(delta_tar(:,3)))), ' [-]']);
disp(['mantle to core moment of inertia ratio:', num2str(ImIc), '+error:', num2str(100*(delta_tar(1,3))/ImIc_meas), '+-', num2str(max(abs(delta_tar(:,3)))/ImIc_meas), ' [%]']);
disp('VALIDATION');
disp(['pressure at CMB ' num2str(Interior(int64(R_mc/Dr)+1,5)/10^9) 'GPa']);
disp(['pressure at core ' num2str(Interior(1,5)/10^9) 'GPa']);