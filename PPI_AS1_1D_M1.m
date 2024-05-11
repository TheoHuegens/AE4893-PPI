clear;
close all;
clc;
figure

R_mc_tests = [1790e3,1860e3,1985e3,2090e3];
cols =  ['b','r','k','y'];
line_type = ["-", "--", ":","-."];

%% input parameters of your planet

% constants
G = 6.6743e-11;

% parameters
M_p = 0.330103e24; %+-0.000021e24
R_p = 2439.4e3; %+-0.1
ImIc_meas = 0.431; %+-0.025
V_p = (4/3)*pi*R_p^3;
rho_bulk = M_p/V_p;
MOI_ratio_meas =  0.346; %+-0.014

for test = 1:length(R_mc_tests)
    R_mc = R_mc_tests(test);
    disp(['Model ' num2str(test) ' with R_m/c=' num2str(R_mc/1e3), ' [km]']);
   
    
    A = [
            [(R_p^3-R_mc^3), R_mc^3],
            [(R_p^5-R_mc^5), R_mc^5]
         ];
    B = [
            [(3/(4*pi))*M_p],
            [(15/(8*pi))*MOI_ratio_meas*M_p*R_p^2]
        ];
    
    X = linsolve(A,B);
    
    density_profile = [ % 2-layer
        [0,X(2)],
        [R_mc,X(1)],
        [R_p,0]
        ];
    disp(['rho_core=' num2str(X(2)), ' [kg/m3]']);
    disp(['rho_mantle=' num2str(X(1)), ' [kg/m3]']);
    disp('----------------------------------------------')
    
    % %% Compute bulk characteristics of density model
    MOI = 0;
    M = 0;
    for layer = 1:length(density_profile)-1
        rho = density_profile(layer,2);
        r = density_profile(layer,1);
        R = density_profile(layer+1,1);
        MOI = MOI + (8/15)*pi* rho * (R^5 - r^5);
        M = M + (4/3)*pi* rho * (R^3 - r^3);
    end
    disp(['Measured Mass is ' num2str(M_p) ' -'])
    disp(['Model Mass is ' num2str(M) ' -'])
    disp(['Difference is ' num2str((M-M_p)/M_p*100) ' percent'])
    disp('----------------------------------------------')
    MOI_ratio = MOI / (M * R_p^2);
    disp(['Measured I/MR^2 is ' num2str(MOI_ratio_meas) ' -'])
    disp(['Model I/MR2 is ' num2str(MOI_ratio) ' -'])
    disp(['Difference is ' num2str((MOI_ratio_meas-MOI_ratio)/MOI_ratio_meas*100) ' percent'])
    disp('----------------------------------------------')
    
    %% Constructing the upward integration loop
    
    % initialise your variables
    Dr = 10;
    M0 = 0;
    r0 = 0;
    g0 = 0;
    MOI = 0;
    MOI_c = 0;
    Interior = zeros(length(1:R_p/Dr),5); % radius, mass, gravity, pressure, density
    
    % perform the upward loop
    layer = 1;
    for ii = 1:R_p/Dr
        r  = r0 + Dr;
        
        if r > density_profile(layer+1,1)
            layer = layer + 1; 
        end
        rho_0 = density_profile(layer,2);
        Interior(ii,5) = rho_0;
    
        M  = M0 + 4*pi*rho_0*r^2*Dr;
        MOI = MOI + (8/15)*pi* rho_0*((r+Dr/2)^5 - (r-Dr/2)^5);
        if r < R_mc
            MOI_c = MOI;
        end

        g = G*M/r^2;
    
        Interior(ii,1) = r;
        Interior(ii,2) = M;
        Interior(ii,3) = g;
    
        M0 = M;
        r0 = r;
    end
    
    % initialise
    p0 = 0;
    r0 = R_p;
    
    % perform the downward loop
    for ii = 1:R_p/Dr
    
        r  = r0 - Dr;
        if r <= density_profile(layer-1,1)
            layer = layer - 1;
        end
        rho_0 = density_profile(layer,2);
        p = p0 + rho_0*Interior(R_p/Dr-(ii-1),3)*Dr;
        p0 = p;
        Interior(R_p/Dr-(ii-1),4) = p;
    end

    disp(['pressure at CMB ' num2str(Interior(int64(R_mc/Dr)+1,4)/10^9) 'GPa']);
    disp(['pressure at core ' num2str(Interior(1,4)/10^9) 'GPa']);
    disp('----------------------------------------------')
    
    %% verification
    disp(['Measured mass is ' num2str(M_p) ' kg'])
    disp(['Numerical mass is ' num2str(M) ' kg'])
    disp(['Difference is ' num2str((M_p-M)/M_p*100) ' percent'])
    disp('----------------------------------------------')
    MOI_ratio = MOI / (M * R_p^2);
    disp(['Measured I/MR^2 is ' num2str(MOI_ratio_meas) ' -'])
    disp(['Numerical I/MR2 is ' num2str(MOI_ratio) ' -'])
    disp(['Difference is ' num2str((MOI_ratio_meas-MOI_ratio)/MOI_ratio_meas*100) ' percent'])
    disp('----------------------------------------------')
    Ic = MOI_c;
    Im = MOI-Ic;
    ImIc = Im/MOI;
    disp(['Measured Im/Ic is ' num2str(ImIc_meas) ' -'])
    disp(['Numerical Im/Ic is ' num2str(ImIc) ' -'])
    disp(['Difference is ' num2str((ImIc-ImIc_meas)/ImIc_meas*100) ' percent'])
    disp('==============================================')
    
    %% plotting the 1D profiles of a planet with homogeneous denisty
    
    bb = 18;
    
    % plot figure
    
    subplot(1,4,1)
    plt1 = plot(Interior(:,2)./1e23,Interior(:,1)./1e3,'linestyle',line_type(test), 'Color', cols(1), 'LineWidth', 2);
    hold on;
    xlabel('Integrated Mass (10^{23} kg)','FontSize',bb)
    ylabel('Radius (km)','FontSize',bb)
    set(gca,'FontSize',bb)

    subplot(1,4,2)
    plt2 = plot(Interior(:,3),Interior(:,1)./1e3,'linestyle',line_type(test), 'Color', cols(2), 'LineWidth', 2);
    hold on;
    xlabel('Gravity (m/s^2)','FontSize',bb)
    ylabel('Radius (km)','FontSize',bb)
    set(gca,'FontSize',bb)

    subplot(1,4,3)
    plt3 = plot((Interior(:,4))./1e9,Interior(:,1)./1e3,'linestyle',line_type(test), 'Color', cols(3), 'LineWidth', 2);
    hold on;
    xlabel('Pressure (GPa)','FontSize',bb)
    ylabel('Radius (km)','FontSize',bb)
    set(gca,'FontSize',bb)

    subplot(1,4,4)
    plt4 = plot((Interior(:,5)),Interior(:,1)./1e3,'linestyle',line_type(test), 'Color', cols(4), 'LineWidth', 2);
    hold on;
    xlabel('Density (kg/m3)','FontSize',bb)
    ylabel('Radius (km)','FontSize',bb)
    set(gca,'FontSize',bb)
end

legend('','Location', 'best','Orientation','vertical', 'LineWidth', 2)
legend('boxoff')

line_names = ["R_CMB=1200km", "R_CMB=1790km", "R_CMB=1860km", "R_CMB=2040km"];
for j =1:length(line_type)
    plot([NaN NaN], [NaN NaN],line_type(j), 'Color', 'k', 'DisplayName', line_names(j))
end
hold off