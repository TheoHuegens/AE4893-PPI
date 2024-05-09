clear;
close all;
clc;
figure

R_mc_tests = [1200e3,1790e3,1860e3,2040e3];
cols =  ['b','r','k','y'];
line_type = ["-", "--", ":","-."];

%% input parameters of your planet

% constants
G = 6.6743e-11;

% validation parameters
M_p = 3.302e23;
rho_bulk = 5427;
MOI_ratio_meas = 0.353;

% input parameters
R_p = 2440e3;

for test = 1:length(R_mc_tests)
    R_mc = R_mc_tests(test);
    disp(['Model ' num2str(test) ' with R_m/c=' num2str(R_mc/1e3), ' [km]']);
    disp('----------------------------------------------')
   
    
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
    
    % %% Compute bulk characteristics of density model
    % MOI = 0;
    % M = 0;
    % for layer = 1:length(density_profile)-1
    %     rho = density_profile(layer,2);
    %     r = density_profile(layer,1);
    %     R = density_profile(layer+1,1);
    %     MOI = MOI + (8/15)*pi* rho * (R^5 - r^5);
    %     M = M + (4/3)*pi* rho * (R^3 - r^3);
    % end
    % disp(['Theoretical Mass is ' num2str(M_p) ' -'])
    % disp(['Model Mass is ' num2str(M) ' -'])
    % disp(['Difference is ' num2str((M-M_p)/M_p*100) ' percent'])
    % disp('----------------------------------------------')
    % MOI_ratio = MOI / (M * R_p^2);
    % disp(['Theoretical I/MR^2 is ' num2str(MOI_ratio_meas) ' -'])
    % disp(['Model I/MR2 is ' num2str(MOI_ratio) ' -'])
    % disp(['Difference is ' num2str((MOI_ratio_meas-MOI_ratio)/MOI_ratio_meas*100) ' percent'])
    % disp('----------------------------------------------')
    
    %% Constructing the upward integration loop
    
    % initialise your variables
    Dr = 10;
    M0 = 0;
    r0 = 0;
    g0 = 0;
    MOI = 0;
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
    
    % perform the upward loop
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
    
    %% verification
    disp(['Measured mass is ' num2str(M_p) ' kg'])
    disp(['Numerical mass is ' num2str(M) ' kg'])
    disp(['Difference is ' num2str((M_p-M)/M_p*100) ' percent'])
    disp('----------------------------------------------')
    MOI_ratio = MOI / (M * R_p^2);
    disp(['Measured I/MR^2 is ' num2str(MOI_ratio_meas) ' -'])
    disp(['Numerical I/MR2 is ' num2str(MOI_ratio) ' -'])
    disp(['Difference is ' num2str((MOI_ratio_meas-MOI_ratio)/MOI_ratio_meas*100) ' percent'])
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

line_names = ["R_m_c=1200 [km]", "R_m_c=1790 [km]", "R_m_c=1860 [km]", "R_m_c=2040 [km]"];
for j =1:length(line_type)
    plot([NaN NaN], [NaN NaN],line_type(j), 'Color', 'k', 'DisplayName', line_names(j))
end
hold off