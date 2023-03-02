%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sample file for Water Waves modelling lab work
%
%% To be modified
%
% G. Ducrozet / F. Bonnefoy - October, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%RUN SECTION BY SECTION PLEASE

function WaterWaves_LabWork()
%
clear all
close all
%
addpath('./functions')
%
%% Define the question/case you solve
test_functions = 1;
irr            = 0;
large_scale    = 0;
wave_structure = 0;
%
% Check if linear solution is studied
if (test_functions == 1)
    %
    % Linear dispersive waves - Airy wave solution
    %
    g = 9.81;              % Define the gravity acceleration
    d = 10;                % Define the water depth
    %
    % Study of dispersion relation k=f(omega)
    %
    omega=linspace(0.01,2,50); % Define the different omega
    %
    k = kfromw(omega,d,g);
    %
    % Draw a figure representing k=f(omega)
    %
    figure;
    plot(omega,k)
    xlabel('\omega (rad/s)')
    ylabel('k (m^-^1)')
    % If you want to change all fonts in current figure
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
    %
    % Phase and group velocity
    %
    c_phase = omega./k;
    c_group = c_phase/2.*(1+2*k*d./sinh(2*k*d));
    c_phase_deep = sqrt(k.^(-1).*g);
    
    mu = k.*d;  %dispersion parameter
    c_phase_ad = c_phase./c_phase_deep;
    c_group_ad = c_group./c_phase_deep;
    %
    figure;
    plot(omega,c_phase,omega,c_group,'g--')
    xlabel('\omega (rad/s)')
    ylabel('Velocity (m/s)')
    legend('c_{phase}','c_{group}')
    % If you want to change all fonts in current figure
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
    
    figure;
    plot(mu, c_phase_ad, mu, c_group_ad)
    xlabel ('Dispersion Parameter')
    ylabel ('Adimensional Velocities')
    legend ('c phase-ad','c group-ad')
    
    %%  ECN basin
    %time span without reflection
    L = 50; %basin length
    B = 30; %basin width
    D = 5; %basin depth
    
    xM = 20; %distance of measuring point M from wavemakers
    T_b = linspace(0.5,4,8); %range of wave periods
    omega_b = 2*pi./T_b;
    
    k_b = kfromw(omega_b,D,g);
    
    c_phase_b = omega_b./k_b;
    c_group_b = c_phase_b./2.*(1+2*k_b*D./sinh(2*k_b*D));
    
    tmin = xM./c_group_b;
    tmax = (2*L-xM)./c_group_b;
    
    %wave amplitudes
    epsilon_b = 0.1;    %fixed steepness parameter
    A_b = epsilon_b./k_b;
    
    %Question 10
    
    epsilon_b1=0.01;
    epsilon_b2=0.02;
    epsilon_b3=0.03;
    epsilon_b4=0.04;
    epsilon_b5=0.05;
    A_b1 = epsilon_b1./k_b;
    A_b2 = epsilon_b2./k_b;
    A_b3 = epsilon_b3./k_b;
    A_b4 = epsilon_b4./k_b;
    A_b5 = epsilon_b5./k_b;
    
    %figure of amplitude as function of time
    plot(T_b,A_b1,T_b,A_b2,T_b,A_b3,T_b,A_b4,T_b,A_b5)
    xlabel('Time period (s)')
    ylabel('Amplitude (m)')
    legend('e=1%','e=2%','e=3%','e=4%','e=5%')
    title('Amplitude as function of time')
    
    %Question 9
    
    plot(T_b,tmin,T_b,tmax)
    xlabel('Time(s)')
    ylabel('Time window(s)')
    legend('tmin','tmax')
    title('Time window plot for ocean water tank')
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Study of linear solution of dispersive waves (Airy wave)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %clear k omega c_phase c_group
    
    %Question 11
    
    % gravity and water depth defined before...
    A   = 10;     % Define amplitude of the wave
    T   = 4;      % Define period of the wave or wavelength
    rho = 1000;   % Define water density
    t_b =linspace(0,100,500);
    figure;
    for i=1:500
    omega  = 2*pi/T;
    k      = kfromw(omega,d,g); % Evaluate the corresponding wavenumber to previous period T
    lambda_b(4) = 2*pi/k_b(4);
    %
    % Free surface elevation as function of space
    %
    x = linspace(0,3*lambda_b(4),100);  % Describe the x-domain
    %
    eta = elevation_linear(A_b(4),k_b(4),D,g,x,t_b(i));    % Describe eta(x)
    %
    
    %
    % First, define the grid on which the quantity has to be constructed
    %
    z = linspace(-D,0,25);         % z-vector (between bottom and FS)
    %
    % We need to construct a full grid (i.e. 2D grid)
    %
    [X,Z] = meshgrid(x,z);  % Note that Matlab is case-sensitive!
    % %
    % % Velocity potential at a given time 
    % %
    %phi = velocity_potential_linear(A,k,h,g,X,Z,t0);
    %
    % Velocity components and pressure
    %
    [U,W,P] = velocity_pressure_linear(A_b(4),k_b(4),D,g,rho,X,Z,t_b(i));
    %
    % Representation of a 2D field
    %
    %figure;
    contourf(X,Z,P)
    xlabel('x (m)')
    ylabel('z (m)')
    h2=colorbar;
    set(get(h2,'ylabel'),'String','Dynamic pressure P (N/m^2)');
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
    %
    % Representation of a 2D vector field with components (U,W)
    %
    %figure;
    hold on;
    quiver(X,Z,U,W)
    xlabel('x (m)')
    ylabel('z (m)')
    %
    % Add free surface elevation if needed
    plot(x,eta)
    ylim([min(z),A_b(4)])
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
    hold off
    pause(0.001)     
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Irregular wave... TO COMPLETE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%QUESTION 12 & 13

if (irr == 1) %test not to compute at every run of the function
    clear all
    %
    % Define parameters
    g = 9.81;
    d = 10;
    rho=1025;
    % Define the frequency vector
    f = linspace(0.1,2,200);
    %
    % Define JONSWAP spectrum parameters
    Hs  = 5;
    Tp  = 8;
    gam = linspace(1,5,5);
    X=2;
    phase = rand(size(f))*2*pi;
    S=zeros(200,5);
    t=linspace(0,50,100);
    eta=zeros(size(gam,2),size(t,2));
    for i=1:size(gam,2)
        S(:,i) = JONSWAP_spectrum(f,Hs,Tp,gam(i));
            for j=2:size(f,2)
                omega=2*pi*f(j);
                k=kfromw(omega,d,g);
                a=sqrt(2*S(j,i)*(f(j)-f(j-1)));
                eta(i,:)=eta(i,:)+elevation_linear(a,k,d,g,X,t,phase(j))
            end

        figure(i)
        subplot(2,1,1)
        plot(f,S(:,i))
        xlabel('Frequency (in Hz)')
        ylabel('Wave spectrum in (m^2/Hz)')
        title('gamma=',gam(i))
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
        subplot(2,1,2)
        plot(t,eta(i,:))
        title('gamma=',gam(i))
    end

    %QUESTION 14
    
    gam=5;
    S = JONSWAP_spectrum(f,Hs,Tp,gam);
    x=linspace(0,100,100);
    x_five=linspace(0,100,200);
    z=linspace(0,-d,100);
    eta_five=zeros(1,200);
    [X,Z] = meshgrid(x,z);
    Ugam=zeros(100,100);
    Wgam=zeros(100,100);
    Pgam=zeros(100,100);
    for i=2:size(f,2)
        omega=2*pi*f(i);
        k=kfromw(omega,d,g);
        a=sqrt(2*S(i)*(f(i)-f(i-1)));
        for j=1:size(x_five,2)
            eta_five(j)=eta_five(j)+elevation_linear(a,k,d,g,x_five(j),t(10),phase(i));
        end
        [u,w,p]=velocity_pressure_linear(a,k,d,g,rho,X,Z,t(10),phase(i));
        Ugam=Ugam+u;
        Wgam=Wgam+w;
        Pgam=Pgam+p;
    end
    
    figure
    subplot(4,1,1)
    plot(x_five,eta_five)
    title('gamma=5')
    subplot(4,1,2)
    contourf(X,Z,Ugam)
    title('Ugam')
    subplot(4,1,3)
    contourf(X,Z,Wgam)
    title('Wgam')
    subplot(4,1,4)
    contourf(X,Z,Pgam)
    title('Pgam')
    

    %%
    % In this function, f=0 should be excluded
    S = JONSWAP_spectrum(f,Hs,Tp,gam);
    %
    figure;
    plot(f,S)
    xlabel('Frequency (in Hz)')
    ylabel('Wave spectrum in (m^2/Hz)')
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
    %
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wave-structure interactions... TO DO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

if wave_structure==1
    %% Large Bodies

    % Define parameters
    g = 9.81;
    d = 50;
    nt = 100;   % number of points in time domain 

    % QUESTION 16

    % To load the transfer functions without forward speed
    TF_0nds = importdata('RAO-0nds.dat',' ',3); % space between columns and 3 headerlines
    TF_0nds = TF_0nds.data;
    ti=(TF_0nds(:,1))';  % periods interval, from data
    
    % frequency domain corresponing to given period
    for i = 1:length(ti)
        for j=2:size(TF_0nds,2)
            TF_0nds(i,1) = 1/ti(length(ti)-i+1);
            TF_0nds(i,j) = TF_0nds(length(ti)-i+1,j);
        end
    end

    f =  TF_0nds(:,1);

    % Defining a spectrum for realistic case
    % swell spectrum parameters
    gams = 10;  % gamma
    Hss = 4;    % significant wave height
    Tps = 15;   % peak period
    % swell spectrum
    Ss = JONSWAP_spectrum(f,Hss,Tps,gams);
    % wind spectrum parameters
    gamw = 2;   %gamma
    Hsw = 3;    % significant wave height
    Tpw = 5;    % peak period
    % wind spectrun
    Sw = JONSWAP_spectrum(f,Hsw,Tpw,gamw);
    
    % total spectrum
    S = Ss;

   % ship response in heave with Wiener-Khichin relation 
   % H = RAO comes from data

   Sy=zeros(18,1);
   for i=1:length(f)
       Sy(i)=((TF_0nds(i,4)))^2*S(i);
   end

   % plots
   figure 
   plot(f,S)
   grid off
   title('Response spectrum without forward speed')
   xlabel('Frequency (in Hz)')
   ylabel('Wave spectrum in (m^2/Hz)')
   set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
   hold on


   %% QUESTION 16 & 17
   L=50;        % length of the ship
   Fr=0.22;

   % RAOs
   chi = 160*pi/180;
   if chi==0
    TF_Fr0p22 = importdata('RAO-Froude0p22.dat',' ',3); %gives you 0 deg (following sea)
    TF_Fr0p22 = TF_Fr0p22.data;
   elseif chi == 135*pi/180
    TF_Fr0p22 = importdata('RAO-Froude0p22.dat',' ',60); %gives you 135 deg (head sea)
    TF_Fr0p22 = TF_Fr0p22.data;
   elseif chi==pi
    TF_Fr0p22 = importdata('RAO-Froude0p22.dat',' ',79); %gives you 180 deg (head sea)
    TF_Fr0p22 = TF_Fr0p22.data;       
   end

   for i = 1:length(ti)
        for j=2:size(TF_0nds,2)
            TF_Fr0p22(i,1) = 1/ti(length(ti)-i+1);
            TF_Fr0p22(i,j) = TF_Fr0p22(length(ti)-i+1,j);
        end
   end

   % ship speed, from Froude similarity
   U=Fr*sqrt(g*L);

   % encounter frequency
   k = kfromw(2*pi*f,d,g)';
   fe = f - U*k*cos(chi);

   % linear interpolation to get RAOs at encounter frequency
   TF_Fr0p22e = interp1(f,TF_Fr0p22(:,4),fe,'pchip');

   % response spectrum from RAOs
   for i=1:length(fe)
       Sh(i) = S(i)/abs(1-4*pi*f(i)*U*cos(chi)/g);
       Sy_sp(i)=((TF_Fr0p22e(i)))^2*Sh(i);
   end
    
   % plot input spectra
   plot(fe,Sh)
   grid on
   xlabel('Frequency (in Hz)')
   ylabel('Wave spectrum in (m^2/Hz)')
   set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
   legend('Without forward speed',strcat('U=', num2str(U),' m/s; \chi=', num2str(chi*180/pi),'°'))
   title('Input Spectrum')
   hold off

   figure
   plot(f,Sy)
   grid on
   xlabel('Frequency (in Hz)')
   ylabel('Ship Response Spectrum Sz (m^2/Hz)')
   set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
   hold on
   plot(fe,Sy_sp)
   grid on
   legend('Without forward speed',strcat('U=', num2str(U),' m/s; \chi=', num2str(chi*180/pi),'°'))
   title('Response Spectrum')
   hold off

   % ship response in time domain
   X = 2;   % fixed location
   phase = 2*pi*rand(length(f));

   t = linspace(0,nt,80);
   eta_w = zeros(1,length(t));
   eta_s0nds = zeros(1,length(t));
   eta_sFr = zeros(1,length(t));

   for i = 2:length(f)

       % wave elevation
       w_w(i) = 2*pi*f(i);
       k_w(i) = kfromw(w_w(i),d,g);
       a_w(i) = sqrt(2*S(i)*(f(i)-f(i-1)));
       eta_w = eta_w + elevation_linear(a_w(i),k_w(i),d,g,X,t,phase(i));

       % ship motion amplitude, no speed
       w_s0nds(i) = 2*pi*f(i);
       k_s0nds(i) = kfromw(w_s0nds(i),d,g);
       a_s0nds(i) = sqrt(2*Sy(i)*(f(i)-f(i-1)));
       eta_s0nds = eta_s0nds + elevation_linear(a_s0nds(i),k_s0nds(i),d,g,X,t,phase(i)); 

       % ship motion amplitude, with speed
       w_sFr0(i) = 2*pi*fe(i);
       k_sFr0(i) = kfromw(w_sFr0(i),d,g);
       a_sFr0(i) = sqrt(2*Sy_sp(i)*(fe(i)-fe(i-1)));
       eta_sFr = eta_sFr + elevation_linear(a_sFr0(i),k_sFr0(i),d,g,X,t,phase(i)); 

   end

   % plot
   figure
   plot(t,eta_w,t,eta_sFr)
   grid on
   title('Temporal signal of the body moment')
   legend('with forward speed in head sea','without forward speed',strcat('U=', num2str(U),' m/s; \chi=', num2str(chi*180/pi),'°'))
   xlabel('Time (s)')
   ylabel('Body displacement (m)')


end
