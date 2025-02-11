%--------------------------------------------------------------------------------------
% Run NEMOH1 and NEMOH2 (QTF) or just plot existing results
%
% User has to specify the output folder in outputdir input below
%--------------------------------------------------------------------------------------
clc
clear all
close all

pathstr = '/home/andreasp/NEMOH/'
addpath(genpath(pathstr)); % Include the subfolders in the Matlab PATH

ID_RUN = 0;
ID_PLOT_RESULTS = 0;
ID_QTF = 1; % Flag to activate QTF computation (0 or 1)
ID_HydrosCal=1; % Flag to compute hydrostatics, 1 computed, 0 Not
ID_WAM = 0; % Write files in WAMIT format
WAMITULEN = 1;

runname = 'Mike';
outputdir = ['/home/andreasp/NEMOH/',runname];   % Update this output files location
projdir=[outputdir,filesep,''];
ofast_path = ['/run/media/data_storage/ONISILOS/OFAST/HydroData/'];

DOFsurge=1; DOFheave=3; DOFpitch=5;

if ID_RUN == 1
  % Check that Nemoh is available
  assert(FindingNemoh(ID_QTF, true))
  %-------Launch Calculation------------
  [Idw,w,A,B,Fe]=Nemoh(projdir,ID_HydrosCal,ID_QTF); % Call the function Nemoh.m

  %% --- Computes QTFs --------------------
  if ID_QTF==1, NemohQTF(projdir);end % Call the function NemohQTF.m
else
  [Idw,w,A,B,Fe,RAO]=read_results(projdir);
end

%% ---------------------------------------
if ID_PLOT_RESULTS==1 %Plot NEMOH1
    if Idw==1
        xlab='$\omega ~[rad/s]$';
    elseif Idw==2
        xlab='f [Hz]';
    elseif Idw==3
        xlab='T [s]';
    end
    % Colors
    dred = [0.7 0.2 0];
    amber = [0.9 0.6 0];
    dgr = [0 0.7 0];

    % Plot Radiation Coefficients
    h = figure();
    subplot(2,1,1)
    a(1,:)=A(1,1,:);
    loglog(w,a,'color',dred); hold on;
    a(1,:)=A(3,3,:);
    loglog(w,a,'color',amber); hold on;
    a(1,:)=A(5,5,:);
    loglog(w,a,'color',dgr)
      title('\textbf{Radiation Coefficients: Added mass (A) and Damping (B)}')
      legend('Surge','Heave','Pitch')
      ylabel('\textbf{A}')
      xlim([6e-2 10])
      grid ON;  plot_properties;

    subplot(2,1,2)
    b(1,:)=B(1,1,:);
    loglog(w,b,'color',dred); hold on;
    b(1,:)=B(3,3,:);
    loglog(w,b,'color',amber); hold on;
    b(1,:)=B(5,5,:);
    loglog(w,b,'color',dgr)
      xlabel(xlab)
      ylabel('\textbf{B}')
      xlim([6e-2 10])
      grid ON;  plot_properties;
    fname = [projdir,filesep,'Figures',filesep,'RadCoeffs'];
    print(h, fname, "-dpng", "-S1920,1080")
    print(h, fname,"-dpdfcrop", "-svgconvert", "-F:11");

    % Plot Excitation Forces
    h = figure(); hold on;
    loglog(w,abs(Fe(:,1)),'color',dred);
    loglog(w,abs(Fe(:,3)),'color',amber);
    loglog(w,abs(Fe(:,5)),'color',dgr);
      legend('Surge [N]', 'Heave [N]', 'Pitch [Nm]')
      xlabel(xlab)
      ylabel('\textbf{Excitation Forces/Moments}')
      xlim([5e-2 10])
      grid ON;  plot_properties; hold off;
    fname = [projdir,filesep,'Figures',filesep,'ExcForces'];
    print(h, fname, "-dpng", "-S1920,1080")
    print(h, fname,"-dpdfcrop", "-svgconvert", "-F:11");

    % Plot RAOs (Hz)
    freq = w./(2*pi);
    h = figure();
    subplot(2,3,1)
    semilogx(freq, RAO(:,DOFsurge))
      ylabel('$|X|$ [m/m]' )
      legend('Surge')
      xlim([1e-2 1])
      grid on; plot_properties;
    subplot(2,3,4)
    semilogx(freq, RAO(:,6+DOFsurge))
      xlabel('f [Hz]')
      ylabel('ang(x) [deg]' )
      xlim([1e-2 1])
      grid on; plot_properties;
    subplot(2,3,2)
    semilogx(freq, RAO(:,DOFheave))
      ylabel('$|Y|$ [m/m]' )
      legend('Heave')
      title('\textbf{RAO - Magnitude}')
      xlim([1e-2 1])
      grid on; plot_properties;
    subplot(2,3,5)
    semilogx(freq, RAO(:,6+DOFheave))
      xlabel('f [Hz]')
      ylabel('ang(y) [deg]' )
      title('RAO - Phase')
      xlim([1e-2 1])
      grid on; plot_properties;
    subplot(2,3,3)
    semilogx(freq, RAO(:,DOFpitch))
      ylabel('$|\theta|$ [deg]' )
      legend('Pitch')
      xlim([1e-2 1])
      grid on; plot_properties;
    subplot(2,3,6)
    semilogx(freq, RAO(:,6+DOFpitch))
      xlabel('f [Hz]')
      ylabel('ang($\theta$) [deg]' )
      xlim([1e-2 1])
      grid on; plot_properties;
    fname = [projdir,filesep,'Figures',filesep,'RAO'];
    print(h, fname, "-dpng", "-S1920,1080")
    print(h, fname,"-dpdfcrop", "-svgconvert", "-F:11");

    % Plot RAOs (T)
    per = 2*pi./w;
    h = figure();
    subplot(3,1,1)
    plot(per, RAO(:,DOFsurge))
      title('\textbf{RAO - Magnitude}')
      ylabel('$|X|$ [m/m]' )
      legend('Surge')
      xlim([0 100])
      grid on; plot_properties;

    subplot(3,1,2)
    plot(per, RAO(:,DOFheave))
      ylabel('$|Y|$ [m/m]' )
      legend('Heave')
      xlim([0 100])
      grid on; plot_properties;

    subplot(3,1,3)
    plot(per, RAO(:,DOFpitch))
      xlabel('T [s]')
      ylabel('$|\theta|$ [deg]' )
      legend('Pitch')
      xlim([0 100])
      grid on; plot_properties;

    fname = [projdir,filesep,'Figures',filesep,'RAO_s'];
    print(h, fname, "-dpng", "-S1920,1080")
    print(h, fname,"-dpdfcrop", "-svgconvert", "-F:11");

    if ID_QTF==1
        NWdatNEM=100; % Number of freq in QTF output
        Idwvar=0; %Freq type: 0 => w [rad/s], 1 => f [Hz], 2 => T [s]
        SwitchBiDir=0; % Bi-directional flag
        NbetaData=1; % Number of beta values (directions)
        betaID=[0 0]; %[beta1 beta2]
        shiftdw1=-1;  ShowLegend=1;
        Idaxisij=0; % Reverse coordinate system (y values increase from top to bottom)
        % QTF data availability in each mode (all==1 => data available for all DOF)
        DOFdatNem=[1 1 1 1 1 1]; % 0: empty, 1: exists

        % Plot QTFs
        qtftype='M'; % Difference-frequencies QTF
        Id_fixed_CB_Axis=1; % 0: fixed colorbar limit, 1: specific minmaxQ limits
        if Id_fixed_CB_Axis==1
          minmaxQR_surge=[-20 20]; minmaxQI_surge=[-20 20];minmaxQMod_surge=[0 50];
          minmaxQR_heave=[-30 30]; minmaxQI_heave=[-5 5];minmaxQMod_heave=[0 60];
          minmaxQR_pitch=[-500 500];minmaxQI_pitch=[-500 500];minmaxQMod_pitch=[0 2000];
        end
        plot_QTF_NEMOH;

        qtftype='P'; % Sum-frequencies QTF
        Id_fixed_CB_Axis=1; % 0: fixed colorbar limit, 1: specific minmaxQ limits
        if Id_fixed_CB_Axis==1
          minmaxQR_surge=[-100 100]; minmaxQI_surge=[-50 50];minmaxQMod_surge=[0 200];
          minmaxQR_heave=[-80 80]; minmaxQI_heave=[-30 30];minmaxQMod_heave=[0 200];
          minmaxQR_pitch=[-600 600];minmaxQI_pitch=[-700 700];minmaxQMod_pitch=[0 3000];
        end
        plot_QTF_NEMOH;
    end
end

%% Write files in WAMIT format
if ID_WAM == 1;
  % 'runname.hst'
  [KHnorm] = write_wamit_hst(WAMITULEN,ofast_path,runname, outputdir);
  % 'runname.1'
  [Anorm, Bnorm] = write_wamit_1(w,A,B,WAMITULEN,ofast_path,runname);
  % 'runname.3'
  [FeNorm] = write_wamit_3(w,Fe,WAMITULEN,ofast_path,runname);
  % 'runname.12d'
  write_wamit_QTF(projdir,NWdatNEM,DOFdatNem,betaID,NbetaData,'M',...
                  SwitchBiDir,WAMITULEN,ofast_path,runname);
  % 'runname.12s'
  write_wamit_QTF(projdir,NWdatNEM,DOFdatNem,betaID,NbetaData,'P',...
                  SwitchBiDir,WAMITULEN,ofast_path,runname);
end
