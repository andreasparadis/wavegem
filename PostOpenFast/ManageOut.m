%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global paths
OFastHome = fullfile("/run","media","data_storage","ONISILOS","OFAST")
addpath(genpath(fullfile(OFastHome,"matlab-toolbox")))
OFastOut = pwd();

% Flags
wout = true    % Write output files
WTid = 2;      % Which wind turbine? (1: Semi-submersible, 2: Mike's FOWT)
WaveCase = 5;  % Case as in Hydrodyn/WAVES (if not 0 or 3, make changes below)
SScase = 'LAE'; % Sea-state case (Depending on Hs,Tp and short or long simul.)
SSrun = '11';  % Sea-state run under specified case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify paths depending on simulation type
if WTid == 1
  fWT = '1_SemiSub'
  OutFile = '5MW_OC4Semi.outb';
elseif WTid == 2
  fWT = '2_Mikes'
  OutFile = 'Mike.outb';
end

%% Which simulated case (where is ~/out/fWT/.../***.outb)?
if WaveCase == 0 % Still water/Initial displacement
    Runfld = 'InitDisp'
elseif WaveCase == 3 % White noise input
    Runfld = 'WhiteNoise'
elseif WaveCase == 5 % Externally generated surface elevation
    Runfld = fullfile('SE',SScase,SSrun)
else
    Runfld = 'JS_8m_10s'
end
% For wave event simulations
if WaveCase == 5
  EEv = ''; % '', MaxFair_#, MaxPitch_#, MaxCoM_#, MaxWave_#
  SimCase = ''; % '', 'event', 'FWG', 'DAM', '2AM', 'SFWG', 'DWG' or 'NW'
  JulPath = fullfile('/home','andreasp','WAVEGEM','library', Runfld, '0','postOFAST',EEv);
else
  EEv = '';
  SimCase = '';
  JulPath = fullfile('/home','andreasp','WAVEGEM','library', 'SE', 'OFAST_FOWT',fWT,Runfld);
end

%% Which outputs?
FASTfilesDesc = {"Time","PtfmSurge", "PtfmHeave", "PtfmPitch", "PtfmYaw",...
              "Wave1Elev", "Wave1Elv1", "FAIRTEN1", "FAIRTEN2",...
              "ANCHTEN1", "ANCHTEN2", "HydroFxi", "HydroFzi", "HydroMyi"}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define necessary paths
FOWTdir = fullfile(OFastOut,fWT);
Casefld = fullfile(Runfld, EEv, SimCase);
OutDir = fullfile(FOWTdir, Casefld);

%% Move simulation output files to the appropriate output folder
suffixes = {".outb", ".out", ".sum", ".ech", ".chkp"};% Suffixes of output files
InDir = fullfile(OFastHome,fWT)
OutFiles = dir(InDir); % List of all files in input directory

% Loop through files and move only those matching suffixes
for i = 1:length(OutFiles)
    filename = OutFiles(i).name;
    if ~OutFiles(i).isdir  % Ignore directories
        for j = 1:length(suffixes)
            if endsWith(filename, suffixes{j})  % Check suffix
                srcFile = fullfile(InDir, filename);
                destFile = fullfile(OutDir, filename);
                movefile(srcFile, destFile);
                break; % Move to the next file
            end
        end
    end
end


FASTfiles = {fullfile(OutDir, OutFile)}
Channels = FASTfilesDesc;
ReferenceFile = 1
ShowLegend = 1
CustomHdr = {[], 2, 1, 2}
PlotPSDs = 1
OnePlot = 0

%% Call the plotting function
[outData, FR, PSDs]=PlotFASToutput(FASTfiles,FASTfilesDesc,ReferenceFile,...
                              Channels,ShowLegend,CustomHdr,PlotPSDs,OnePlot);

t = outData{1,:};

%% Plot RAOs
% Colors
dred = [0.7 0.2 0];
amber = [0.9 0.6 0];
dgr = [0 0.7 0];

figure();
subplot(1,3,1)
loglog(1./FR, PSDs(:,2),'LineWidth',2,'color',dred)
xlabel("T [s]")
ylabel("RAO")
xlim([1 1000])
legend("Surge")
grid on

subplot(1,3,2)
loglog(1./FR, PSDs(:,3),'LineWidth',2,'color',amber)
xlabel("T [s]")
ylabel("RAO")
xlim([1 1000])
legend("Heave")
grid on

subplot(1,3,3)
loglog(1./FR, PSDs(:,4),'LineWidth',2,'color',dgr)
xlabel("T [s]")
ylabel("RAO")
xlim([1 1000])
legend("Pitch")
grid on

if WaveCase == 0
  eigfreq = zeros(1,3);
  % Damping ratios
  [pm, ipm] = max(PSDs(2:30,2));
  omn = 2*pi*FR(ipm+1); eigfreq(1) = FR(ipm+1);
  zeta = 0.005;
  env = max(abs(outData{1,2}(:,2))) * exp(-zeta*omn*t);

  figure()
  plot(t, outData{1,2}(:,2),'k-','LineWidth',2); hold on;
  plot(t, env,'r--','LineWidth',2);

  [pm, ipm] = max(PSDs(2:end,3));
  omn = 2*pi*FR(ipm+1); eigfreq(2) = FR(ipm+1);
  zeta = 0.02;
  env = max(abs(outData{1,2}(:,3))) * exp(-zeta*omn*t);

  figure()
  plot(t, outData{1,2}(:,3),'k-','LineWidth',2); hold on;
  plot(t, env,'r--','LineWidth',2);

  [pm, ipm] = max(PSDs(2:end,4));
  omn = 2*pi*FR(ipm+1); eigfreq(3) = FR(ipm+1);
  zeta = 0.03;
  env = max(abs(outData{1,2}(:,4))) * exp(-zeta*omn*t);

  figure()
  plot(t, outData{1,2}(:,4),'k-','LineWidth',2); hold on;
  plot(t, env,'r--','LineWidth',2);

  fname = fullfile(OutDir,'FR');
  save("-text", fname, "FR");

  fname = fullfile(OutDir,'PSDs');
  save("-text", fname, "PSDs");

  fid = fopen(fullfile(OutDir,"eigen"), 'w')
  fprintf(fid, '%f \t %f \t %f \t \n', eigfreq);
  fprintf(fid, '%f \t %f \t %f \t \n', 1./eigfreq);
  fprintf(fid, '%f \t %f \t %f \t \n', 2*pi*eigfreq);
  fclose(fid)
end

if wout
  % Save output data to txt
  results = outData{1,2};
  fname = fullfile(JulPath, "outD_", SimCase)
  save("-text", fname, "results");

  fid = fopen(fullfile(OutDir,"dataHdr"), 'w')
  heads = Channels.';
  fprintf(fid, '%s \t', heads{:});
  fclose(fid)
end
