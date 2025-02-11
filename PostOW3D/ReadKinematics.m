% This script reads the unformatted binary kinematics output file from the
% OceanWave3D code.
clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the number of bits and the file name to correspond to your system.
Nbits = 32; % Set Nbits=32 or 64

% Paths
dirp = 'SE/';
casep = 'LAE/';
runp = '20/';

jullib = '/home/andreasp/WAVEGEM/library/';
post = '/postOW3D/';
decomp = 'Decomposition/';

WriteOut = 'yes'; % 'yes','no': Save variables to text files
Pressure='yes'; %'yes','no': Compute time-derivatives and pressures
Plots='no'; %'yes','no': Plot data
%%%%%%%%%%%%%%%%%%%%%%% END OF INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idn = 1:5 % Set idn= kinematics file # to read
  for phsp = 0:3
    run_folder = [dirp,casep,runp,num2str(phsp)]; % Folder of OCW3D run results
    save_fld = [jullib,run_folder,post]; % Path to Julia lib folder
    save_fld2 = [jullib,dirp,casep,runp,decomp]; % Path to Decomposition


    %% INPUT CHECKS %%
    if exist('Nbits')==0 | exist('idn')==0
        error('Set "Nbits"=32 or 64 and "idn"= kinematics file number to read');
    end

    % Open the file
    if idn<10
        fname=[run_folder,'/Kinematics0',num2str(idn),'.bin'];
    else
        fname=[run_folder,'/Kinematics',num2str(idn),'.bin'];
    end
    fid1 = fopen(fname); % File ID number for kinematics data

    % Check for problems
    if fid1 == -1
        error(['Could not open file ',fname,', returned fid1 = - 1']);
    end

    % Choose number of bits
    if Nbits == 32
        int_nbit = 'int';
        njunkreak = 2;
    elseif Nbits == 64
        int_nbit = 'int64';
        njunkread = 2;
    else
        disp(['Illegal value for Nbits, Nbits = ' int2str(Nbits)]);
    end

    %% SCRIPT %%
    % Read the data from the file
    % These read statements must correspond exactly to what appears in the
    % Fortran subroutine: <top dir>/src/IO/StoreKinematicData.f90
    %
    junk = fread(fid1,1,int_nbit); % Read as a junk parameter
    xbeg = fread(fid1,1,'int'); %
    xend = fread(fid1,1,'int'); %
    xstride = fread(fid1,1,'int'); %
    ybeg = fread(fid1,1,'int'); %
    yend = fread(fid1,1,'int'); %
    ystride = fread(fid1,1,'int'); %
    tbeg = fread(fid1,1,'int'); %
    tend = fread(fid1,1,'int'); %
    tstride = fread(fid1,1,'int'); %
    dt = fread(fid1,1,'double'); % Time step size
    nz = fread(fid1,1,'int'); %

    junk = fread(fid1,2,int_nbit); % Junk read statements are necessary for eol markers

    nx=floor((xend-xbeg)/xstride)+1; ny=floor((yend-ybeg)/ystride)+1; nt=floor((tend-tbeg)/tstride)+1;
    %
    % A scratch vector for reading the data
    %
    tmp=zeros(nx*ny*max(nz,5),1);
    %
    % The x-y grid, the depth and bottom gradients for this slice of data
    %
    tmp(1:5*nx*ny)=fread(fid1,5*nx*ny,'double');
    junk = fread(fid1,2,int_nbit);
    %
    x=zeros(nx,ny); x(:)=tmp(1:5:5*nx*ny);
    y=zeros(nx,ny); y(:)=tmp(2:5:5*nx*ny);
    h=zeros(nx,ny); h(:)=tmp(3:5:5*nx*ny);
    hx=zeros(nx,ny); hx(:)=tmp(4:5:5*nx*ny);
    hy=zeros(nx,ny); hy(:)=tmp(5:5:5*nx*ny);
    %
    % The sigma coordinate
    %
    for i=1:nz
        sigma(i)=fread(fid1,1,'double');
    end
    junk = fread(fid1,2,int_nbit);
    %
    % Initialize arrays for the solution on this slice
    %
    eta=zeros(nt,nx,ny); etax=zeros(nt,nx,ny); etay=zeros(nt,nx,ny);
    phi=zeros(nt,nz,nx,ny); w=zeros(nt,nz,nx,ny);
    u=zeros(nt,nz,nx,ny); uz=zeros(nt,nz,nx,ny); ux=u; uy=u;
    v=zeros(nt,nz,nx,ny); vz=zeros(nt,nz,nx,ny); vx=v;vy=v;
    wz=zeros(nt,nz,nx,ny); wx=w; wy=w;
    t=[0:nt-1]*dt*tstride;   % The time axis
    %
    % Read in the solution variables eta, gradeta, phi, u, v, w, dudz, dvdz.
    %
    for it=1:nt
        try
            tmp(1:nx*ny)=fread(fid1,nx*ny,'double');
            eta(it,:)=tmp(1:nx*ny);
            junk = fread(fid1,2,int_nbit);
        catch
            warning(['Read failed at time step ',num2str(it)]);
            break
        end
        %
        tmp(1:nx*ny)=fread(fid1,nx*ny,'double');
        etax(it,:)=tmp(1:nx*ny);
        junk = fread(fid1,2,int_nbit);
        %
        tmp(1:nx*ny)=fread(fid1,nx*ny,'double');
        etay(it,:)=tmp(1:nx*ny);
        junk = fread(fid1,2,int_nbit);
        %
        tmp=fread(fid1,nx*ny*nz,'double');
        phi(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        tmp=fread(fid1,nx*ny*nz,'double');
        u(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        tmp=fread(fid1,nx*ny*nz,'double');
        v(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        tmp=fread(fid1,nx*ny*nz,'double');
        w(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        wz(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        wx(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        wy(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        uz(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        ux(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        uy(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        vz(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        vx(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
        %
        [tmp,count]=fread(fid1,nx*ny*nz,'double');
        % Check for an incomplete run.
        if count==0, it=it-1; break; end
        vy(it,:)=tmp;
        junk = fread(fid1,2,int_nbit);
    end
    display(['Read ',num2str(it),' data points out of ',num2str(nt)])
    %%
    % Write final surface elevation
    %% Surface elevation (spatial at t(end)) %%
    if idn == 5
      fout = 'eta_x@tend';
      fileID = fopen([save_fld,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:xend
          fprintf(fileID, '%.6f\t%.6f\n', x(i), eta(i));
      end
      fclose(fileID);
      Pressure='no';
    endif
    %%
    switch Pressure
        case 'yes'
            %
            % Compute the pressure and fluid acceleration from the standard
            % output kinematics.  This is only done along the first slice in y
            % for 3D problems.
            %
            % Build the 4th-order, even grid time-differentiation matrix
            %
            alpha=3; r=2*alpha+1; c=BuildStencilEven(alpha,1);
            % Put in the centered schemes everywhere:
            Dt=spdiags([ones(nt,1)*c(:,alpha+1)'],[-alpha:alpha],nt,nt);
            % Correct the boundary points with one-sided schemes:
            for j=1:alpha
                Dt(j,:)=0; Dt(j,1:r)=c(:,j)';
                Dt(nt-j+1,:)=0; Dt(nt-j+1,nt-r+1:nt)=c(:,r-j+1)';
            end
            Dt=Dt/dt; % Scale with dt.
            %
            % Compute the Eulerian time-derivatives of eta, phi, u, w, the
            % second time-derivative of eta, and the pressure.
            %
            etat=zeros(nt,nx); etatt=etat; phit=zeros(nt,nz,nx); p=phit; ut=p;
            wt=p;

            for ip=1:nx
                etat(:,ip)=Dt*eta(:,ip); etatt(:,ip)=Dt*etat(:,ip);
                %
                for j=1:nz
                    phit(:,j,ip)=Dt*phi(:,j,ip)-w(:,j,ip).*sigma(j).*etat(:,ip);
                    p(:,j,ip)=-(phit(:,j,ip)+1/2*(u(:,j,ip).^2+w(:,j,ip).^2));
                    ut(:,j,ip)=Dt*u(:,j,ip)-uz(:,j,ip).*sigma(j).*etat(:,ip);
                    wt(:,j,ip)=Dt*w(:,j,ip)-wz(:,j,ip).*sigma(j).*etat(:,ip);
                end
            end
    end
    %%
    %test = zeros(nt,1);
    %for i = 2:nt-1
    %  test(i) = (u(i+1,end,1)-u(i-1,end,1))/(t(i+1)-t(i-1));
    %end
    %figure()
    %plot(t, test); hold on;
    %plot(t,ut(:,end,1),'r--')
    % A block for plotting the results
    %
    switch Plots
        case 'yes'
            %
            % ip=round(nx/2)+1; jp=1; % The horizontal grid point position to plot out
            % itp=[nt-33:2:nt]; % The time slice to focus on
            ip=1; jp=1; % The horizontal grid point position to plot out
            itp=[it-20:it]; % The time slice to focus on
            %
            % Plot the elevation as a function of time, then space at the last time
            % step.
            %
            ifig=1; figure(ifig);clf;
            plot(t,eta(:,ip,jp)); xlabel('t');
            ylabel(['eta(t) at x=',num2str(x(ip,jp))]);
            ifig=ifig+1; figure(ifig);clf;
            plot(x(:,1),eta(end,:,1),x(:,1),etax(end,:,1));
            xlabel('x'); title(['Free surface at t=',num2str(t(end))]);
            legend('eta(x)','eta_x(x)');
            %
            % This is the vertical coordinate for a linear problem:
            %
            z=sigma*h(ip,jp)-h(ip,jp);
            %
            %
            ifig=ifig+1; figure(ifig);clf;
            plot(phi(itp,:,ip,jp),sigma,'-+');
            title(['Profiles of phi(sigma) at x=',num2str(x(ip,jp))] );
            ylabel('sigma');
            %
            ifig=ifig+1; figure(ifig);clf;
            plot(u(itp,:,ip,jp),sigma,'-+');
            title(['Profiles of u(sigma) at x=',num2str(x(ip,jp))] );
            ylabel('\sigma');
            %
            ifig=ifig+1; figure(ifig);clf;
            plot(w(itp,:,ip,jp),sigma,'-+');
            title(['Profiles of w(sigma) at x=',num2str(x(ip,jp))] );
            ylabel('sigma');
            %
            ifig=ifig+1; figure(ifig);clf;
            plot(uz(itp,:,ip,jp),sigma,'-+');
            title(['Profiles of du/dz(sigma) at x=',num2str(x(ip,jp))] );
            ylabel('sigma');
    end

    %% Save variables to text files
    switch WriteOut
        case 'yes'
    %% Write to Decomposition folder %%
      if ~exist(save_fld2, 'dir')
        mkdir(save_fld2)
      endif

      %% Surface elevation for decomposition %%
      fout = ['/eta',num2str(phsp)];
      if ~exist([save_fld2,fout], 'file') && idn==4
        fileID = fopen([save_fld2,fout], 'w');
        if fileID == -1 % Check if file opened successfully
            error('Could not open file for writing');
        end
        for i = 1:tend % Write vector to the file
            fprintf(fileID, '%.6f\t%.6f\n', t(i), eta(i,1));
        end
        fclose(fileID); % Close the file
      endif

    %% Write to postprocessing folder %%
    % Wave probes
    % Old runs
    %  if idn == 1
    %    fldpost = [save_fld,'x0'];
    %  elseif idn == 3
    %    fldpost = [save_fld,'x1000'];
    %  elseif idn == 4
    %    fldpost = [save_fld,'x1600'];
    %  elseif idn == 7
    %    fldpost = [save_fld,'x100'];
    %  elseif idn == 8
    %    fldpost = [save_fld,'x400'];
    %  elseif idn == 9
    %    fldpost = [save_fld,'x700'];
    %  end
    % New runs
        if idn == 1
        fldpost = [save_fld,'x0'];
      elseif idn == 2
        fldpost = [save_fld,'x100'];
      elseif idn == 3
        fldpost = [save_fld,'x700'];
      elseif idn == 4
        fldpost = [save_fld,'x1000'];
      end

      if ~exist(fldpost, 'dir')
        mkdir(fldpost)
      endif

      %% Surface elevation %%
      fout = '/eta';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), eta(i,1));
      end
      fclose(fileID);

      %% Time derivative of surface elevation %%
      fout = '/etat';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), etat(i,1));
      end
      fclose(fileID);

      %% 2nd time derivative of surface elevation %%
      fout = '/etatt';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), etatt(i,1));
      end
      fclose(fileID);

      %% Slope %%
      fout = '/etax';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), etax(i,1));
      end
      fclose(fileID);

      %% Horizontal velocity %%
      fout = '/u';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), u(i,end,1));
      end
      fclose(fileID);

      %% Horizontal acceleration %%
      fout = '/ut';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), ut(i,end,1));
      end
      fclose(fileID);

      %% Vertical acceleration %%
      fout = '/wt';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), wt(i,end,1));
      end
      fclose(fileID);

      %% Pressure %%
      fout = '/p';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), p(i,end,1));
      end
      fclose(fileID);

      %% Potential %%
      fout = '/phi';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), phi(i,end,1));
      end
      fclose(fileID);

      %% Potential flux %%
      fout = '/phit';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end
      for i = 1:tend
          fprintf(fileID, '%.6f\t%.6f\n', t(i), phit(i,end,1));
      end
      fclose(fileID);

      %% Profile of horizontal velocity %%
      fout = '/u_prf';
      fileID = fopen([fldpost,fout], 'w');
      if fileID == -1
          error('Could not open file for writing');
      end

      for i = 2:length(sigma)
        for j = tend-20:tend
          fprintf(fileID, '%.6f\t', u(j,i,1));
        end
        fprintf(fileID, '\n');
      end
      fclose(fileID);
    end
  end
end

