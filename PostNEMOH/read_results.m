% Inputs :
% - projdir     : Project directory path
% Outputs :
% - Idw: A frequency type: 1,2,3= rad/s,Hz,s
% - w  : Vector length(w) of wave frequencies (rad/s)
% - A  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of added mass coefficients
% - B  : Matrix (6xnBodies)x(6xnBodies)xlength(w) of radiation damping coefficients
% - Fe : Matrix (6xnBodies)xlength(w) of exciation forces (complex
% values)
%
function [Idw,w,A,B,Fe,RAO]=read_results(projdir)
% Read Nemoh.cal
fid=fopen([projdir,filesep,'Nemoh.cal'],'r');
for i=1:6
    ligne=fgetl(fid);
end
nBodies=fscanf(fid,'%g',1);
for i=1:nBodies
    for ii=1:4
        ligne=fgetl(fid);
    end
    Ndof=fscanf(fid,'%g',1);
    for j=1:Ndof
        ligne=fgetl(fid);
    end
    ligne=fgetl(fid);
    Nforc=fscanf(fid,'%g',1);
    for j=1:Nforc
        ligne=fgetl(fid);
    end
    ligne=fgetl(fid);
end
ligne=fgetl(fid);
ligne=fgetl(fid);
ligne=fscanf(fid,'%g',2);
Idw=ligne(1);
nw=ligne(2);
fclose(fid);

% Read ExcitationForce.tec
fid=fopen([projdir,filesep,'results',filesep,'ExcitationForce.tec'],'r');
ligne=fgetl(fid);
for c=1:Nforc*nBodies
    ligne=fgetl(fid);
end
ligne=fgetl(fid);
for k=1:nw
    ligne=fscanf(fid,'%f',1+2*Nforc*nBodies);
    w(k)=ligne(1);
    for j=1:Nforc*nBodies
        Famp(k,j)=ligne(2*j);
        Fphi(k,j)=ligne(2*j+1);
    end
end
status=fclose(fid);

% Read RadiationCoefficients.tec
fid=fopen([projdir,filesep,'results',filesep,'RadiationCoefficients.tec'],'r');
ligne=fgetl(fid);
for i=1:Ndof*nBodies
    ligne=fgetl(fid);
end
for i=1:nBodies*Ndof
    ligne=fgetl(fid);
    for k=1:nw
        ligne=fscanf(fid,'%f',1+2*Ndof*nBodies);
        for j=1:Ndof*nBodies
            A(i,j,k)=ligne(2*j);
            B(i,j,k)=ligne(2*j+1);
        end
        ligne=fgetl(fid);
    end
end
status=fclose(fid);
% Expression des efforts d excitation de houle sous la forme Famp*exp(i*Fphi)
Fe=Famp.*(cos(Fphi)+1i*sin(Fphi));

% Read RAO.dat
fid = fopen([projdir,filesep,'Motion',filesep,'RAO.dat'],'r');
ligne=fgetl(fid);
ligne=fgetl(fid);
for i = 1:nw
  ligne=fscanf(fid,'%f',1+2*Ndof*nBodies);
  RAO(i,:) = ligne(2:end);
end

end
