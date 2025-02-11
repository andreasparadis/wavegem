function [KH_norm] = write_wamit_hst(WAMITULEN,fpath,runname, outdir)
  % Non-dimensionalisation of radiation coefficient matrices
  % Output in WAMIT format for OpenFAST
  rho = 1025;
  g=9.81;

  rfname = [outdir,filesep,'Mechanics',filesep,'Kh.dat']
  KH = dlmread(rfname)

  fname = [fpath,runname,'.hst']
  fid = fopen(fname, 'w')

  for i = 1:6
    for j = 1:6
        if and(le(i,3),le(j,3))
          KH_norm(i,j) = KH(i,j)/(rho*g*WAMITULEN^2);
        elseif and(ge(i,4),ge(j,4))
          KH_norm(i,j) = KH(i,j)/(rho*g*WAMITULEN^4);
        else
          KH_norm(i,j) = KH(i,j)/(rho*g*WAMITULEN^3);
        endif
        line = [i j KH_norm(i,j)];
        fprintf(fid, '%d \t %d \t %e \n', line);
    endfor
  endfor

  fclose(fid)
end

