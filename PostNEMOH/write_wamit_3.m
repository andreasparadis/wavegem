function [FeNorm] = write_wamit_3(w,Fe,WAMITULEN,fpath,runname)
  % Non-dimensionalisation of excitation forces
  % Output in WAMIT format for OpenFAST
  rho = 1025; g = 9.81; Amp = 1;
  Nw = length(w);

  fname = [fpath,runname,'.3']
  fid = fopen(fname, 'w')

  for k = 1:Nw
    for i = 1:6
      if le(i,3);
        FeNorm(k,i) = Fe(k,i)/(rho*g*Amp*WAMITULEN^2);
      elseif ge(i,4);
        FeNorm(k,i) = Fe(k,i)/(rho*g*Amp*WAMITULEN^3);
      endif
      line = [2*pi/w(k) 0 i abs(FeNorm(k,i)) angle(FeNorm(k,i))*180/pi real(FeNorm(k,i)) imag(FeNorm(k,i))];
      fprintf(fid, '%e \t %e \t %d \t %e \t %e \t %e \t %e \n', line);
    endfor
  endfor

  fclose(fid)
end

