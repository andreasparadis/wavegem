function [Anorm, Bnorm] = write_wamit_1(w,A,B,WAMITULEN,fpath,runname)
  % Non-dimensionalisation of radiation coefficient matrices
  % Output in WAMIT format for OpenFAST
  rho = 1025; Nw = length(w);

  % Anorm: Normalized Added mass matrix
  % Bnorm: Normalized Damping matrix

  fname = [fpath,runname,'.1']
  fid = fopen(fname, 'w')

  for i = 1:6
    for j = 1:6
      for k = 1:Nw
        if and(le(i,3),le(j,3))
          Anorm(i,j,k) = A(i,j,k)/(rho*WAMITULEN^3);
          Bnorm(i,j,k) = B(i,j,k)/(w(k)*rho*WAMITULEN^3);
        elseif and(ge(i,4),ge(j,4))
          Anorm(i,j,k) = A(i,j,k)/(rho*WAMITULEN^5);
          Bnorm(i,j,k) = B(i,j,k)/(w(k)*rho*WAMITULEN^5);
        else
          Anorm(i,j,k) = A(i,j,k)/(rho*WAMITULEN^4);
          Bnorm(i,j,k) = B(i,j,k)/(w(k)*rho*WAMITULEN^4);
        endif
      endfor
    endfor
  endfor

  A0 = Anorm(:,:,1);
  Ainf = Anorm(:,:,end);

  % Write infinite-frequency matrix
  for i = 1:6
    for j = 1:6
      if or((i==j),and((i==1),(j==5)),and((i==5),(j==1)),and((i==2),(j==4)),and((i==4),(j==2)));
        line = [-1 i j Ainf(i,j)];
        fprintf(fid, '%e \t %d \t %d \t %e \n', line);
      endif
    endfor
  endfor

  % Write zero-frequency matrix
  for i = 1:6
    for j = 1:6
      if or((i==j),and((i==1),(j==5)),and((i==5),(j==1)),and((i==2),(j==4)),and((i==4),(j==2)));
        line = [0 i j A0(i,j)];
        fprintf(fid, '%e \t %d \t %d \t %e \n', line);
      endif
    endfor
  endfor

  % Write "non-zero" (significant) radiation coeffciients
  for k = 1:Nw
    for i = 1:6
      for j = 1:6
        if or((i==j),and((i==1),(j==5)),and((i==5),(j==1)),and((i==2),(j==4)),and((i==4),(j==2)));
          line = [2*pi/w(k) i j Anorm(i,j,k) Bnorm(i,j,k)];
          fprintf(fid, '%e \t %d \t %d \t %e \t %e \n', line);
        endif
      endfor
    endfor
  endfor

  fclose(fid)
end

