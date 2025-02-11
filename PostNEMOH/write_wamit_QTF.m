function [] = write_wamit_QTF(projdir,Nw,dof_data,betaId,Nbetadat,qtftype,...
                              SwitchBiDir,WAMITULEN,fpath,runname)
  % Output of QTFs in WAMIT format for OpenFAST

  if qtftype == 'M'
    fname = [fpath,runname,'.12d']
  elseif qtftype == 'P'
    fname = [fpath,runname,'.12s']
  endif
  fid = fopen(fname, 'w')
  fprintf(fid,'%s \n',['NEMOH Numeric Output -- Andreas Paradeisiotis -- ',date]);

  % QTF file name (depending on qtftype)
  QTFfname = [projdir,filesep,'results',filesep,'QTF',filesep,'OUT_QTF',qtftype,'_N.dat'];

  for DOF = 1:6
    [Qreal(:,:,DOF),Qimag(:,:,DOF),Qmod(:,:,DOF),w1]=fun_readNEMOHQTFOutput_N(QTFfname,...
                           Nw,dof_data,DOF,betaId,Nbetadat,qtftype,SwitchBiDir);
  endfor

  for k1 = 1:Nw
    for k2 = k1:Nw
        for DOF = 1:2:5
          if le(DOF,3);
            RQTF = Qreal(k1,k2,DOF)/WAMITULEN;
            IQTF = Qimag(k1,k2,DOF)/WAMITULEN;
            MQTF = Qmod(k1,k2,DOF)/WAMITULEN;
          elseif ge(DOF,4);
            RQTF = Qreal(k1,k2,DOF)/WAMITULEN^2;
            IQTF = Qimag(k1,k2,DOF)/WAMITULEN^2;
            MQTF = Qmod(k1,k2,DOF)/WAMITULEN^2;
          endif
          line = [2*pi/w1(k1) 2*pi/w1(k2) betaId(1) betaId(2)...
                  DOF MQTF atan2(IQTF,RQTF)*180/pi RQTF IQTF];
          fprintf(fid, '%e \t %e \t %e \t %e \t %d \t %e \t %e \t %e \t %e \n', line);
        endfor
        for DOF = 2:2:6
          if le(DOF,3);
            RQTF = Qreal(k1,k2,DOF)/WAMITULEN;
            IQTF = Qimag(k1,k2,DOF)/WAMITULEN;
            MQTF = Qmod(k1,k2,DOF)/WAMITULEN;
          elseif ge(DOF,4);
            RQTF = Qreal(k1,k2,DOF)/WAMITULEN^2;
            IQTF = Qimag(k1,k2,DOF)/WAMITULEN^2;
            MQTF = Qmod(k1,k2,DOF)/WAMITULEN^2;
          endif
          line = [2*pi/w1(k1) 2*pi/w1(k2) betaId(1) betaId(2)...
                  DOF MQTF atan2(IQTF,RQTF)*180/pi RQTF IQTF];
          fprintf(fid, '%e \t %e \t %e \t %e \t %d \t %e \t %e \t %e \t %e \n', line);
        endfor
    endfor
  endfor

  fclose(fid)
end

