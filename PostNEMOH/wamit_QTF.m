function [Qmod,Qreal,Qimag] = wamit_QTF(projdir,Nw,dof_data,betaId,Nbetadat,qtftype,...
                              SwitchBiDir,WAMITULEN,fpath,runname)

  % QTF file name (depending on qtftype)
  QTFfname = [projdir,filesep,'results',filesep,'QTF',filesep,'OUT_QTF',qtftype,'_N.dat'];

  for DOF = 1:6
    [Qreal(:,:,DOF),Qimag(:,:,DOF),Qmod(:,:,DOF),w1]=fun_readNEMOHQTFOutput_N(QTFfname,...
                           Nw,dof_data,DOF,betaId,Nbetadat,qtftype,SwitchBiDir);
  endfor

end

