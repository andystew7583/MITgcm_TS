%%% 
%%% backupAverages.m
%%%
%%% Creates a .mat backup file of all of the time-averaged model fields,
%%% for convenient storage.
%%%

loadexp;
avg_t;
sdir = './backups';
mkdir(sdir);
sfname = fullfile(sdir,[expname,'_backup.mat']);
if (exist(sfname,'file'))
  delete(sfname);
end
save(sfname,'expname','tmin','tmax');
for m=1:length(diag_fields)
  if (strcmp(diag_fields(m),'UVEL'));
    save(sfname,'uu','-append');
  end
  if (strcmp(diag_fields(m),'VVEL'));
    save(sfname,'vv','-append');
  end
  if (strcmp(diag_fields(m),'WVEL'));
    save(sfname,'ww','-append');
  end
  if (strcmp(diag_fields(m),'THETA'));
    save(sfname,'tt','-append');
  end
  if (strcmp(diag_fields(m),'SALT'));
    save(sfname,'ss','-append');
  end
  if (strcmp(diag_fields(m),'UVELTH'));
    save(sfname,'ut','-append');
  end
  if (strcmp(diag_fields(m),'VVELTH'));
    save(sfname,'vt','-append');
  end
  if (strcmp(diag_fields(m),'WVELTH'));
    save(sfname,'wt','-append');
  end
  if (strcmp(diag_fields(m),'UVELSLT'));
    save(sfname,'us','-append');
  end
  if (strcmp(diag_fields(m),'VVELSLT'));
    save(sfname,'vs','-append');
  end
  if (strcmp(diag_fields(m),'WVELSLT'));
    save(sfname,'ws','-append');
  end
  if (strcmp(diag_fields(m),'THETASQ'));
    save(sfname,'tsq','-append');
  end
  if (strcmp(diag_fields(m),'SALTSQ'));
    save(sfname,'ssq','-append');
  end
  if (strcmp(diag_fields(m),'THSLT'));
    save(sfname,'ts','-append');
  end
  if (strcmp(diag_fields(m),'UVELSQ'));
    save(sfname,'usq','-append');
  end
  if (strcmp(diag_fields(m),'VVELSQ'));
    save(sfname,'vsq','-append');
  end
  if (strcmp(diag_fields(m),'WVELSQ'));
    save(sfname,'wsq','-append');
  end
end