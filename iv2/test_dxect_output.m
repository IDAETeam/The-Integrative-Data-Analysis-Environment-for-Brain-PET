 disp('To return IDAE home directory of the user:')
out           = dxetc4dfw('idx','Bhavin');
 disp(out)
  disp('To return segment structures of PET and MRI paths in data sercer:')
out           = dxetc4dfw('str',[]);
 disp(out)
 disp('To generate PET path for data server from PET path in image server')
 out        	= dxetc4dfw('pet_i2d','full/path/PET/in/image_server');
 disp(out)
 disp('To return PET path of image server from a PET file in data server')
 out        	= dxetc4dfw('pet_d2i','full/path/PET/in/data_server/file.ext');
 disp(out)
 disp('To retu >> rn variables for identification of PET source files')
 out        	= dxetc4dfw('get_pet',[]);
 disp(out)
disp(' To generate MRI path for data server from MRI path in image server')
 out        	= dxetc4dfw('mri_i2d','full/path/PET/in/image_server');
 disp(out)
 disp('To return variables for identification of MRI source files')
 out        	= dxetc4dfw('get_mri',[]);
 disp(out)
 disp('To convert a PET/MRI path between Windows and Linux systems')
 out           = dxetc4dfw('w2d','full/path/in/Linux_or_Windows/whatever');
 disp(out)
 disp('To return user-specific PET path in the data server:')
 out           = dxetc4dfw('ezr','full/path/PET/in/data_server/pet/folder');
 disp(out)
 disp('To return user-specific MRI path in the data server:')
 out         	= dxetc4dfw('ezr','full/path/MRI/in/data_server/pet/folder');
 disp(out)
 disp('To return VOI set registry file:')
 out           = dxetc4dfw('vst','Bhavin');
 disp(out)
 disp('To return Freesurfer-related paths:')
 out           = dxetc4dfw('fsd','Bhavin');
 disp(out)
 disp('To return local radioactivity unit')
 out          = dxetc4dfw('unit',[]);
 disp(out)  
 