basename = 'movie_output/frame';
n = 1;
nimages = 10;
ext = '.png';
% M = moviein(nimages);


aviobj = avifile('~/Desktop/test.avi','fps',10,'quality',100);

for i = 1:1:nimages
  
  imgname = strcat(basename,num2str(n),ext);
  
  if (exist(imgname,'file') ~= 2)
    'Ran out of images!'
    break;
  else
    imgdata = imread(imgname);    
    aviobj = addframe(aviobj,imgdata);    
    n = n+1;
  end
  
end  

aviobj = close(aviobj);
