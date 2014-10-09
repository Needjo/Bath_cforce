
% Working reader

fid = fopen('Q4_mesh_SEN_200x20.txt');
nC = fscanf(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d', 1);

nEL = fscanf(fid, '%*s %*s %*s %*s %*s %d', 1);

xy = zeros ( 2 , nC );

xy ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %f' , 1);
xy ( 2 , 1 ) = fscanf ( fid , '%f' , 1 );

for ii = 1 : nC - 1 
    
    xy ( 1 , ii + 1 ) = fscanf ( fid , '%*s %*s %*s %f' , 1 );
    xy ( 2 , ii + 1 ) = fscanf ( fid , '%f' , 1 );
    
end;

EL = zeros ( 5 , nEL );

EL ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %d' , 1);
EL ( 2 , 1 ) = fscanf ( fid , '%d' , 1);
EL ( 3 , 1 ) = fscanf ( fid , '%d' , 1);
EL ( 4 , 1 ) = fscanf ( fid , '%d' , 1);
EL ( 5 , 1 ) = fscanf ( fid , '%d' , 1);

for ii = 1 : nEL - 1
    
    EL ( 1 , ii + 1 ) = fscanf ( fid , '%*s %d' , 1 );
    EL ( 2 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    EL ( 3 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    EL ( 4 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    EL ( 5 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    
end;

fclose(fid);
    
    