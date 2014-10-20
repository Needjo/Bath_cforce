function [ xy , EL , nC , nEL , nSEG ] = Q4_mesh_reader_extended ( textfilename )
% Reader function for reading .txt mesh files
% The mesh files need to be created with the appropriate functions (e.g.
% Q4_mesh_rectangular_writer_extended or Q4_mesh_SEN_writer_extended
% Procedure author: Nikola Lustig, mag.ing.aedif.

fid = fopen(textfilename);
nC = fscanf(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d', 1);

nEL = fscanf(fid, '%*s %*s %*s %*s %*s %d', 1);

nSEG = fscanf(fid, '%*s %*s %*s %*s %*s %d', 1);

xy = zeros ( 2 , nC );

xy ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %f' , 1);
xy ( 2 , 1 ) = fscanf ( fid , '%f' , 1 );

for ii = 1 : nC - 1 
    
    xy ( 1 , ii + 1 ) = fscanf ( fid , '%*s %*s %*s %f' , 1 );
    xy ( 2 , ii + 1 ) = fscanf ( fid , '%f' , 1 );
    
end;

EL = zeros ( 5 , nEL );

EL ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %d' , 1);
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
    
    