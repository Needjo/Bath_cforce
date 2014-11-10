function [ xy , EL , nC , nEL , nSEG , SEG , EL_SEG , EL_NEIGHBOUR , SEG_NEIGHBOUR , Node_occurence_EL , Node_occurence_SEG , Node_EL , Node_SEG ] = Q4_mesh_reader_extended ( textfilename )
% Reader function for reading .txt mesh files
% The mesh files need to be created with the appropriate functions (e.g.
% Q4_mesh_rectangular_writer_extended or Q4_mesh_SEN_writer_extended)
% Procedure author: Nikola Lustig, mag.ing.aedif.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nEL = numer of elements in the mesh
% nC = number of nodes in the mesh
% nSEG = number of segments in the whole mesh
% xy = node ccordinates
% EL(2:5,:) = nodal connectivity for each element
% SEG = nodal connectivity for each segment (segment equivalent to EL)
% EL_SEG = element connectivity for segments (segment equivalent to EL)
% EL_NEIGHBOUR = neighbouring elements of element EL
% SEG_NEIGHBOUR = neighbouring elements of segment SEG
% Node_occurence_EL = number of occurences for each node
% Node_occurence_SEG = number of occurences for each segment
% Node_EL = elements having the node 
% Node_SEG = segments having the node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Opening the file and reading the header
fid = fopen(textfilename);
nC = fscanf(fid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %d', 1);

nEL = fscanf(fid, '%*s %*s %*s %*s %*s %d', 1);

nSEG = fscanf(fid, '%*s %*s %*s %*s %*s %d', 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xy = zeros ( 2 , nC );

xy ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %f' , 1);
xy ( 2 , 1 ) = fscanf ( fid , '%f' , 1 );

for ii = 1 : nC - 1 
    
    xy ( 1 , ii + 1 ) = fscanf ( fid , '%*s %*s %*s %f' , 1 );
    xy ( 2 , ii + 1 ) = fscanf ( fid , '%f' , 1 );
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EL = zeros ( 5 , nEL );

EL ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %*s %d' , 1);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SEG = zeros ( 2 , nSEG );

SEG ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %*s %d' , 1);
SEG ( 2 , 1 ) = fscanf ( fid , '%d' , 1);

for ii = 1 : nSEG - 1
    
    SEG ( 1 , ii + 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    SEG ( 2 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EL_SEG = zeros ( 4 , nEL );

EL_SEG ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %*s %d' , 1);
EL_SEG ( 2 , 1 ) = fscanf ( fid , '%d' , 1);
EL_SEG ( 3 , 1 ) = fscanf ( fid , '%d' , 1);
EL_SEG ( 4 , 1 ) = fscanf ( fid , '%d' , 1);

for ii = 1 : nEL - 1
    
    EL_SEG ( 1 , ii + 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    EL_SEG ( 2 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    EL_SEG ( 3 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    EL_SEG ( 4 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EL_NEIGHBOUR = zeros ( 4 , nEL );

EL_NEIGHBOUR ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %d' , 1);
EL_NEIGHBOUR ( 2 , 1 ) = fscanf ( fid , '%d' , 1);
EL_NEIGHBOUR ( 3 , 1 ) = fscanf ( fid , '%d' , 1);
EL_NEIGHBOUR ( 4 , 1 ) = fscanf ( fid , '%d' , 1);

for ii = 1 : nEL - 1
    
    EL_NEIGHBOUR ( 1 , ii + 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    EL_NEIGHBOUR ( 2 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    EL_NEIGHBOUR ( 3 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    EL_NEIGHBOUR ( 4 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SEG_NEIGHBOUR = zeros ( 2 , nSEG );

SEG_NEIGHBOUR ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %d' , 1);
SEG_NEIGHBOUR ( 2 , 1 ) = fscanf ( fid , '%d' , 1);

for ii = 1 : nSEG - 1
    
    SEG_NEIGHBOUR ( 1 , ii + 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    SEG_NEIGHBOUR ( 2 , ii + 1 ) = fscanf ( fid , '%d' , 1 );
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Node_occurence_EL = zeros ( nC , 1 );

Node_occurence_EL ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %*s %*s %d' , 1);

for ii = 1 : nC - 1
    
    Node_occurence_EL ( ii + 1 , 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Node_occurence_SEG = zeros ( nC , 1 );

Node_occurence_SEG ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %*s %*s %d' , 1);

for ii = 1 : nC - 1
    
    Node_occurence_SEG ( ii + 1 , 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Node_EL = zeros ( nC , max ( Node_occurence_EL ) );

Node_EL ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %*s %*s %d' , 1);
Node_EL ( 1 , 2 ) = fscanf ( fid , '%d' , 1);
Node_EL ( 1 , 3 ) = fscanf ( fid , '%d' , 1);
Node_EL ( 1 , 4 ) = fscanf ( fid , '%d' , 1);

for ii = 1 : nC - 1
    
    Node_EL ( ii + 1 , 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    Node_EL ( ii + 1 , 2 ) = fscanf ( fid , '%d' , 1 );
    Node_EL ( ii + 1 , 3 ) = fscanf ( fid , '%d' , 1 );
    Node_EL ( ii + 1 , 4 ) = fscanf ( fid , '%d' , 1 );
 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Node_SEG = zeros (nC ,  max ( Node_occurence_SEG ) );

Node_SEG ( 1 , 1 ) = fscanf ( fid , '%*s %*s %*s %*s %*s %*s %*s %*s %*s %d' , 1);
Node_SEG ( 1 , 2 ) = fscanf ( fid , '%d' , 1);
Node_SEG ( 1 , 3 ) = fscanf ( fid , '%d' , 1);
Node_SEG ( 1 , 4 ) = fscanf ( fid , '%d' , 1);

for ii = 1 : nC - 1
    
    Node_SEG ( ii + 1 , 1 ) = fscanf ( fid , '%*s %*s %d' , 1 );
    Node_SEG ( ii + 1 , 2 ) = fscanf ( fid , '%d' , 1 );
    Node_SEG ( ii + 1 , 3 ) = fscanf ( fid , '%d' , 1 );
    Node_SEG ( ii + 1 , 4 ) = fscanf ( fid , '%d' , 1 );
 
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fid);
    
    