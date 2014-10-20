function Q4_mesh_rectangular_writer_extended(lgrede,hgrede,nlelemenata,nhelemenata)
% Write a .txt file with regular rectangular geometry
% A plate of length lgrede and height hgrede is meshed. The mesh 
% consists of a .txt file with node coordinates, element nodal 
% connectivity and several other data structures required for crack 
% propagation according to Miehe, C. and Gurses,E.: "A robust algorithm for
% configurational force driven brittle crack propagation with R-adaptive 
% mesh alignment" The file can then be read into MatLAB by using an 
% appropriate reader. Procedure author: Nikola Lustig, mag.ing.aedif.

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

dl = lgrede / nlelemenata;
dh = hgrede / nhelemenata;
nC = ( nlelemenata + 1 ) * ( nhelemenata + 1 );
nEL = nlelemenata * nhelemenata;
nSEG = ( nlelemenata + 1 ) * nhelemenata + nlelemenata * ( nhelemenata + 1 ); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regular rectangular node distribution

n = 1;
for ii = 1:(nhelemenata+1)
    for j = 1:(nlelemenata+1)
        if j==1
        xy ( 1 , n ) = 0;
        else
        xy ( 1 , n ) = xy ( 1, n - 1 ) + dl;
        end;
        xy ( 2 , n ) = ( ii - 1 ) * dh;
        n = n + 1;
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the node connectivity - regular quadrilateral elements

h = 0;
g = 1;
for ii = 1:nEL
    EL ( 1 , ii ) = ii;
    EL ( 2 , ii ) = h * ( nlelemenata + 1 ) + g;
    EL ( 3 , ii ) = h * ( nlelemenata + 1 ) + g + 1;
    EL ( 5 , ii ) = ( h + 1 ) * ( nlelemenata + 1 ) + g;
    EL ( 4 , ii ) = ( h + 1 ) * ( nlelemenata + 1 ) + g + 1;
     if rem( ii , nlelemenata ) == 0
    h = h + 1;
    g=1;
     else g=g+1;
    end;
end;
for ii = 1:nEL
    ELX ( 1 , ii ) = EL ( 1 , ii );
    ELX ( 2 , ii ) = xy ( 1, EL ( 2 , ii ) );
    ELX ( 3 , ii ) = xy ( 1, EL ( 3 , ii ) );
    ELX ( 5 , ii ) = xy ( 1, EL ( 5 , ii ) );
    ELX ( 4 , ii ) = xy ( 1, EL ( 4 , ii ) );
end;
for ii = 1:nEL
    ELY ( 1 , ii ) = EL ( 1 , ii );
    ELY ( 2 , ii ) = xy ( 2, EL ( 2 , ii ) );
    ELY ( 3 , ii ) = xy ( 2, EL ( 3 , ii ) );
    ELY ( 5 , ii ) = xy ( 2, EL ( 5 , ii ) );
    ELY ( 4 , ii ) = xy ( 2, EL ( 4 , ii ) );
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the segment node connectivity - regular quadrilateral elements
SEG = [ EL(2:3,:) , EL(5:-1:4,nEL-nlelemenata+1:nEL) ];

for ii = 1 : nhelemenata
    
    SEG = [ SEG , [EL(2,(ii-1)*nlelemenata+1);EL(5,(ii-1)*nlelemenata+1)] ,  EL(3:4,(ii-1)*nlelemenata+1:(ii-1)*nlelemenata+nlelemenata) ];
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the element segment connectivity - regular quadrilateral elements
EL_SEG = zeros ( 4 , nEL );

for ii = 1 : nEL
    
    kk = 1;
    
    for jj = 1 : nSEG
        
        if sum(ismember(SEG(:,jj),EL(2:5,ii))) == 2
            
            EL_SEG ( kk , ii ) = jj;
            kk = kk + 1;
            
        end;
        
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the element neighbouring elements - regular quadrilateral
% elements

EL_NEIGHBOUR = zeros ( 4 , nEL );

for ii = 1 : nEL
    
    kk = 1;
    
    for jj = 1 : nEL
        
        if sum(ismember(EL_SEG(:,ii),EL_SEG(:,jj))) == 1 && jj ~= ii
            
            EL_NEIGHBOUR ( kk , ii ) = jj;
            kk = kk + 1;
            
        end;
        
    end;
    
end;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the segment neighbouring elements - regular quadrilateral
% elements

SEG_NEIGHBOUR = zeros ( 2 , nSEG );

for ii = 1 : nSEG
    
    kk = 1;
    
    for jj = 1 : nEL
        
        if sum(ismember(ii,EL_SEG(:,jj))) == 1 
            
            SEG_NEIGHBOUR ( kk , ii ) = jj;
            kk = kk + 1;
            
        end;
        
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the node occurence in elements and segments - regular quadrilateral
% elements

Node_occurence_EL = zeros ( nC , 1 );

for ii = 1 : nC
    
    Node_occurence_EL ( ii , 1 ) = sum ( sum ( ii == EL ( 2:5 , : ) ) );
    
end;

Node_occurence_SEG = zeros ( nC , 1 );

for ii = 1 : nC
    
    Node_occurence_SEG ( ii , 1 ) = sum ( sum ( ii == SEG  ) );
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the segments which contain the node - regular quadrilateral
% elements

Node_SEG = zeros ( max ( Node_occurence_SEG ) );

for ii = 1 : nC
    
    temp = find ( sum ( ii == SEG ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_SEG ( ii , jj ) = temp ( jj );
        
    end;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the elements which contain the node - regular quadrilateral
% elements
Node_EL = zeros ( max ( Node_occurence_EL ) );

for ii = 1 : nC
    
    temp = find ( sum ( ii == EL(2:5,:) ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_EL ( ii , jj ) = temp ( jj );
        
    end;

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing the mesh out into a .txt file

textfilename = ['Q4_mesh_extended_' num2str(lgrede) '_' num2str(hgrede) '_EL_' num2str(nlelemenata) 'x' num2str(nhelemenata) '.txt'];
fid = fopen(textfilename, 'w');

if fid == -1
    error('Cannot open file for writing.');
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['% Q4 rectangular mesh for regular geometry' nl '% Total number of nodes: ' num2str(nC) nl ...
     '% Total number of elements: ' num2str(nEL) nl '% Total number of segments: ' num2str(nSEG) nl ]);

% Loop over nodes for node coordinates
fwrite(fid, ['% Node coordinates ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(xy(1,ii)) ' ' num2str(xy(2,ii)) ' 0.0' ' ;' nl ]);
    
end;

% Loop over elements for nodal connectivity
fwrite(fid, ['% Element node connectivity ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(EL(2,ii)) ' ' num2str(EL(3,ii)) ' ' num2str(EL(4,ii)) ' ' num2str(EL(5,ii)) ' ;' nl ]);
    
end;

% Loop over segments for nodal connectivity
fwrite(fid, ['% Segment node connectivity ' nl ]);
for ii = 1 : nSEG
    
    fwrite (fid, [num2str(ii) ' ' num2str(SEG(1,ii)) ' ' num2str(SEG(2,ii)) ' ;' nl ]);
    
end;

% Loop over elements for segment connectivity
fwrite(fid, ['% Element segment connectivity ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(EL_SEG(1,ii)) ' ' num2str(EL_SEG(2,ii)) ' ' num2str(EL_SEG(3,ii)) ' ' num2str(EL_SEG(4,ii)) ' ;' nl ]);
    
end;

% Loop over elements for element neighbours 
fwrite(fid, ['% Element neighbours ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(EL_NEIGHBOUR(1,ii)) ' ' num2str(EL_NEIGHBOUR(2,ii)) ' ' num2str(EL_NEIGHBOUR(3,ii)) ' ' num2str(EL_NEIGHBOUR(4,ii)) ' ;' nl ]);
    
end;

% Loop over segments for segment neighbouring elements
fwrite(fid, ['% Segment neighbours ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(SEG_NEIGHBOUR(1,ii)) ' ' num2str(SEG_NEIGHBOUR(2,ii)) ' ;' nl ]);
    
end;

% Loop over nodes for node occurence in elements
fwrite(fid, ['% Node occurence in elements ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_occurence_EL(ii,1)) ' ;' nl ]);
    
end;

% Loop over nodes for node occurence in segments
fwrite(fid, ['% Node occurence in segments ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_occurence_SEG(ii,1)) ' ;' nl ]);
    
end;

% Loop over nodes for elements containing the node
fwrite(fid, ['% Nodes appearing in elements ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_EL(ii,1)) ' ' num2str(Node_EL(ii,2)) ' ' num2str(Node_EL(ii,3)) ' ' num2str(Node_EL(ii,4)) ' ;' nl ]);
    
end;

% Loop over nodes for segments containing the node
fwrite(fid, ['% Nodes appearing in segments ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_SEG(ii,1)) ' ' num2str(Node_SEG(ii,2)) ' ' num2str(Node_SEG(ii,3)) ' ' num2str(Node_SEG(ii,4)) ' ;' nl ]);
    
end;



% Loop over elements for element types - Q4 --> EL_TIP = 1
fwrite(fid, ['% Element type ' nl]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(1) ' ;' nl]);
    
end;


% Close the file.
fclose(fid);

