function Q4_mesh_rectangular(lgrede,hgrede,nlelemenata,nhelemenata)
% Write a .txt file with regular rectangular geometry
% A plate of length lgrede and height hgrede with is meshed. The mesh 
% consists of a .txt file with node coordinates and element nodal 
% connectivity. The file can then be read into MatLAB by using an 
% appropriate reader. Procedure author: Nikola Lustig, mag.ing.aedif.


dl = lgrede / nlelemenata;
dh = hgrede / nhelemenata;
nC = ( nlelemenata + 1 ) * ( nhelemenata + 1 );
nEL = nlelemenata * nhelemenata;


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
% Printing the mesh out into a .txt file

textfilename = ['Q4_mesh_rectangular_' num2str(nlelemenata) 'x' num2str(nhelemenata) '.txt'];
fid = fopen(textfilename, 'w');

if fid == -1
    error('Cannot open file for writing.');
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['% Q4 rectangular mesh for regular geometry' nl '% Total number of nodes: ' num2str(nC) nl ...
     '% Total number of elements: ' num2str(nEL) nl ]);

% Loop over nodes for node coordinates
fwrite(fid, ['% Node coordinates ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(xy(1,ii)) ' ' num2str(xy(2,ii)) ' 0.0' ' ;' nl ]);
    
end;

% Loop over elements for nodal connectivity
fwrite(fid, ['% Element topology ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(EL(2,ii)) ' ' num2str(EL(3,ii)) ' ' num2str(EL(4,ii)) ' ' num2str(EL(5,ii)) ' ;' nl ]);
    
end;

% Loop over elements for element types - Q4 --> EL_TIP = 1
fwrite(fid, ['% Element type ' nl]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(1) ' ;' nl]);
    
end;


% Close the file.
fclose(fid);

