function [xy,EL]=Q4_mesh_SEN(lgrede,hgrede,notch,nlelemenata,nhelemenata)
% Write a .txt file with single-edge notched geometry
% A plate of length lgrede and height hgrede with a verticaal notch on the
% bottom centre part of the beam is meshed. The mesh consists of a .txt
% file with node coordinates and element nodal connectivity. The notch is
% formed as double nodes with a one-node crack tip at the notch top.
% The file can then be read into MatLAB by using an appropriate reader.
% Procedure author: Nikola Lustig, mag.ing.aedif.


dl = lgrede / nlelemenata;
dh = hgrede / nhelemenata;
nC = ( nlelemenata + 1 ) * ( nhelemenata + 1 );
nEL = nlelemenata * nhelemenata;
nl_pukotine = nlelemenata/2 + 1;
nh_pukotine = notch / dh;


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
% Adding an extra column of nodes in the middle strip of the plate

ii = 1;
n = 0;
while ii <= nC + nhelemenata + 1
    if ii > length ( xy ) || xy ( 2 , ii ) >= notch
        break
    end;
    if round ( xy ( 1 , ii ) * 1000 ) / 1000 == round ( lgrede / 2 * 1000 ) / 1000
        xy = [ xy(:,1:ii),xy(:,ii:end) ];
        ii = ii + 1;
    end;
    ii = ii + 1;
end;

nC = length ( xy );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the node connectivity - different approaches used for the zone
% under the crack tip, the crack tip strip of elements and the zone above
% the crack tip.

EL = zeros ( 4 , nEL );

for jj = 1 : nhelemenata
    
    for ii = 1 : nlelemenata
        
        if jj < nh_pukotine
            
            if ii < nl_pukotine
                
                EL ( 1 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii;
                EL ( 2 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii + 1;
                EL ( 3 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii + 1;
                EL ( 4 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii;
                
            else
                
                EL ( 1 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii + 1;
                EL ( 2 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii + 2;
                EL ( 3 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii + 2;
                EL ( 4 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii + 1;
                
            end;
            
        elseif jj == nh_pukotine
            
            if ii < nl_pukotine
                
                EL ( 1 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii;
                EL ( 2 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii + 1;
                EL ( 3 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii + 1;
                EL ( 4 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii;
                
            else
                
                EL ( 1 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii + 1;
                EL ( 2 , ( jj - 1 ) * nlelemenata + ii ) = ( jj - 1 ) * ( nlelemenata + 2 ) + ii + 2;
                EL ( 3 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii + 1;
                EL ( 4 , ( jj - 1 ) * nlelemenata + ii ) = jj * ( nlelemenata + 2 ) + ii;
                
            end;
            
        else
            
            EL ( 1 , ( jj - 1 ) * nlelemenata + ii ) = nh_pukotine * ( nlelemenata + 2 ) + ( jj - nh_pukotine - 1 ) * ( nlelemenata + 1 ) + ii;
            EL ( 2 , ( jj - 1 ) * nlelemenata + ii ) = nh_pukotine * ( nlelemenata + 2 ) + ( jj - nh_pukotine - 1 ) * ( nlelemenata + 1 ) + ii + 1;
            EL ( 3 , ( jj - 1 ) * nlelemenata + ii ) = nh_pukotine * ( nlelemenata + 2 ) + ( jj - nh_pukotine ) * ( nlelemenata + 1 ) + ii + 1;
            EL ( 4 , ( jj - 1 ) * nlelemenata + ii ) = nh_pukotine * ( nlelemenata + 2 ) + ( jj - nh_pukotine ) * ( nlelemenata + 1 ) + ii;
            
        end;
        
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing the mesh out into a .txt file

textfilename = ['Q4_mesh_SEN_' num2str(nlelemenata) 'x' num2str(nhelemenata) '.txt'];
fid = fopen(textfilename, 'w');

if fid == -1
    error('Cannot open file for writing.');
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['% Q4 rectangular mesh for SEN geometry' nl '% Total number of nodes: ' num2str(nC) nl ...
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

% % Loop over elements for element types - Q4 --> EL_TIP = 1
% fwrite(fid, ['% Element type ' nl]);
% for ii = 1 : nEL
%
%     fwrite (fid, [num2str(ii) ' ' num2str(1) ' ;' nl]);
%
% end;



% Close the file.
fclose(fid);

end

