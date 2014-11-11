function Q4_Conf_force_to_VTK_force_only(nC,nEL,xy_popis,EL,SLIKA_pomaka,broj_koraka,Conf_force_popis,Scale_factor,Internal_node)


for k = 1 : broj_koraka
    
    trenutni_xy = xy_popis(:,:,k);

textfilename = ['vizualizacija\Cforce_Mueller1_' num2str(k) '.vtk'];
fid = fopen(textfilename, 'w');

% fid = fopen('uncoupled_Q4_INT4.vtk', 'w');

if fid == -1
    error('Cannot open file for writing.');
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['# vtk DataFile Version 2.0' nl  'Configurational force plotting' nl 'ASCII' nl ...
    'DATASET UNSTRUCTURED_GRID' nl 'FIELD FieldData 1 ' nl 'TIME 1 1 double' nl ' ' num2str(k-1) nl]);

% Loop over nodes for point coordinates
fwrite(fid, ['POINTS ' num2str(nC) ' float' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(trenutni_xy(1,ii)) ' ' num2str(trenutni_xy(2,ii)) ' 0.0' nl ]);
    
end;

% Loop over elements for nodal connectivity
fwrite(fid, ['CELLS ' num2str(nEL) ' ' num2str(nEL*5) nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(4) ' ' num2str(EL(2,ii)-1) ' ' num2str(EL(3,ii)-1) ' ' num2str(EL(4,ii)-1) ' ' num2str(EL(5,ii)-1) nl ]);
    
end;

% Loop over elements for cell types - 9 = VTK_QUAD
fwrite(fid, ['CELL_TYPES ' num2str(nEL) nl]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(9) nl]);
    
end;

% Loop over nodes for point data - displacements
fwrite(fid, [ 'POINT_DATA ' num2str(nC) nl ]);
fwrite(fid, ['VECTORS DISPLACEMENT float' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(SLIKA_pomaka(2*ii-1,k)) ' ' num2str(SLIKA_pomaka(2*ii,k)) ' 0.0' nl ]);
    
end;

% Loop over nodes for point data - configurational forces
% fwrite(fid, [ 'POINT_DATA ' num2str(nC) nl ]);
fwrite(fid, ['VECTORS CONFIGURATIONAL_FORCE float' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(Conf_force_popis(2*ii-1,k)) ' ' num2str(Conf_force_popis(2*ii,k)) ' 0.0' nl ]);
    
end;

% Loop over nodes for point data - scaled configurational forces for
% plotting glyphs
fwrite(fid, ['VECTORS CONFIGURATIONAL_FORCE_SCALED float' nl ]);
for ii = 1 : nC
    
    if sum(ii == Internal_node) == 0
    
        fwrite (fid, [num2str(Conf_force_popis(2*ii-1,k)*Scale_factor) ' ' num2str(Conf_force_popis(2*ii,k)*Scale_factor) ' 0.0' nl ]);
     
    else
        
        fwrite (fid, [num2str(Conf_force_popis(2*ii-1,k)) ' ' num2str(Conf_force_popis(2*ii,k)) ' 0.0' nl ]);
        
    end;
    
end;



% Close the file.
fclose(fid);


end;