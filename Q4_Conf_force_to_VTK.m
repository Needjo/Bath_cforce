function Q4_Conf_force_to_VTK(nC,nEL,xy_popis,EL,SLIKA_pomaka,broj_koraka,Stress_point,Conf_popis)


for k = 1 : broj_koraka
    
    trenutni_xy = xy_popis(:,:,k);

textfilename = ['vizualizacija\Cforce_Mueller2_' num2str(k) '.vtk'];
fid = fopen(textfilename, 'w');

% fid = fopen('uncoupled_Q4_INT4.vtk', 'w');

if fid == -1
    error('Cannot open file for writing.');
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['# vtk DataFile Version 2.0' nl  'Configurational force printout' nl 'ASCII' nl ...
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
fwrite(fid, ['VECTORS CONFIGURATIONAL_FORCE float' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(Conf_popis(2*ii-1,k)) ' ' num2str(Conf_popis(2*ii,k)) ' 0.0' nl ]);
    
end;




% Loop over nodes for point data - stress at nodes

fwrite(fid, ['SCALARS STRESS_xx float' nl ]);
fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(Stress_point(ii,1,k)) nl ]);
    
end;

fwrite(fid, ['SCALARS STRESS_yy float' nl ]);
fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(Stress_point(ii,2,k)) nl ]);
    
end;

fwrite(fid, ['SCALARS STRESS_xy float' nl ]);
fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(Stress_point(ii,3,k)) nl ]);
    
end;

fwrite(fid, ['SCALARS STRESS_11 float' nl ]);
fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(Stress_point(ii,4,k)) nl ]);
    
end;

fwrite(fid, ['SCALARS STRESS_22 float' nl ]);
fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
for ii = 1 : nC
    
     fwrite (fid, [num2str(Stress_point(ii,5,k)) nl ]);
    
end;

% % % Loop over elements for stress quantities - Sigma 11 at centroid
% % fwrite(fid, [ 'CELL_DATA  ', num2str(nEL), nl ]);
% % fwrite(fid, [ 'SCALARS STRESS_11_centre float 1' nl ]);
% % fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
% % for ii = 1 : nEL
% %     
% %     fwrite (fid, [num2str(Naprezanja_element(ii,3)) nl]);
% %     
% % end;
% % 
% % 
% % % Loop over elements for stress quantities - Sigma 22 at centroid
% % fwrite(fid, [ 'SCALARS STRESS_22_centre float 1' nl ]);
% % fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
% % for ii = 1 : nEL
% %     
% %     fwrite (fid, [num2str(Naprezanja_element(ii,4)) nl]);
% %     
% % end;

% % % % fwrite(fid, [ 'CELL_DATA 1' nl ]);
% % % % fwrite(fid, [ 'SCALARS POMAK_GLOBALNI float 1' nl ]);
% % % % fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
% % % % fwrite(fid, [ num2str(SLIKA_Stanja_vrh(k,2)) nl]);
% % % % 
% % % % fwrite(fid, [ 'SCALARS SILA_GLOBALNO float 1' nl ]);
% % % % fwrite(fid, [ 'LOOKUP_TABLE default ' nl ]);
% % % % fwrite(fid, [ num2str(SLIKA_Stanja_vrh(k,1)) nl]);


% Close the file.
fclose(fid);


end;