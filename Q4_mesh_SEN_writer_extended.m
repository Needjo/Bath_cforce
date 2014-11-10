function Q4_mesh_SEN_writer_extended(lgrede,hgrede,x_notch,y_notch,notch_length,nlelemenata,nhelemenata)
% Write a .txt file with SEN rectangular geometry. Input from a regular
% mesh is needed which is then updated to include the notch.
% A plate of length lgrede and height hgrede is meshed with a notch on 
% position (x_notch,y_notch) of length notch_length. The mesh 
% consists of a .txt file with node coordinates, element nodal 
% connectivity and several other data structures required for crack 
% propagation according to Miehe, C. and Gurses,E.: "A robust algorithm for
% configurational force driven brittle crack propagation with R-adaptive 
% mesh alignment" The file can then be read into MatLAB by using an 
% appropriate reader. Procedure author: Nikola Lustig, mag.ing.aedif.
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
% File and location from where the rectangular .txt mesh is loaded
textfilename = ['Mesh files\Q4_mesh_extended_' num2str(lgrede) '_' num2str(hgrede) '_EL_' num2str(nlelemenata) 'x' num2str(nhelemenata) '.txt'];
[ xy , EL , nC , nEL , nSEG , SEG , EL_SEG , EL_NEIGHBOUR , SEG_NEIGHBOUR , Node_occurence_EL , Node_occurence_SEG , Node_EL , Node_SEG ] = Q4_mesh_reader_extended ( textfilename );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Depending on the location of the notch the direction of notch propagation
% is assumed. The notch is assumed to be positioned either vertically of
% horizontally. The tolerance is needed to find the node coordinates. The
% length of the notch must be represented by a finite number of elements.
dl = lgrede / nlelemenata;
dh = hgrede / nhelemenata;

tol = 10^-4;

if x_notch == 0 || x_notch == lgrede
    
    broj_koraka = notch_length / dl;
    
elseif y_notch == 0 || y_notch == hgrede
    
    broj_koraka = notch_length / dh;
    
else
    
    error ( 'The notch position can not be correctly identified!' )
    
end;


for kk = 1 : broj_koraka
% Loop over notch elements
% Depending on the location of the notch the crack direction is assumed.
% The notch always propagates on a straight line perpendicular to the edge.
% The Notch_print command classifies the notch direction as: 1-horizontal
% from left edge, 2-horizontal from right edge, 3-vertical from bottom
% edge, 4-vertical from top edge.
if x_notch == 0
    
    Crack_tip = find( abs(xy ( 1 , : ) - x_notch - (kk-1) * dl) <= tol & abs ( xy ( 2 , : ) - y_notch ) <= tol );
    Seg_tip = find(sum(repmat([Crack_tip;Crack_tip+1],1,nSEG)==SEG)==2);
    Notch_print = 1;
    
elseif x_notch == lgrede
    
    Crack_tip = find( abs(xy ( 1 , : ) - x_notch + (kk-1) * dl) <= tol & abs ( xy ( 2 , : ) - y_notch ) <= tol );
    Seg_tip = find(sum(repmat([Crack_tip-1;Crack_tip],1,nSEG)==SEG)==2);
    Notch_print = 2;
    
elseif y_notch == 0
    
    Crack_tip = find( abs(xy ( 1 , : ) - x_notch ) <= tol & abs ( xy ( 2 , : ) - y_notch - (kk-1) * dh ) <= tol );
    Seg_tip = find(sum(repmat([Crack_tip;Crack_tip+nlelemenata+1],1,nSEG)==SEG)==2);
    Notch_print = 3;
    
elseif y_notch == hgrede
    
    Crack_tip = find( abs(xy ( 1 , : ) - x_notch ) <= tol & abs ( xy ( 2 , : ) - y_notch + (kk-1) * dh ) <= tol );
    Seg_tip = find(sum(repmat([Crack_tip-nlelemenata-1;Crack_tip],1,nSEG)==SEG)==2);
    Notch_print = 4;
    
else
    
    error ('The notch position is not correctly defined!')
    
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting the crack propagation by invoking segment and node
% doubling and appropriate data structure updates. Crack_tip node
% is the node to be doubled while Seg_tip is the segment along
% which the crack propagation occurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the elements connected by the crack tip
Element_tip_1 = SEG_NEIGHBOUR(1,Seg_tip);
Element_tip_2 = SEG_NEIGHBOUR(2,Seg_tip);

% Defining the new segment, positioning it at the last place in the
% segment list), increasing the global number of segments
Seg_new = nSEG + 1;
SEG = [ SEG , SEG(:,Seg_tip) ];
nSEG = nSEG + 1;

% Changing the neighbouring elements of the new and critical segment
SEG_NEIGHBOUR ( : , Seg_new ) = [ SEG_NEIGHBOUR(2,Seg_tip) 0 ];
SEG_NEIGHBOUR ( 2 , Seg_tip ) = 0;

% Changing the element segment connectivity - the right element is now
% connected to the new segment. The left element remains as it was.
% Changing the neighbouring elements for the two elements at the crack tip
for ii = 1 : 4
    
    if EL_SEG ( ii , Element_tip_2 ) == Seg_tip
        
        EL_SEG ( ii , Element_tip_2 ) = Seg_new;
        
    end;
    
    if EL_NEIGHBOUR ( ii , Element_tip_1 ) == Element_tip_2
        
        EL_NEIGHBOUR ( ii , Element_tip_1 ) = 0;
        
    end;
    
    if EL_NEIGHBOUR ( ii , Element_tip_2 ) == Element_tip_1
        
        EL_NEIGHBOUR ( ii , Element_tip_2 ) = 0;
        
    end;
    
end;

% Node doubling, the new node is added next to the old node to
% retain a better global stiffness matrix structure. Increasing the
% total number of nodes.
xy = [ xy(:,1:Crack_tip) , xy(:,Crack_tip:nC) ];
nC = nC + 1;
New_node = Crack_tip+1;


% Renumbering of the nodes due to the addition of a new node in the
% mesh. Updating the boundary conditions and force vector
Cracker = EL>Crack_tip;
Cracker ( 1 , : ) = false;
EL(Cracker)= EL(Cracker) + 1;
SEG(SEG>Crack_tip)= SEG(SEG>Crack_tip) + 1;


% Finding the element which retains the old node Crack_tip
Keeper = find ( sum( repmat(Seg_tip,4,nEL) ==  EL_SEG ) > 0 & sum ( repmat(Crack_tip,4,nEL) ==  EL(2:5,:) ) );
Keeper = [ Keeper EL(1,sum(EL_NEIGHBOUR == Keeper) & sum(EL(2:5,:) == Crack_tip))];


% Logical operators for changing the element topology of the
% elements which contain the new node New_node
logical_aid = true ( 4 , nEL );
logical_aid ( : , Keeper ) = false;
logical_aid = repmat(Crack_tip,4,nEL) ==  EL(2:5,:) & logical_aid;
logical_aid = [ false(1,nEL) ; logical_aid ];

EL(logical_aid) = New_node;

% Finding the segments which get the new node New_node
NonKeeper = find ( sum( repmat(Seg_new,4,nEL) ==  EL_SEG ) > 0 & sum ( repmat(New_node,4,nEL) ==  EL(2:5,:) ) );

% Logical operators for changing the segment topology of the
% segments which contain the new node New_node
logical_aid = SEG_NEIGHBOUR == NonKeeper;
logical_aid = (SEG==Crack_tip) & (repmat(logical_aid(1,:) | logical_aid(2,:) ,2,1));

SEG(logical_aid) = New_node;

for ii = 1 : 2
    
    if SEG_NEIGHBOUR ( ii , NonKeeper ) ~= Element_tip_1 & SEG_NEIGHBOUR ( ii , NonKeeper ) ~= Element_tip_2
        
        NonKeeper_element = SEG_NEIGHBOUR ( ii , NonKeeper );
        
    end;
    
end;

for ii = 1 : 4
    
    if NonKeeper_element == 0
        
        break
        
    end;
    
    for jj = 1 : 2
        
        if SEG(jj,EL_SEG ( ii , NonKeeper_element ) ) == Crack_tip & EL_SEG ( ii , NonKeeper_element ) ~= NonKeeper
            
            SEG ( jj,EL_SEG ( ii , NonKeeper_element ) ) = New_node;
            
        end;
        
    end;
    
end;

% Defining the node occurence in elements and segments - regular quadrilateral
% elements

Node_occurence_EL = zeros ( nC , 1 );
Node_occurence_SEG = zeros ( nC , 1 );

for ii = 1 : nC
    
    Node_occurence_EL ( ii , 1 ) = sum ( sum ( ii == EL ( 2:5 , : ) ) );
    Node_occurence_SEG ( ii , 1 ) = sum ( sum ( ii == SEG  ) );
    
end;

% Defining the elements and segments which contain the node
Node_EL = zeros ( nC , max ( Node_occurence_EL ) );
Node_SEG = zeros ( nC , max ( Node_occurence_SEG ) );

for ii = 1 : nC
    
    temp = find ( sum ( ii == EL(2:5,:) ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_EL ( ii , jj ) = temp ( jj );
        
    end;
    
    temp = find ( sum ( ii == SEG ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_SEG ( ii , jj ) = temp ( jj );
        
    end;
    
end;

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing the SEN mesh out into a .txt file

textfilename = ['Mesh files\Q4_mesh_SEN_extended_' num2str(lgrede) '_' num2str(hgrede) '_' num2str(notch_length) '_EL_' num2str(nlelemenata) 'x' num2str(nhelemenata) '_DIR_' num2str(Notch_print) '.txt'];
fid = fopen(textfilename, 'w');

if fid == -1
    error('Cannot open file for writing.');
end

% New line.
nl = sprintf('\n'); % Stupid matlab doesn't interpret \n normally.

% Write the file header.
fwrite(fid, ['% Q4 SEN mesh for regular geometry' nl '% Total number of nodes: ' num2str(nC) nl ...
     '% Total number of elements: ' num2str(nEL) nl '% Total number of segments: ' num2str(nSEG) nl ]);

% Loop over nodes for node coordinates
fwrite(fid, ['% Node coordinates - xy' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(xy(1,ii)) ' ' num2str(xy(2,ii)) ' 0.0' ' ;' nl ]);
    
end;

% Loop over elements for nodal connectivity
fwrite(fid, ['% Element node connectivity - EL ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(EL(2,ii)) ' ' num2str(EL(3,ii)) ' ' num2str(EL(4,ii)) ' ' num2str(EL(5,ii)) ' ;' nl ]);
    
end;

% Loop over segments for nodal connectivity
fwrite(fid, ['% Segment node connectivity - SEG ' nl ]);
for ii = 1 : nSEG
    
    fwrite (fid, [num2str(ii) ' ' num2str(SEG(1,ii)) ' ' num2str(SEG(2,ii)) ' ;' nl ]);
    
end;

% Loop over elements for segment connectivity
fwrite(fid, ['% Element segment connectivity - EL_SEG ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(EL_SEG(1,ii)) ' ' num2str(EL_SEG(2,ii)) ' ' num2str(EL_SEG(3,ii)) ' ' num2str(EL_SEG(4,ii)) ' ;' nl ]);
    
end;

% Loop over elements for element neighbours 
fwrite(fid, ['% Element neighbours - EL_NEIGHBOUR ' nl ]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(EL_NEIGHBOUR(1,ii)) ' ' num2str(EL_NEIGHBOUR(2,ii)) ' ' num2str(EL_NEIGHBOUR(3,ii)) ' ' num2str(EL_NEIGHBOUR(4,ii)) ' ;' nl ]);
    
end;

% Loop over segments for segment neighbouring elements
fwrite(fid, ['% Segment neighbours - SEG_NEIGHBOUR ' nl ]);
for ii = 1 : nSEG
    
    fwrite (fid, [num2str(ii) ' ' num2str(SEG_NEIGHBOUR(1,ii)) ' ' num2str(SEG_NEIGHBOUR(2,ii)) ' ;' nl ]);
    
end;

% Loop over nodes for node occurence in elements
fwrite(fid, ['% Node occurence in elements - Node_occurence_EL ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_occurence_EL(ii,1)) ' ;' nl ]);
    
end;

% Loop over nodes for node occurence in segments
fwrite(fid, ['% Node occurence in segments - Node_occurence_SEG ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_occurence_SEG(ii,1)) ' ;' nl ]);
    
end;

% Loop over nodes for elements containing the node
fwrite(fid, ['% Nodes appearing in elements - Node_EL ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_EL(ii,1)) ' ' num2str(Node_EL(ii,2)) ' ' num2str(Node_EL(ii,3)) ' ' num2str(Node_EL(ii,4)) ' ;' nl ]);
    
end;

% Loop over nodes for segments containing the node
fwrite(fid, ['% Nodes appearing in segments - Node_SEG ' nl ]);
for ii = 1 : nC
    
    fwrite (fid, [num2str(ii) ' ' num2str(Node_SEG(ii,1)) ' ' num2str(Node_SEG(ii,2)) ' ' num2str(Node_SEG(ii,3)) ' ' num2str(Node_SEG(ii,4)) ' ;' nl ]);
    
end;



% Loop over elements for element types - Q4 --> EL_TIP = 1
fwrite(fid, ['% Element type - EL_TIP ' nl]);
for ii = 1 : nEL
    
    fwrite (fid, [num2str(ii) ' ' num2str(1) ' ;' nl]);
    
end;


% Close the file.
fclose(fid);

