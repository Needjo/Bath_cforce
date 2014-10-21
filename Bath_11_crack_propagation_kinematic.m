clear;clc;close all;

lgrede = 10;
hgrede = 10;


nlelemenata = 4;
nhelemenata = 4;



textfilename = ['Mesh files\Q4_mesh_extended_' num2str(lgrede) '_' num2str(hgrede) '_EL_' num2str(nlelemenata) 'x' num2str(nhelemenata) '.txt'];
[ xy , EL , nC , nEL , nSEG , SEG , EL_SEG , EL_NEIGHBOUR , SEG_NEIGHBOUR , Node_occurence_EL , Node_occurence_SEG , Node_EL , Node_SEG ] = Q4_mesh_reader_extended ( textfilename );


for i = 1:nEL
    ELX ( 1 , i ) = EL ( 1 , i );
    ELX ( 2 , i ) = xy ( 1, EL ( 2 , i ) );
    ELX ( 3 , i ) = xy ( 1, EL ( 3 , i ) );
    ELX ( 5 , i ) = xy ( 1, EL ( 5 , i ) );
    ELX ( 4 , i ) = xy ( 1, EL ( 4 , i ) );
end;
for i = 1:nEL
    ELY ( 1 , i ) = EL ( 1 , i );
    ELY ( 2 , i ) = xy ( 2, EL ( 2 , i ) );
    ELY ( 3 , i ) = xy ( 2, EL ( 3 , i ) );
    ELY ( 5 , i ) = xy ( 2, EL ( 5 , i ) );
    ELY ( 4 , i ) = xy ( 2, EL ( 4 , i ) );
end;


ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);

ELXPLOT = [ELX(2:5,:);ELX(2,:)];
ELYPLOT = [ELY(2:5,:);ELY(2,:)];



figure

plot(ELXPLOT,ELYPLOT,'-k')





Crack_tip = 3;
Seg_tip = 23;

Element_tip_1 = SEG_NEIGHBOUR(1,Seg_tip);
Element_tip_2 = SEG_NEIGHBOUR(2,Seg_tip);

Seg_new = nSEG + 1;

% Adding the new segment on the last line and updating the number of
% segments
SEG = [ SEG , SEG(:,Seg_tip) ];
nSEG = nSEG + 1;

% Changing the neighbouring elements of the new and critical segment
SEG_NEIGHBOUR ( : , Seg_new ) = [ SEG_NEIGHBOUR(2,Seg_tip) 0 ]; 
SEG_NEIGHBOUR ( 2 , Seg_tip ) = 0;

% Changing the element segment connectivity - the right element is now
% connected to the new segment. The left element remains as it was.
for ii = 1 : 4
    
    if EL_SEG ( ii , Element_tip_2 ) == Seg_tip
        
        EL_SEG ( ii , Element_tip_2 ) = Seg_new;
        
    end;
    
end;

% Changing the neighbouring elements for the two elements at the crack tip
for ii = 1 : 4
    
    if EL_NEIGHBOUR ( ii , Element_tip_1 ) == Element_tip_2
        
        EL_NEIGHBOUR ( ii , Element_tip_1 ) = 0;
        
    end;
    
    if EL_NEIGHBOUR ( ii , Element_tip_2 ) == Element_tip_1
        
        EL_NEIGHBOUR ( ii , Element_tip_2 ) = 0;
        
    end;
    
end;


% Node doubling

xy = [ xy(:,1:Crack_tip) , xy(:,Crack_tip:nC) ];
nC = nC + 1;
New_node = Crack_tip+1;

EL(EL>Crack_tip)= EL(EL>Crack_tip) + 1;
SEG(SEG>Crack_tip)= SEG(SEG>Crack_tip) + 1;

Keeper = find ( sum( repmat(Seg_tip,4,nEL) ==  EL_SEG ) > 0 & sum ( repmat(Crack_tip,4,nEL) ==  EL(2:5,:) ) );

% logical_aid = false ( 4 , nEL );
logical_aid = repmat(Crack_tip,4,nEL) ==  EL(2:5,:) & [ true(4,Keeper-1) false(4,1) true(4,nEL-Keeper) ];
logical_aid = [ false(1,nEL) ; logical_aid ];

EL(logical_aid) = New_node;

NonKeeper = find ( sum( repmat(Seg_new,4,nEL) ==  EL_SEG ) > 0 & sum ( repmat(New_node,4,nEL) ==  EL(2:5,:) ) );


% logical_aid = false ( 2 , nSEG );
logical_aid = SEG_NEIGHBOUR == NonKeeper;
logical_aid = (SEG==Crack_tip) & (repmat(logical_aid(1,:) | logical_aid(2,:) ,2,1));

SEG(logical_aid) = New_node;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the elements which contain the node - regular quadrilateral
% elements
Node_EL = zeros ( nC , max ( Node_occurence_EL ) );

for ii = 1 : nC
    
    temp = find ( sum ( ii == EL(2:5,:) ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_EL ( ii , jj ) = temp ( jj );
        
    end;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the segments which contain the node - regular quadrilateral
% elements

Node_SEG = zeros ( nC , max ( Node_occurence_SEG ) );

for ii = 1 : nC
    
    temp = find ( sum ( ii == SEG ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_SEG ( ii , jj ) = temp ( jj );
        
    end;

end;





ELX = zeros ( 5 , nEL );
ELY = zeros ( 5 , nEL );
ELXPLOT = zeros ( 4 , nEL );
ELYPLOT = zeros ( 4 , nEL );

for i = 1:nEL
    ELX ( 1 , i ) = EL ( 1 , i );
    ELX ( 2 , i ) = xy ( 1, EL ( 2 , i ) );
    ELX ( 3 , i ) = xy ( 1, EL ( 3 , i ) );
    ELX ( 5 , i ) = xy ( 1, EL ( 5 , i ) );
    ELX ( 4 , i ) = xy ( 1, EL ( 4 , i ) );
end;
for i = 1:nEL
    ELY ( 1 , i ) = EL ( 1 , i );
    ELY ( 2 , i ) = xy ( 2, EL ( 2 , i ) );
    ELY ( 3 , i ) = xy ( 2, EL ( 3 , i ) );
    ELY ( 5 , i ) = xy ( 2, EL ( 5 , i ) );
    ELY ( 4 , i ) = xy ( 2, EL ( 4 , i ) );
end;


ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);

ELXPLOT = [ELX(2:5,:);ELX(2,:)];
ELYPLOT = [ELY(2:5,:);ELY(2,:)];
figure
plot(ELXPLOT,ELYPLOT,'-r')




