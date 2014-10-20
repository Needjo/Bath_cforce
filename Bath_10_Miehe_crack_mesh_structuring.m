% Crack propagation model - data structures and mesh adjustment for
% the critical node and critical segment - Bath 20.10.2014

clear; 
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters

lgrede = 10;
hgrede = 10;
tgrede = 1;

nlelemenata = 2;
nhelemenata = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finite element mesh creation

xy = pravokutni_mesh4kutni( lgrede , hgrede , nlelemenata , nhelemenata );
[ELX , ELY , EL ] = formiranje_4kutnihelemenata( xy, nlelemenata, nhelemenata );

ELX1 = ELX(2:5,:);
ELY1 = ELY(2:5,:);

ELXPLOT = [ELX(2:5,:);ELX(2,:)];
ELYPLOT = [ELY(2:5,:);ELY(2,:)];


nC = length ( xy );
nEL = nlelemenata * nhelemenata;

plot(ELXPLOT,ELYPLOT,'k-')

% nEL = numer of elements in the mesh
% nC = number of nodes in the mesh
% xy = node ccordinates
% EL(2:5,:) = nodal connectivity for each element
% nSEG = number of segments in the whole mesh
% SEG = nodal connectivity for each segment (segment equivalent to EL)
% EL_SEG = element connectivity for segments (segment equivalent to EL)
% EL_NEIGHBOUR = neighbouring elements of element EL
% SEG_NEIGHBOUR = neighbouring elements of segment SEG
% Node_occurence_EL = number of occurences for each node
% Node_occurence_SEG = number of occurences for each segment
% Node_EL = elements having the node 
% Node_SEG = segments having the node

nSEG = ( nlelemenata + 1 ) * nhelemenata + nlelemenata * ( nhelemenata + 1 ); 

SEG = [ EL(2:3,:) , EL(5:-1:4,nEL-nlelemenata+1:nEL) ];

for ii = 1 : nhelemenata
    
    SEG = [ SEG , [EL(2,(ii-1)*nlelemenata+1);EL(5,(ii-1)*nlelemenata+1)] ,  EL(3:4,(ii-1)*nlelemenata+1:(ii-1)*nlelemenata+nlelemenata) ];
    
end;

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

Node_occurence_EL = zeros ( nC , 1 );

for ii = 1 : nC
    
    Node_occurence_EL ( ii , 1 ) = sum ( sum ( ii == EL ( 2:5 , : ) ) );
    
end;

Node_occurence_SEG = zeros ( nC , 1 );

for ii = 1 : nC
    
    Node_occurence_SEG ( ii , 1 ) = sum ( sum ( ii == SEG  ) );
    
end;

Node_SEG = zeros ( max ( Node_occurence_SEG ) );

for ii = 1 : nC
    
    temp = find ( sum ( ii == SEG ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_SEG ( ii , jj ) = temp ( jj );
        
    end;

end;

Node_EL = zeros ( max ( Node_occurence_EL ) );

for ii = 1 : nC
    
    temp = find ( sum ( ii == EL(2:5,:) ) > 0  );
    
    for jj = 1 : length ( temp )
        
        Node_EL ( ii , jj ) = temp ( jj );
        
    end;

end;



          


