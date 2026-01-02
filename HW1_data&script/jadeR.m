function B =  jadeR(X,m)
%   B = jadeR(X, m) is an m*n matrix such that Y=B*X are separated sources
%    extracted from the n*T data matrix X.
%   If m is omitted,  B=jadeR(X)  is a square n*n matrix (as many sources as sensors)
%
% Blind separation of real signals with JADE.  Version 1.9.   August 2013
%
% Usage: 
%   * If X is an nxT data matrix (n sensors, T samples) then
%     B=jadeR(X) is a nxn separating matrix such that S=B*X is an nxT
%     matrix of estimated source signals.
%   * If B=jadeR(X,m), then B has size mxn so that only m sources are
%     extracted.  This is done by restricting the operation of jadeR
%     to the m first principal components. 
%   * Also, the rows of B are ordered such that the columns of pinv(B)
%     are in order of decreasing norm; this has the effect that the
%     `most energetically significant' components appear first in the
%     rows of S=B*X.
%
% Quick notes (more at the end of this file)
%
%  o this code is for REAL-valued signals.  An implementation of JADE
%    for both real and complex signals is also available from
%    http://perso.telecom-paristech.fr/~cardoso/guidesepsou.html
%
%  o This algorithm differs from the first released implementations of
%    JADE in that it has been optimized to deal more efficiently
%    1) with real signals (as opposed to complex)
%    2) with the case when the ICA model does not necessarily hold.
%
%  o There is a practical limit to the number of independent
%    components that can be extracted with this implementation.  Note
%    that the first step of JADE amounts to a PCA with dimensionality
%    reduction from n to m (which defaults to n).  In practice m
%    cannot be `very large' (more than 40, 50, 60... depending on
%    available memory and CPU time)
%
%  o See more notes, references and revision history at the end of
%    this file and more stuff on the WEB
%    http://perso.telecom-paristech.fr/~cardoso/guidesepsou.html
%
%  o This code is supposed to do a good job!  Please report any
%    problem to cardoso@sig.enst.fr


% Copyright (c) 2013, Jean-Francois Cardoso
% All rights reserved.
%
%
% BSD-like license.
% Redistribution and use in source and binary forms, with or without modification, 
% are permitted provided that the following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
%
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
%
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS 
% OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
% AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
% IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
% OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


verbose	= 1 ;	% Set to 0 for quiet operation

% Finding the number of sources
[n,T]	= size(X);
if nargin==1, m=n ; end; 	% Number of sources defaults to # of sensors
if m>n ,    fprintf('jade -> Do not ask more sources than sensors here!!!\n'), return,end
if verbose, fprintf('jade -> Looking for %d sources\n',m); end ;


% to do: add a warning about complex signals

% Mean removal
%=============
if verbose, fprintf('jade -> Removing the mean value\n'); end 
X	= X - mean(X')' * ones(1,T);


%%% whitening & projection onto signal subspace
%   ===========================================
if verbose, fprintf('jade -> Whitening the data\n'); end

[U,D]     = eig((X*X')/T) ; %% An eigen basis for the sample covariance matrix
[Ds,k]    = sort(diag(D)) ; %% Sort by increasing variance
PCs       = n:-1:n-m+1    ; %% The m most significant princip. comp. by decreasing variance

%% --- PCA  ----------------------------------------------------------
B         = U(:,k(PCs))'    ; % At this stage, B does the PCA on m components

%% --- Scaling  ------------------------------------------------------
scales    = sqrt(Ds(PCs)) ; % The scales of the principal components .
B         = diag(1./scales)*B  ; % Now, B does PCA followed by a rescaling = sphering


%% --- Sphering ------------------------------------------------------
X         = B*X;  %% We have done the easy part: B is a whitening matrix and X is white.

clear U D Ds k PCs scales ;

%%% NOTE: At this stage, X is a PCA analysis in m components of the real data, except that
%%% all its entries now have unit variance.  Any further rotation of X will preserve the
%%% property that X is a vector of uncorrelated components.  It remains to find the
%%% rotation matrix such that the entries of X are not only uncorrelated but also `as
%%% independent as possible'.  This independence is measured by correlations of order
%%% higher than 2.  We have defined such a measure of independence which
%%%   1) is a reasonable approximation of the mutual information
%%%   2) can be optimized by a `fast algorithm'
%%% This measure of independence also corresponds to the `diagonality' of a set of
%%% cumulant matrices.  The code below finds the `missing rotation ' as the matrix which
%%% best diagonalizes a particular set of cumulant matrices.

 
%%% Estimation of the cumulant matrices.
%   ====================================
if verbose, fprintf('jade -> Estimating cumulant matrices\n'); end

%% Reshaping of the data, hoping to speed up things a little bit...
X = X';

dimsymm 	= (m*(m+1))/2;	% Dim. of the space of real symm matrices
nbcm 		= dimsymm  ; 	% number of cumulant matrices
CM 		= zeros(m,m*nbcm);  % Storage for cumulant matrices
R 		= eye(m);  	%% 
Qij 		= zeros(m);	% Temp for a cum. matrix
Xim		= zeros(m,1);	% Temp
Xijm		= zeros(m,1);	% Temp
Uns		= ones(1,m);    % for convenience


%% I am using a symmetry trick to save storage.  I should write a short note one of these
%% days explaining what is going on here.
%%
Range     = 1:m ; % will index the columns of CM where to store the cumulant matrices.

for im = 1:m
  Xim = X(:,im) ;
  Xijm= Xim.*Xim ;
  %% Note to myself: the -R on next line can be removed: it does not affect
  %% the joint diagonalization criterion
  Qij           = ((Xijm(:,Uns).*X)' * X)/T - R - 2 * R(:,im)*R(:,im)' ;
  CM(:,Range)	= Qij ; 
  Range         = Range  + m ; 
  for jm = 1:im-1
    Xijm        = Xim.*X(:,jm) ;
    Qij         = sqrt(2) *(((Xijm(:,Uns).*X)' * X)/T - R(:,im)*R(:,jm)' - R(:,jm)*R(:,im)') ;
    CM(:,Range)	= Qij ;  
    Range       = Range  + m ;
  end ;
end;
%%%% Now we have nbcm = m(m+1)/2 cumulants matrices stored in a big m x m*nbcm array.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The inefficient code below does the same as above: computing the big CM cumulant matrix.
%% It is commented out but you can check that it produces the same result.
%% This is supposed to help interested people understand the (rather obscure) code above.
%% See section 4.2 of the Neural Comp paper referenced below.  It can be found at
%% "http://www.tsi.enst.fr/~cardoso/Papers.PS/neuralcomp_2ppf.ps",
%%
%% 
%%  
%%  if 1,
%%  
%%    %% Step one: we compute the sample cumulants
%%    Matcum = zeros(m,m,m,m) ;
%%    for i1=1:m,
%%      for i2=1:m,
%%        for i3=1:m,
%%  	for i4=1:m,
%%  	  Matcum(i1,i2,i3,i4) = mean( X(:,i1) .* X(:,i2) .* X(:,i3) .* X(:,i4) ) ...
%%  	      - R(i1,i2)*R(i3,i4) ...
%%  	      - R(i1,i3)*R(i2,i4) ...
%%  	      - R(i1,i4)*R(i2,i3) ;
%%  	end
%%        end
%%      end
%%    end
%%    
%%    %% Step 2; We compute a basis of the space of symmetric m*m matrices
%%    CMM = zeros(m, m, nbcm) ;  %% Holds the basis.   
%%    icm = 0                 ;  %% index to the elements of the basis
%%    vi          = zeros(m,1);  %% the ith basis vetor of R^m
%%    vj          = zeros(m,1);  %% the jth basis vetor of R^m
%%    Id          = eye  (m)  ;  %%  convenience
%%    for im=1:m,
%%      vi             = Id(:,im) ;
%%      icm            = icm + 1 ;
%%      CMM(:, :, icm) = vi*vi' ;
%%      for jm=1:im-1,
%%        vj             = Id(:,jm) ;
%%        icm            = icm + 1 ;
%%        CMM(:, :, icm) = sqrt(0.5) * (vi*vj'+vj*vi') ;
%%      end
%%    end
%%    %% Now CMM(:,:,i) is the ith element of an orthonormal basis for_ the space of m*m symmetric matrices
%%    
%%    %% Step 3.  We compute the image of each basis element by the cumulant tensor and store it back into CMM.
%%    mat = zeros(m) ; %% tmp
%%    for icm=1:nbcm
%%      mat = squeeze(CMM(:,:,icm)) ;
%%      for i1=1:m
%%        for i2=1:m
%%  	CMM(i1, i2, icm) = sum(sum(squeeze(Matcum(i1,i2,:,:))  .* mat )) ;
%%        end
%%      end
%%    end;
%%    %% This is doing something like  \sum_kl [ Cum(xi,xj,xk,xl) * mat_kl ] 
%%    
%%    %% Step 4.  Now, we can check that CMM and CM are equivalent
%%    Range = 1:m ;
%%    for icm=1:nbcm,
%%      M1    = squeeze( CMM(:,:,icm)) ;
%%      M2    = CM(:,Range) ; 
%%      Range = Range  + m ; 
%%      norm (M1-M2, 'fro' ) , %% This should be a numerical zero.
%%    end;
%%  
%%  end;  %%  End of the demo code for the computation of cumulant matrices
%%  
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Joint diagonalization of the cumulant matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Init
if 0, 	%% Init by diagonalizing a *single* cumulant matrix.  It seems to save
	%% some computation time `sometimes'.  Not clear if initialization is really worth
	%% it since Jacobi rotations are very efficient.  On the other hand, it does not
	%% cost much...

	if verbose, fprintf('jade -> Initialization of the diagonalization\n'); end
	[V,D]	= eig(CM(:,1:m)); % Selectng a particular cumulant matrix.
	for u=1:m:m*nbcm,         % Accordingly updating the cumulant set given the init
		CM(:,u:u+m-1) = CM(:,u:u+m-1)*V ; 
	end;
	CM	= V'*CM;

else,	%% The dont-try-to-be-smart init
	V	= eye(m) ; % la rotation initiale
end;

%% Computing the initial value of the contrast 
Diag    = zeros(m,1) ;
On      = 0 ;
Range   = 1:m ;
for im = 1:nbcm,
  Diag  = diag(CM(:,Range)) ;
  On    = On + sum(Diag.*Diag) ;
  Range = Range + m ;
end
Off = sum(sum(CM.*CM)) - On ;


seuil	= 1.0e-6 / sqrt(T) ; % A statistically scaled threshold on `small' angles
encore	= 1;
sweep	= 0; % sweep number
updates = 0; % Total number of rotations
upds    = 0; % Number of rotations in a given seep
g	= zeros(2,nbcm);
gg	= zeros(2,2);
G	= zeros(2,2);
c	= 0 ;
s 	= 0 ;
ton	= 0 ;
toff	= 0 ;
theta	= 0 ;
Gain    = 0 ;

%% Joint diagonalization proper
if verbose, fprintf('jade -> Contrast optimization by joint diagonalization\n'); end

while encore, encore=0;   

  if verbose, fprintf('jade -> Sweep #%3d',sweep); end
  sweep = sweep+1;
  upds  = 0 ; 
  Vkeep = V ;
  
  for p=1:m-1,
    for q=p+1:m,

      Ip = p:m:m*nbcm ;
      Iq = q:m:m*nbcm ;
      
      %%% computation of Givens angle
      g	    = [ CM(p,Ip)-CM(q,Iq) ; CM(p,Iq)+CM(q,Ip) ];
      gg    = g*g';
      ton   = gg(1,1)-gg(2,2); 
      toff  = gg(1,2)+gg(2,1);
      theta = 0.5*atan2( toff , ton+sqrt(ton*ton+toff*toff) );
      Gain  = (sqrt(ton*ton+toff*toff) - ton) / 4 ;
      
      %%% Givens update
      if abs(theta) > seuil,
%%      if Gain > 1.0e-3*On/m/m ,
	encore  = 1 ;
	upds    = upds    + 1;
	c	= cos(theta); 
	s	= sin(theta);
	G	= [ c -s ; s c ] ;
	
	pair 		= [p;q] ;
	V(:,pair) 	= V(:,pair)*G ;
	CM(pair,:)	= G' * CM(pair,:) ;
	CM(:,[Ip Iq]) 	= [ c*CM(:,Ip)+s*CM(:,Iq) -s*CM(:,Ip)+c*CM(:,Iq) ] ;
	

	On   = On  + Gain;
	Off  = Off - Gain;
	
	%% fprintf('jade -> %3d %3d %12.8f\n',p,q,Off/On);
      end%%of the if
    end%%of the loop on q
  end%%of the loop on p
  if verbose, fprintf(' completed in %d rotations\n',upds); end
  updates = updates + upds ;
  
end%%of the while loop
if verbose, fprintf('jade -> Total of %d Givens rotations\n',updates); end


%%% A separating matrix
%   ===================
B	= V'*B ;


%%% Permut the rows of the separating matrix B to get the most energetic components first.
%%% Here the **signals** are normalized to unit variance.  Therefore, the sort is
%%% according to the norm of the columns of A = pinv(B)

if verbose, fprintf('jade -> Sorting the components\n',updates); end
A           = pinv(B) ;
[Ds,keys]   = sort(sum(A.*A)) ;
B           = B(keys,:)       ;
B           = B(m:-1:1,:)     ; % Is this smart ?


% Signs are fixed by forcing the first column of B to have non-negative entries.

if verbose, fprintf('jade -> Fixing the signs\n',updates); end
b	= B(:,1) ;
signs	= sign(sign(b)+0.1) ; % just a trick to deal with sign=0
B	= diag(signs)*B ;



return ;