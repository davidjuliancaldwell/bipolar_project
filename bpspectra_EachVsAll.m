function [M,Mbpdist,frx]=bpspectra_EachVsAll(d,sfx,frxrange,em,nchtocheck)
% for each window in d (samples x channels x windows), take channel c1
% minus channel c2 and computes spectrum on that bipolar signal, placing
% into a confusion matrix. Also does referential spectra (on the diagonal).
% OUTPUTS:
%   M
%   Mbpdist is a confusion matrix of the euclideandistances for every bipolar pair (channel x channel)
%   frx is the frrequency index for maz

%just getting frequency index
 [~,frx]=spectrogramjk_chronuxmtfft(zeros(1,size(d,1)),sfx,frxrange,[.5,1],0); 
 
 nwindtocheck=size(d,3);
 
 % Create matrix M:  channels X channels X frequencies X windows 
 M=nan(nchtocheck,nchtocheck,length(frx),nwindtocheck); 
 MMbpdist=nan(nchtocheck,nchtocheck); % will log the euclidean distance for each of these pairs
 tic
 for w = 1:nwindtocheck 
     parfor c1=1:nchtocheck 
         for c2=1:nchtocheck 
             trc=sq(d(:,c1,w))-sq(d(:,c2,w));
             if ~any(isnan(trc))
                 M(c1,c2,:,w)=spectrogramjk_chronuxmtfft(trc,sfx,frxrange,[.5,1],0);
             end  
         end
     end; disp([num2str(round(w/nwindtocheck*100,1)) '% done, elapsed time: ' num2str(toc)])
 end; clear c1 c2 trc
 % include referential signal power on the diagonal for storage/plotting
 for w = 1:nwindtocheck 
     for c1=1:nchtocheck
             trc=sq(d(:,c1,w));
             if ~any(isnan(trc))
                 M(c1,c1,:,w)=spectrogramjk_chronuxmtfft(trc,sfx,frxrange,[.5,1],0);
             end  
     end;
 end; 

 % euclidean distance for each pair
 for c1=1:nchtocheck; for c2=1:nchtocheck; Mbpdist(c1,c2)=distance3D(em(c1,:),em(c2,:)); end; end
    Mbpdist(Mbpdist==0)=nan; %remove spurious zero values, then
    for c1=1:nchtocheck; Mbpdist(c1,c1)=0; end %referential signal denoted as "0" distance for coding/indexing purposes below
 % remove mirror image in confusion matrix
 for c1=1:nchtocheck; for c2=c1+1:nchtocheck; M(c1,c2,:)=nan; Mbpdist(c1,c2)=nan; end; end

      
