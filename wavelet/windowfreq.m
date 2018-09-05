% window the frequency of the signal with cutoff and taper
% using plank taper window
% all values should be given in terms of *actual* frequency ratio
% the unmentioned ratio will be considered to be untapered
function z = windowfreq(w, leftcut, lefttaper, righttaper, rightcut)

if (leftcut+ lefttaper + righttaper + rightcut > 1)
	disp('total exceeds 1');
else
	n = length(w);
	N = floor(n/2);	% max frequency possible



	% get the actual number of time points based on ratio
	leftcutn = floor(N*leftcut);
	lefttapern = floor(N*lefttaper);
	righttapern = floor(N*righttaper);
	rightcutn = floor(N*rightcut);
	extran = N - (leftcutn + lefttapern + rightcutn + righttapern); 	% ratio that is not mentioned	;

	% prepare the left side
	lc = zeros(1,leftcutn); 
	lt= planktaper(2*lefttapern,0.5)(1:lefttapern);
       xt = (zeros(1,extran)+1);
      rt =  planktaper(2*righttapern,0.5)((righttapern+1):(2*righttapern));
     rc= zeros(1,rightcutn);

     halfcutoff = [lc lt xt rt rc];


     if (2*N == n)	% i.e even length
	cutoff = [halfcutoff flip(halfcutoff)];
else

	cutoff = [halfcutoff 1 flip(halfcutoff)];
end

     %% plot the cutoff function
     %%%% this helps plot then cutoff in anita freq
     plotfanita(ifft(cutoff));


	newf = fft(w).*cutoff;

	z =real(ifft(newf));
end





