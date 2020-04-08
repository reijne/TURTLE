function x = LGSIter(A, b)
	D = inv(eye(size(A,1)).*A);
	
	delta = [1];
	x = D*b;
	
	while ( abs(max(delta))>1e-10 )
		r = b - A*x;
		delta = D*r;
		x = x + delta;
		max(delta)
	endwhile
	
	LGSIter = x;
endfunction

dim = 10
A=rand(dim);
A=A+A';
A= A - (5*rand(dim) + 100).*eye(dim);
A
b=rand(dim, 1);
b
inv(A)*b
LGSIter(A, b)
