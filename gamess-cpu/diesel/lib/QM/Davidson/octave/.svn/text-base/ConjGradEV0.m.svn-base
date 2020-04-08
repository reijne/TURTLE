function x = ConjGrad(A, b)

	D = inv(eye(size(A,1)).*A);
	x = D*b;

	normb = norm(b);	
	r = b - A*x;
	rho_1 = 0;
	if ( normb == 0 )
		normb = 1;
	endif
	
	i = 1;
	while ( i<1000 )
	
		z = D*r;
		rho = r'*z;
		
		if ( i==1 )
			p = z;
		else
			beta = rho/rho_1;
			p = z + beta*p;
		endif

		q = A*p;

		alpha = rho / (p'*q);
		
		x = x + alpha*p;
		r = r - alpha*q;

		resid = norm(r);
		resid/normb
		if ( resid/normb<=1e-5 )
			ConjGrad = x;
			break;
		endif
		rho_1 = rho;
		i++
	endwhile
endfunction

dim = 5
A=rand(dim);
A=A+A';
A= A - (5*rand(dim) + 5).*eye(dim);
A(dim,:) = 1*A(dim-1,:)+5*A(dim-2,:)+3*A(dim-3,:)+2*A(dim-4,:);
A
[V,l]=eig(A)
ev0=V(:,dim)

B=A(1:dim-1,1:dim-1)
eig(B)

b=rand(dim-1, 1);
b;
x=inv(B)*b;
xs=ConjGrad(B, b);
max(x-xs)
