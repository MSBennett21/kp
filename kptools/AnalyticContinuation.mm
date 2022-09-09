
# This procedure computes the analytic continuation of 
# a set of function values on the Riemann surface, 
# along a given path. The output is a map from the initial
# points to the end points, given in the form of a table. 
AnalyticContinuation:=proc(curve,x,y,path,t,DI)
	local tstep,i,fibre,t1,dig,t_old,R,CL;
	global very_close;
	options remember, 
	    `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	# Not including the option remember results in going
	# through these slow calculations more than necessary.
	
	userinfo(6,'algcurves',`following path`,path);
	fibre:=[fsolve(subs(x= subs(t=0,path) ,curve),y,complex)];
	R:=[0,fibre];
	dig:=Digits;
	CL:=round(sqrt(max(1,dig-8)));
	if assigned(_Env_algcurves_CL) then CL:=_Env_algcurves_CL fi;
	tstep:=Initial_tstep / CL^2;
	t_old:=0;
	while t_old<1 do
		t1:=t_old + tstep;
		userinfo(9,'algcurves',`Digits, tstep, t1 `,Digits,tstep,t1);
		fibre:=Continue( CL,curve,x,y,path,t,t_old,t1,fibre);
		R:=R,fibre;
		fibre:=[fibre][-1][2];
		if very_close=1 and tstep<1/(4+CL) then
			tstep:=tstep*2;
			if Digits>dig then
				Digits:=Digits-3
			fi
		elif very_close=-1 then
			tstep:=tstep/2;
			if tstep < Increase_Digits / CL^2 then
				Digits:=Digits + 3
			fi
		fi;
		tstep:=min(tstep,1-t1);
		t_old:=t1
	od;
    [R]
end;

# This procedure continues a given vector along line
#og path pt1*(1-t0)
Continue:=proc(CL,curve,x,y,path,t,t0,t1,fibre,nested)
	local newy,neary,d,der,i,j,x0,x1,m,outy,
	      closeenough,root,mm,newfibre,Fx1;
	global very_close;
option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
if nargs=9 then
	return procname(args,0)
fi;
	#print(curve);
    Digits:=Digits+EXTRA1;
	x0:=evalf(subs(t=t0,path));
	x1:=evalf(subs(t=t1,path));
	Fx1:=subs(x=x1,curve);
    Digits:=Digits-EXTRA1;
if nargs>10 then
	newy:=args[11]
else
#print(Fx1);
	newy:={fsolve(Fx1,y,complex)};
fi;
	newfibre:=newy;
    Digits:=Digits + EXTRA2;
	der:=(x1-x0)*subs(x=x0,Der(curve,x,y,Digits));
	neary:=[seq(i+subs(y=i,der),i=fibre)];
    Digits:=Digits - EXTRA2;
	# neary should be a nearby permutation of newy. It will
	# be used to reorder newy. 
	if nops(neary)<>nops(newy) then 
		error "Cannot determine all roots"
	fi;
	outy:=[];
	closeenough:=true;
	very_close:=1;
	for i in neary do
		# d:=min(seq(evalf(abs(i-j)),j={op(neary)} minus {i}));
		m:=infinity;
		for j in newy do
			mm:=abs(i-j);
			if mm<m then
				d:=m;   # d = second closest
				m:=mm;  # m = closest
				root:=j;
			elif mm<d then
				d:=mm
			fi
		od;
		outy:=[op(outy),root];
		newy:=newy minus {root};
		if m>d/CLOSE2 / CL^2 then
			very_close:=0
		fi;
		if m>d/CLOSE / CL^2 then closeenough:=false;
		fi;
	od;
	newy:=fibre;
	if nested > 0 then
		very_close:=-1;
		Digits:=Digits + 3*nested
	fi;
	if closeenough then
		return [t1,outy]
	fi;
	if nested > 1 then
		newfibre:=NULL;
		newy:=[fsolve(subs(x=x0,curve),y,complex)];
		
		newy:=[seq(newy[i], i=  getPermutation(newy,fibre))]
	else
		newy:=fibre
	fi;
	# The roots are not close enough
	# We'll repeat the continuation, but with
	# smaller stepsize
	userinfo(6,'algcurves',`nested, Digits, x0, x1`,nested,Digits,x0,x1);
	# i = intermediate
	i:=procname(args[1..7],(t0+t1)/2,newy,nested+1);
	i, procname(args[1..6],(t0+t1)/2,t1,[i][-1][2],nested+1,newfibre)
end:

# This procedure computes the first derivative y'(x) on the 
# Riemann surface, using implicit differentiation.

Der:=proc(curve,x,y,DI)
	options remember, 
	    `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	# unapply(-diff(curve,x)/diff(curve,y),x,y);
	evalf(collect(normal(-diff(curve,x)/diff(curve,y)),y),DI)
end:

# This procedure sorts a list of complex numbers 
# according to their distance from a given point x0.
# If two distances are equal, smallest argument comes
# first.

Distsort:=proc(points,x0)
local dec;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	dec:=(s,t)->evalb(
		abs(s-x0)<abs(t-x0) or
		(Im(s-x0)*Re(t-x0)<Im(t-x0)*Re(s-x0) and 
		 abs(s-x0)=abs(t-x0)
		)
		  	 );
	sort(points,dec)
end:

# This procedure gives as output a table of paths, one around each branchpoint, 
# given in the form of a table.

