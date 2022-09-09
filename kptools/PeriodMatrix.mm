# PeriodMatrix
#-make it work with unprotect
# DESCRIPTION
# working idea: compute mono to get hurwitz, get linear combos from homology, get the associated hurwitz eles and parametrize them. Then use intgrate comand against the differentials and Y as is.

PeriodMatrix := proc(f, x::name, y::name)
local pm,alpha,beta,g,F,A,i,j,m,DataMono, DataHomo,dx ;option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved. Authors: B. Deconinck and M. van Hoeij.`;
    #$include "C:/Users/msben/Monodromy-Voronoi/ParametricProduct.mm"	
    #step1: get the needed info
    #print("get monodromy, homology, and differentials");
    F:=collect(primpart(f,{x,y}),{x,y}, 'distributed', Normalizer);
    g:=algcurves[genus](F,x,y);
	if g<=0 then error "Input must have positive genus" fi;
    if indets(f,'name') <> {x,y} or not type(f,polynom(anything, [x, y])) then
		error "Curve should be a polynomial in both variables %1 and no other variables", {x,y}
	elif indets(f,float)<>{} then
		error "No floating point coefficients allowed"
	fi;
	_Env_GlobalOptions:=[args[4..-1]];
	pm:=Pmatrix(Digits,F,x,y);    
	if member('Riemann', _Env_GlobalOptions) or member('normalized',_Env_GlobalOptions) then
		alpha:=delcols(pm,g+1..2*g,outputoptions=[]);
		beta:=delcols(pm,1..g,outputoptions=[]);
		pm:= linsolve(alpha,beta,method=none,free=_t0,
			conjugate=true,
			inplace=false,
			outputoptions=[]); # works better than evalm(alpha^(-1))
		if g>1 then
			m:=max(seq(seq(abs(pm[i,j]-pm[j,i]),i=j+1..g),j=1..g-1));
			i:=cat("About: ",convert(evalf(-log10(m),3),string),
				" accurate Digits");
			if m>10^(3-Digits) then
				WARNING( i ) # Shouldn't happen.
			else
				userinfo(2,'algcurves',i)
			fi
		fi;
		# Test if Imaginary part is positive definite.
		m:=Im(pm);
        #print(m);
		if not LinearAlgebra[IsDefinite](m,query='positive_definite') then
			error "Imaginary part not positive definite" 
		fi;
	fi;
		pm

	
end:
ErrorEstimate:=proc(e)
gammae:=e[1]*(1-t)/2 + e[2]*(1+t)/2;
gamma_0:=diff(gammae,t);
#given NN, h, abiscuss  <-edge data
#1) get your AC data, y=gamma analytic continuation 
#2) get your abiscussis=gamma(G(t))
#2) construct the error functions
#for each sheet compute the sum
#if all of them satisfy criteria good and pass error estimate
#otherwise add more values
tvals:=computedTValues;
h:=tvals[1];
l:=2;
while abs(subs(t=l,G))>= 10^(-Digits) do

for j from N tp h

nN:=NN;
tvals:=computedTValues;
h:=tvals[1];

yabiscus:=1/(e^(Pi/2*sinh(t))*cosh(Pi/2*sinh(t)));


end proc;
Pmatrix:=proc(DI,curve,x,y)
	local h,g,pm,kappa,i,j,k,contours,diffs,M, Diffs:=[];
	options remember;        
	h:=Homology(curve,x,y,`give paths`);
    
    kappa:=h[1]['linearcombination'];
    contours:=h[1]['cycles'];
	g:=nops(diffs);
	pm:=Matrix(g,2*g);
    #print(M);
    #print(contours);
    #print(h[2]);
    for i to 2*g do		   
		:=+kappa[j,k]*Contourintegrate(
			        curve,x,y,contours[k],
			        h[1]['basepoint'],h[1]['sheets'],h[2],h[3]);
		        fi;
            od;
	    od;
    od;
	forget(Contourintegrate);
	forget(Pathintegral);
	forget(Partintegral);
	forget(integrand2);
	pm, Diffs, [h[1]['basepoint'],h[1]['sheets']]
end:

Contourintegrate:=proc(curve,x,y,contour,x0,sheets,sites,paths)
	local n,m,i,integ,bpoint,firstsheet,lastsheet,perm,cycle;
	options remember, 
	    `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	integ:=0;
	n:=nops(contour)/2;
	cycle:=[op(contour),contour[1]];
    for i from 1 to n do
	    bpoint:=cycle[2*i][1];
	    perm:=cycle[2*i][2];
	    firstsheet:=cycle[2*i-1];
	    lastsheet:=cycle[2*i+1];
	    perm:=Putfirst(perm,firstsheet);
	for m from 1 to Whereis(perm,lastsheet)-1 do
		# Go around the branchpoint, and add the integral
		# over that path to integ:
	    integ:=integ+
		Pathintegral(curve,x,y,x0
		  ,sheets[perm[m]],sheets[perm[m+1]],bpoint,sites,paths);		
	od;
od;
	integ
end:

Pathintegral:=proc(curve,x,y,x0,y0,y1,bpoint,sites,paths)
	local i,n,integ,yi,yj,parint,rpoints,ipoints,
	      rmax,imin,imax,v,path;
	options remember, 
	    `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;

	integ, yi, yj := 0, y0, y1;
	# go to branchpoint
    #print(sites);
    #print(bpoint);
    #print(searchList(bpoint,sites));
    path:=paths[searchList(bpoint,sites)];
    n:=nops(path);
	for i to n do
	    parint:=Partintegral(curve,x,y,path[i][1]*(1-t)+path[i][2]*t ,t,yi);
	    yi:=parint[2];
	    integ:=integ+parint[1];
	od;
	if abs(yi-yj)>(1+abs(yi))*10^(-iquo(Digits,2)) then 
		error "Path was not followed correctly during the integration"
	fi;
	integ
end:

Partintegral:=proc(curve,x,y,path,t,y0,acc)
	local k,integrand,N;
	#Add the edges on the path together using the one starting at y0.
	#return the vector plus what you should get for y1
end:



    
Whereis:=proc(list,element)
	local i;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	for i from 1 to nops(list)+1 do if list[i]=element then
		break
	fi od;
	i
end:
Putfirst:=proc(L::list,element)
	local pos,i;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	pos:=Whereis(L,`if`(nargs=1,Smallest(L),element));
	[seq(L[i],i=pos..nops(L)),seq(L[i],i=1..pos-1)]
end: