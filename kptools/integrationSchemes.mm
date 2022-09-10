








lineInt:=proc(curve,differential,x,y,contour,x0,sheets)
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
	    integ:=integ+
			lineInt1(curve,differential,x,y,x0
		  		,sheets[perm[m]],sheets[perm[m+1]],bpoint);		
		od;
	od;
	integ
end:

lineInt1:=proc(curve,differential,x,y,x0,y0,y1,bpoint)
	local i,n,integ,yi,yj,parint,rpoints,ipoints,
	      rmax,imin,imax,v,path;
	options remember;
	
	integ, yi, yj := 0, y0, y1;
	
    ppath:=pathData[bpoint]["path"];
	
	#print(yj);
	
	for e in ppath do 
		
	    parint:=lineInt2(curve,differential,op(e),yi,bpoint);
	    yi:=parint[2];#the y value we end at
	
		#print("yi update");
		#print(yi);
	    integ:=integ+parint[1];		
	od;
	
	if abs(yi-yj)> 0 then 
		error "Path was not followed correctly during the integration"
	fi;
	integ
end:

lineInt2:=proc(curve,differential,pt1,pt2,y0,bp)
	local k,integrand,n;
	
	# userinfo(4,'algcurves',`Integrating over path`,path);
	
	Digits:=2*max(10,Digits);
	A1:=pathData[bp][[pt1,pt2]];
	
	sheetnum := op(getPermutation(A1[1][2],[y0])); # Number of current sheet.
	
	sumVal:=0;
	for i from 1 to nops(tInterval) do
		xi:= evalf( pt1*(1- tInterval[i])/2+pt2*(1+ tInterval[i])/2);

		for ac in A1 do	
			#print(ac);		
			if ac[1]= tInterval[i] then
				yval:= ac[2][sheetnum];
				break;
			fi;
		od;
		if i=1 or i=nops(tInterval) then
			f:= ( pt2-pt1)/2*subs(x=xi,y=yval, differential);
		elif i mod 3 = 0 then
			f:= 2* ( pt2-pt1)/2*subs(x=xi,y=yval, differential);
		else 
			f:= 3* ( pt2-pt1)/2*subs(x=xi,y=yval, differential);
		fi;
	sumVal:= sumVal + f;
	od;
	sumVal := evalf(3/8*2/(nops(tInterval)-1)*sumVal);
	sumVal, A1[-1][2][sheetnum];
end: