
# This procedure computes the analytic continuation of 
# a set of function values on the Riemann surface, 
# along a given path. The output is a map from the initial
# points to the end points, given in the form of a table. 
matchUpAC:=proc(perm, edge)  
	return [ edgeData[{op(edge)}][edge]["perm"], map(a-> [a[1], permuteFibs(perm, a[2])],edgeData[{op(edge)}][edge]["AC"])];
end;

	

aContEdge:=proc(edge,x::name,y::name)
	local i,fibre,s1,R,sstep, e,path,sGrid,s_old,s_new,zeroIndex;
	global very_close;	
	# Not including the option remember results in going
	# through these slow calculations more than necessary
	Digits:=UI_LDigits + EXTRA1;
	
	sGrid := [-1,0,1];
	e:=convert~(edgeData[edge]["normDir"],rational);
	path:=edgeData[edge]["path"];

	fibre := getOrderedRoots( x,y,F[1],subs(s=sGrid[1], path));
	 
	R := [ sGrid[1], fibre ]; #seq of lists [t, fiber at t]

	if member('bugAC',_Env_GlobalOptions )  then 
        print("~~~AC:");
		print("digits ", Digits);
		print("normal form of edge", toFloat(e));
		print("order of inital fibre is lex order:"); 
		print(R);               
    fi;
	for i from 2 to nops(sGrid) do
		sstep:=1/2^(BaseLvl);
		s_old := sGrid[ i-1 ];
		s_new := sGrid[ i ];
		if member('bugAC',_Env_GlobalOptions )  then 
			print("continuing from:",s_old);
			print("to ", s_new);        
    	fi;
		while s_old < s_new do
			s1 := s_old + sstep;
			if member('bugAC',_Env_GlobalOptions )  then 
				print("~~~AC");
				print(cat(
					"s_old is updated to ", convert( s_old, string)));
				print("step being used is: ", convert( sstep, string));
				print("digits are: ", Digits);
				print("now continuing");                   
    		fi;	
			fibre :=		
				ContinueDK( x, y,path, s_old, s1, fibre ); 
			R := R, fibre; #add to R	
			fibre := [fibre][-1][2]; #this is the last fiber element [t, firber t]			
			if s1=0 and s_new=0 then
				zeroIndex:=nops([R]);
			fi;
			if member('bugAC',_Env_GlobalOptions )  then 
				print("~~~AC");
				print(cat(
					"continued and now have ",
				convert(nops([R]),string), " data points on edge"));
				print("current step: ", sstep);
				print("current digits: ", Digits);
				print("very close? ", very_close);
				print("minimal digits needed", minDigits ); 
				             
    		fi;		
			if very_close=1 and sstep < 1/2^(max(minDigits,EXTRA1)-EXTRA1+1) then
				sstep:=sstep*2;
				if Digits>minDigits+EXTRA1 then
					Digits:=Digits-EXTRA1
				fi;
			elif very_close=-1 then
				sstep:=sstep/2;
				if sstep < 1/10^(minDigits+1) then
					Digits:=Digits + 3;
				fi
			fi;
			sstep:=min(sstep,s_new-s1);
			if sstep >0 then
				while minStep > sstep do
					minStep:=minStep/2;
					BLvl:=BLvl+1;
				od;
			fi;
			s_old:=s1
		od;	
	od;
	edgeData[edge]["AC"] := [R];
	edgeData[edge]["perm"] := Perm( getPermutation(
			getOrderedRoots(  x,y,F[1],subs(s=sGrid[-1], path)),
			[R][-1][2]));
		
	edgeData[edge]["zeroIndex"]:= zeroIndex;
	
end;
#Continue:=proc(CL,curve,x,y,path,t,t0,t1,fibre,nested)
ContinueDK:=proc(x,y::name,path,s0::rational,s1::rational,fibre,nested:=0)
	local j,newy,neary,dOld,der,derder,i,x0,x1,mm,newfibre;
	global very_close;
	#print(curve);	
	x0:=normal(subs(s=s0,path)); 
	x1:=normal(subs(s=s1,path));
	Fx1:=convert(subs(x=x1,F[1]),horner,y);
	if member('bugAC2',_Env_GlobalOptions )  then 
		print("~~~AC2");
		print(cat(
			"nested ", convert(nested, string)));
		print("x0 ", x0,s0);
		print("x1 ", x1,s1);

		print(normal(x1-x0));
		print("deltax", toFloat(normal(x1-x0)));
		e1:=subs(s=-1,path);
		e2:=subs(s=1,path);
		print(plots[display]({plot(complexLToPtL([e1,e2]),color="black"), 
		plots[pointplot](complexToPt(x1), color="green"),
		plots[pointplot](complexToPt(x0), color="blue")}));
		print("function at x1", Fx1);                   
	fi;	
	if nargs>7 then
		newy:=args[8]
	else
		newy:={};
		DI:=Digits;
		Digits:=UI_HDigits;
			
		newy:=toFloat~(convert(polySolve(toFloat~(Fx1),y),set));
		
		if member('bugAC2',_Env_GlobalOptions )  then 
			print("~~~AC2");
			print("new fibre computed at digits=",Digits);
			print(newy);                 
		fi;
		Digits:=DI;
	fi;
	newfibre:=newy;
	Digits:=Digits + EXTRA1;
	der:=(x1-x0)*subs(x=x0, F["Yx"]); #there is cancelation here
	
	 #there is cancelation her
	der:=convert((y+der),horner);
	neary:=evalf([seq(subs(y=i,der),i=fibre)]); #this is the order
	if member('bugACErr',_Env_GlobalOptions )  then 
		print("~~~AC2");
		print("approx in rationals");
		print(der);
		print("in floats");       
		print(neary);
		print("should match up with");
		print(newy);
		print("fibre used");
		print(fibre);
		der:=evalf((x1-x0)*subs(x=x0, F["Yx"]));
		print(toFloat~([seq(i+subs(y=i,der),i=fibre)]));
	fi;
	Digits:=Digits - EXTRA1;
	# neary should be a nearby permutation of newy. It will
	# be used to reorder newy. 
	if nops(neary)<>nops(newy) then 
		error "Cannot determine all roots"
	fi;
	mm := DK(Fx1,y,neary,1/(2*numSheets),newfibre);
	very_close:=1;
	if mm <> FAIL then #it must be close enough to determine it is the right one?
		if member('bugAC',_Env_GlobalOptions )  then 
			print("~~~AC2");
			print("x1",evalf(x1));
			print("edge,",evalf~([subs(s=-1,path),subs(s=1,path)]));
			print("values from dk ");
			print(mm);
			print(newfibre); 			              
		fi;	
		if mm[2]<1/5 then
			very_close:=-1;
		elif mm[2]<1 then
			very_close:=0;
		fi;
		#if min digits of dis between fiber is on the order of min digits or less
		return [ s1,mm[1] ]
	fi;
	if member('bugACF',_Env_GlobalOptions )  then 
		print("~~~AC2FAIL");
		print("digits",Digits);
	fi;
	if nested > 0 then
		very_close:=-1;
		Digits:=Digits + 3*nested
	fi;
	if nested > 1 then
		newfibre:=NULL;
		DI:=Digits;
		Digits:=max(UI_HDigits,DI);
		newy:=evalf(convert(polySolve(subs(x=x0,F[1]),y),list));
		newy:=[seq(newy[i], i=  getPermutation(newy,fibre))]; #list 1 to list 2
		Digits:=DI;
	else
		newy:=fibre
	fi;
	
	i:=procname(args[1..4],(s0+s1)/2,fibre,nested+1); #cut interval into two and apply again
	i, procname(args[1..3],(s0+s1)/2,s1,[i][-1][2],nested+1,newfibre)
end:







aContEdgeInt:=proc(e)
	local i,fibre,s1,dig,t_old,R,CL, edge, lvl,N,completeSGrid,k,s0;
	global very_close;
	    `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	# Not including the option remember results in going
	# through these slow calculations more than necessary.
	lvl := intPrimitives["Glvl"]-intPrimitives["Blvl"];
	N := intPrimitives["N"];
	completeSGrid :=[ 
		intPrimitives["sVals"][1],
		seq(intPrimitives["sVals"][i*2^lvl+1],i=1..floor(N/2^lvl))];
	
	edge:=edgeData[e]["normDir"];
	for k from 1 to numSheets do
		edgeData[e]["int"][k] := table();
		edgeData[e]["int"][k]["term"] := -1;
		edgeData[e]["int"][k]["vec"] := Vector(genusCurve);
	od;
	Digits := UI_HDigits;
	fibre := toFloat~( subs( x= edge[1], F[3] ) ); #this may be incorrect if 1 is not in
	R := [ -1, fibre ]; 
	for i from 1 to 2*N do
		if i < N then		
			s0 := evalf(-completeSGrid[-i]);	
			s1 := evalf(-completeSGrid[-i-1]);
		elif i=N  then
			s0 := evalf(-completeSGrid[-i]);	
			s1 := evalf(completeSGrid[-i-1]);
		elif i=N+1 then
			s0 := evalf(completeSGrid[i-N]);	
			s1 := evalf(completeSGrid[i-N+1]);
		fi;
		fibre :=	ContinueDK(x,y, path, s0,s1,fibre); 
		R:=R,fibre; #add to R		
		fibre:=[fibre][-1][2]; 
		if i<> 2*N then 	
			updateVal1(edge[1],edge[2],s1,-N+1-i,fibre);
		fi;
	od;	
	edgeData[e]["AC"]:=[R];
	edgeData[e]["perm"] := Perm( getPermutation(
			toFloat~(subs(x=edge[2],F[3])),
			[R][-1][2]));
	for k from 1 to numSheets do
		edgeData[e]["int"][k]["term"]:=edgeData[e]["perm"][k];
		edgeData[e]["int"][k]["vec"]:=evalf( edgeData[e]["gam"]*2^(-intPrimitives["Blvl"])*edgeData[e]["int"][k]["vec"]);
	od;
end;
reverseEdgeAC:= proc(e)
	local invPerm;
	invPerm:=pInv(edgeData[{op(e)}]["perm"]);

	return [invPerm, reverseList(
			map(a-> [-a[1], permuteFibs(invPerm,a[2])],
    		edgeData[{op(e)}]["AC"])
		)];
end proc;



	#if evalf(1-s0) <= -7/8  then 
	#	x0:= evalf(evalf(e1*s0/2+e2)- evalf(e2*s0/2));
	#	x1:= evalf(evalf(e1*s1/2+e2)- evalf(e2*s1/2));
	#else
	#	x0:= evalf(evalf(e1*s0/2)+ evalf(e2*(2-s0)/2));
	#	x1:= evalf(evalf(e1*s1/2)+ evalf(e2*(2-s0)/2));
	#fi;
	


Continue:=proc(x,y,ePath,s0,s1,fibre)
	local newy,neary,d,der,i,j,x0,x1,m,outy,nested,
	      closeenough,root,mm,newfibre,Fx1,derder;
	global very_close;

	if nargs=6 then #if there are no extra args then it is not nested so return it as not nested
		return procname(args,0) #bc
	fi;
	nested:=args[7];
     #increase digits as we integrate along this path by 3 just to be accurate in terms of x0,x1
	x0:=evalf(subs(s=s0,ePath)); 
	x1:=evalf(subs(s=s1,ePath));
	if nargs>7 then #if args >10 then this is nested and there is added to a known fiber
		newy:=args[8]
	else #otherwise we need to compute the roots of x1
		newy:=scrub~( { op(evalf(subs(x=x1,F[3])))});	
	fi;
	newfibre:=newy;
	Digits:= Digits+EXTRA1; 
	der:=(x1-x0)*subs(x=x0, F["Yx"]); #there is cancelation here
	derder:=(x1-x0)^2*subs(x=x0, F["Yxx"]/2); 
	neary:=evalf~([seq(i+subs(y=i,der+derder),i=fibre)]); #this is the order
	# take elements of the fiber which may be a seq of lists
	# and compute y+der this is the being near portion
    Digits:= Digits-EXTRA1; 
	# neary should be a nearby permutation of newy. It will
	# be used to reorder newy. 
	if nops(neary)<>nops(newy) then #the points collapsed and so too close to branch point
		error "Cannot determine all roots"
	fi;
	outy:=[];
	closeenough:=true; #close enough means we have n distinct fibers
	very_close:=1; 
	for i in neary do #each approx 
		# d:=min(seq(evalf(abs(i-j)),j={op(neary)} minus {i}));
		m:=infinity;
		for j in newy do #get closest y vector 
			mm:=abs(i-j);
			if mm<m then #mm< closet
				d:=m;   # d = second closest
				m:=mm;  # m = closest
				root:=j;
			elif mm<d then #mm< second closest
				d:=mm
			fi
		od;
		outy:=[op(outy),root]; #store it in order
		newy:=newy minus {root}; #remove from consideration
		#idea is to measure how close the two closest possible choices are:
		# as we have better accuracy the correct on will converge to zero and the otherone will converge to a number
		# provided we are close enough so that the second order behavior becomes negliable
		if 2*m>d then # the closest dist is larger than the second/20 / (dig-8), then it is not very close m 20(Digits)> d #not a magnitude smaller
			very_close:=0 #measures if we are close enough to say that this is the right choice
		fi;
		if (3/2)*m > d then closeenough:=false; # even a fifth of the above quantity is to large 
		#definitely implies #the distance between the approx and solution is not zero and
		#m is not even a magnituge smaller than d
		fi;
	od;
	newy:=fibre;
	if nested > 0 then # if nested then it is not very close becuase we needed to cut the interval
		very_close:=-1;
		#print("not close enough");
		#print(nested);
		Digits:=Digits + 3*nested #increase digits according to nesting
	fi;
	if closeenough then #it must be close enough to determine it is the right one?
		return [s1,outy]
	fi;
	if nested > 1 then  #first try spliting then increase digits too 
		newfibre:=NULL;
		newy:=
			[ op( evalf( polySolve(subs(x=x0,F[1]),y)) )]; #need more accuracy here now
		newy:=[seq(newy[i], i=  getPermutation(newy,fibre))]
	else
		newy:=fibre
	fi;
	# The roots are not close enough
	# We'll repeat the continuation, but with
	# smaller stepsize
	#userinfo(6,'algcurves',`nested, Digits, x0, x1`,nested,Digits,x0,x1);
	# i = intermediate
	i:=procname(args[1..4],(s0+s1)/2,newy,nested+1); #cut interval into two and apply again
	i, procname(args[1..3],(s0+s1)/2,s1,[i][-1][2],nested+1,newfibre)
end:

# This procedure computes the first derivative y'(x) on the 
# Riemann surface, using implicit differentiation.

# This procedure sorts a list of complex numbers 
# according to their distance from a given point x0.
# If two distances are equal, smallest argument comes
# first.

#get rep in yVec of neary
extractFib:=proc(neary,yVec)
	local i, mm,II,m:=infinity;	
	for i from 1 to nops(yVec) do 
		
		mm:=evalf(abs(yVec[i]-neary));

		if mm<m then 
			m:=mm;  # m = closest
			II:=i;			
		fi;
	od;
	return yVec[II],II;
end:
# This procedure gives as output a table of paths, one around each branchpoint, 
# given in the form of a table.


