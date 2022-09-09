#fix center problem
#clean up code
#find better way to do intersections
# Do very quality test on Hurwitz

Monodromy := proc(curve, x, y) 
    options remember;
    local F,sigma,C,i,b,mono;	
	
	if indets( evalf(curve),'name') minus {x,y}<>{} then	#if other symbols are used we got a problem	
		error "Only 2 variables allowed"
	elif indets(curve,float)<>{} then #if its a float we got a problem	
		error "No floating point coefficients allowed"
	elif member(group,[args[4..-1]]) then  #if group is mentioned, return the permgroup
		C:=procname(args[1..3],op(subs(group=NULL,[args[4..-1]])));
		'permgroup'(nops(C[2]),{seq(i[2],i=C[3])})
	elif args[nargs] <> `give paths` then #if not give paths is mentioned get the paths we pass it back in and ommit the paths produced
		C:=procname(args,`give paths`);
		[C[1],C[2],C[3]];
    else #otherw
	    F := collect(primpart(curve,{x,y}),{x,y},`distributed`,Normalizer);
	    Digits := max( 10,Digits);
	    _EnvA1 := [args[ 4..nargs ], _Env_algcurves_opt ];		
        Monodromy_Processed( Digits,F,x,y);   
    fi
end:

Monodromy_Processed:=proc(DI,F,x,y)
	if degree(F,y)=0 then
		error "curve is independent of the second variable"
	elif degree(F,x)=0 then
		error "curve is independent of the first variable"
	elif genus(F,y)=1 then
		Linearmonodromy(F,x,y)
	else 
        if degree(F,x)<3 then
		    userinfo(1,'algcurves',
		    `Consider switching the roles of the variables. The curve has degree <=2 in the first variable.`);
	    fi;
		Generalmonodromy(args[2..nargs])
	fi
end:
# This procedure gives the monodromy information of the Riemann surface
# of a single-valued function.

Linearmonodromy:=proc(curve,x,y)
	local xvalue,preimages;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	xvalue := 0;
	while subs(x = xvalue, lcoeff( curve, y))=0 do
		xvalue := xvalue - 1
	od;
	preimages := [  simplify(fnormal(evalf( polySolve( subs( x=xvalue, curve ), y ) )),zero) ];
	[xvalue,preimages,[]]
end:


solveMono:= proc(expression,var)
	return op(simplify(fnormal(evalf( SolveTools[Polynomial]( expression,var))),zero));
end:
Generalmonodromy:=proc(curve,x,y)
local digRequired,vToEData,sToCData, eToSData,preimages,bps,theMono, ACDATA,F, NumOfProblemPts, E:={}, P:=[], sites_ordered, listOfSites, canonicalVertex, problemPoints, Y_v, distRequired,nxtPt,sIndex,sVector,s,sitesOfEdge,
		directionVector,neighbor,possibleNxtPts,crtNeighbors, reachedEnd,pstPt,crtPt,si,p1,sigma,C,i,b, mono,infinityIsBp;
	uses GroupTheory; #fix this  
    Digits := max( 10, Digits);
    Digits := Digits+10;
    #get discrim pts;
    
	F := ProcessPoly( curve, x, y );  
    problemPoints :=   [ op( { solveMono(  F[2] , x ) ,  solveMono(  F[3] , x )} ) ];
    NumOfProblemPts := nops(problemPoints);
	
	listOfSites := [ op(sortComplex( problemPoints,'normal') ), ComputeBoundary( problemPoints, NumOfProblemPts) ]; #maybe check if embedded real
    if member('bug',{args}) then
            print("List of sites ");
            print( listOfSites ); 
    fi;
    Digits:= Digits - 10;    
    
    digRequired := round(-log10(mindist(map( a-> [Re(a), Im(a)],listOfSites))[1]))+1; 
    if digRequired > 0 then
        Digits := Digits + digRequired;
    fi;
    if member('bug',{args}) then
        print("digits are resolved with");
	    print( 10^(-digRequired) );
    fi;

    Dual := voronoiDual( map( a-> [ Re(a), Im(a)], listOfSites ) );
	if member('bug',{args}) then
        print("Voronoi dual computed");
    fi;
	E, vToEData, sToCData, eToSData := ProduceVoronoi( Dual, listOfSites, NumOfProblemPts, digRequired );
    
    canonicalVertex := sortComplex( op~(E),'normal')[1];
    preimages := [ solveMono( subs(x=canonicalVertex,curve),y) ];

	sites_ordered := sortComplex( listOfSites[1..NumOfProblemPts], `Arg`, canonicalVertex );

	for s in sites_ordered do
        sToCData[ s ] := [ sToCData[ s ], ProduceCycle( sToCData[ s ], s, vToEData) ];  
        #print(sToCData[ s ]);
    od;
	
	if member('bug', {args}) then
		print( plots[display]( plot( map(a->[[Re(a[1]),Im(a[1])],[Re(a[2]),Im(a[2])]],E ), color ="red" ), 
		plots[pointplot]( map(a->[Re(a),Im(a)],[canonicalVertex, op(listOfSites[1..NumOfProblemPts])])),axes=NORMAL,scaling=CONSTRAINED, 
		title=typeset("Paths chosen for the analytic continuation of y(x)"), 	
        labels=[Re, Im]) );
	fi;
	P := ProduceHurwitz( canonicalVertex, sites_ordered, NumOfProblemPts, vToEData, sToCData, eToSData, E ); 
	
	#gc();
	if member('showpaths', _EnvA1) and not member('bug', {args}) then
		print( plots[display]( plot( map(a->[[Re(a[1]),Im(a[1])],[Re(a[2]),Im(a[2])]],E ), color ="red" ), 
		plots[pointplot]( map(a->[Re(a),Im(a)],[canonicalVertex, op(listOfSites[1..NumOfProblemPts])])),axes=NORMAL,scaling=CONSTRAINED, 
		title=typeset("Paths chosen for the analytic continuation of y(x)"), 	
        labels=[Re, Im]) );
	fi;
	#can delete a lot here
	#print("generating mondromy");
	infinityIsBp, bps, theMono, ACDATA := GenerateMonodromy(F[1],E, sites_ordered, canonicalVertex, preimages, NumOfProblemPts, P);
	
	if infinityIsBp then 
		 	sites_ordered:=[op(sites_ordered),infinity];
			#Sprint("inf");   
            p1 := [];
                     
            si := bps[-1];
            crtPt := canonicalVertex;
            pstPt := canonicalVertex;
			reachedEnd:=false;
            while not reachedEnd do    
				#print(vToEData[crtPt]);           
            	crtNeighbors :=  op~( vToEData[crtPt] ) minus { crtPt, pstPt };
            	#print(crtNeighbors);
            	possibleNxtPts := crtNeighbors;               
            	for neighbor in crtNeighbors do
                	#print("considering neighbor");   
            		directionVector := < Re(neighbor - crtPt), Im(neighbor - crtPt), 0 >;
                	#print("site associate");
                	sitesOfEdge := eToSData[ { crtPt, neighbor } ];
                  	#print(sitesOfEdge);
                 	#print(possibleNxtPts);
                	if nops( sitesOfEdge ) = 2 then
                    	possibleNxtPts := possibleNxtPts minus {neighbor};
                  	else #s1,s2 are new so not neighbor
                    	s := sitesOfEdge[1];	
                     	sVector := < Re(s - crtPt), Im(s - crtPt), 0 >;	
                     	sIndex := searchList( s, sites_ordered );
                    	#print(sIndex);
                    	if  cross(directionVector, sVector)[3] > 0 then #site is to the right                         
                        	possibleNxtPts := possibleNxtPts minus {neighbor};
                    	fi; #sIndex< i so this site is old. This means this is the neighbor
                	fi;                  
            	od;
            	nxtPt :=possibleNxtPts[1]; 
				if nxtPt = canonicalVertex then
				reachedEnd:=true;
			   fi;			       
               p1 := [ op(p1) , [ crtPt , nxtPt]];
               pstPt := crtPt;
               crtPt := nxtPt;              
               #print( plots[display](plot(p1, color="black"),plot(map(a->[op(a)],E), color="yellow")));               
               #print(display([ plot(G minus {op(p1)}, color="red"),pointplot(crtNeighbors, color="blue"), plot(p1, color="black"), pointplot(canonicalVertex, color="red", symbol=solidbox),pointplot(nxtPt, color="green"), pointplot({si},color="red", symbol=soliddiamond)]));

            od;
			
            
            P := [ op(P), [ op(p1) ] ];  
			fi;
	
	return [canonicalVertex, preimages, convert(theMono, list), sites_ordered,P,ACDATA];
end proc;
#IS this the best way to do this for the polynomial
ProcessPoly := proc(curve, x, y)
   local F1 := curve, F2 := evala(normal( resultant(F1, diff(F1,y), y) ) ), F3 := evala( normal( lcoeff(F1,y)));    	  
   if degree(F2,x) > 0 then     
      F2 := evala(quo(F2 , gcd( F2, diff(F2, x) ),x));
    fi;              
    if degree(F3,x) > 0 then
       F3 := evala( quo( F3 , gcd(F3,diff(F3,x) ),x));
    fi;
    F2 := evala(quo(F2, gcd(F2,F3), x));  	
   return [F1,F2,F3];
end proc:


ComputeBoundary := proc(ll::list,size::integer)
   local center,k,r;
   
	
	
	center:=  simplify(fnormal(evalf( add( ll[k], k=1..size )/size)),zero);
	r := simplify(fnormal(evalf( 2 * max( {  op( map( a -> abs(  a - center ), ll) ), 1 } ) )),zero); 
	
    if member('bug',{args}) then
        print("Boundary computed, with center");
        print(center);
        print("with radius");
        print(r);
    fi;
	#print(seq(simplify(fnormal( center + r * exp( (1/3*Pi*k+Pi/6)*I )) ,zero),k=1..6 ));
	return seq(simplify(fnormal(evalf(
        center + r * exp( (1/3*Pi*k+Pi/6)*I)  
    )),zero),k=1..6 );	
end proc:

#Do a better job getting verts here
#Input: dualnay triangles: list of lists comprised of 3 integers and the integers rep the site indicies
#sites
#output: Edges
# sites to pre cell, verts to edges, edges and sites associated to them --some of this is beyond ness but since we can expect num of sites to not get to large there is
# more of a constraint on speed then on space and so we should worry less about space (i.e. extra objects).
ProduceVoronoi := proc(triangleSet,sitelist, size,req)
    local triangleToVert:=[], siteToTriangles := table(),edgeSet, siteToCellData := table(), vertexToEdgeData := table(),
	edgeToSiteData := table(),i, j, l,s, Ti,sj, si, ci, indexSet,sameCenters,cj,c,E,triangles,Tk,neighborsInTk,pt1, possibleTj, Tl,pt2;
    ########GET THE VERTICIES
    #print(triangleSet);
    for i from 1 to nops(triangleSet) do	
		Ti := triangleSet[i]; #get a triangle
		for j in Ti do #for each site in the triangle see if we allredy have it
			sj := sitelist[j];
			if assigned( siteToTriangles[ sj ] ) then #if we do, make sure to keep track of the triangle it was apart of
				siteToTriangles[ sj ] := siteToTriangles[ sj ] union { i }; 
	 		else #otherwise keep track of the triangle after assignment
				siteToTriangles[ sj ] := { i };
	 		fi;
	 	od;
        #for Ti, find the appropriate lines       
        l1,l2:=getLine(Ti, sitelist);        
		XInt:=fsolve(l1=l2);
        YInt:=subs(x=XInt,l2);	
		ci := [XInt,YInt];
		ci := simplify(fnormal(Complex(op(ci))),zero);
		triangleToVert := [ op( triangleToVert ), ci  ];	 					
	od;    
    ########process centers of triangles
	triangleToVert := Array(1..nops(triangleToVert),triangleToVert); #may be so many triangles that an array must be used
    indexSet := { seq(i, i=1..nops(triangleSet)) };
	while indexSet <> {} do #while we still have triangles
		i := min(indexSet); #chose the smallest
		ci := triangleToVert[i]; # get the center of the triangle: the vert for the edge
		sameCenters := []; 
		for j from 1 to nops(triangleSet) do #go through the remaining edges and log which ones should be the same
			cj := triangleToVert[j];
            
			if fnormal( evalf( abs( ci-cj ) ) ) < 10^(-req-1) then
				sameCenters:=[ op(sameCenters), j];
			fi;
		od;		
        #Once this is done, average to get a better approximation to what the true "center" of these triangles are.
		c := simplify(fnormal(evalf(
            add( 
                triangleToVert[ sameCenters[j] ], j=1..nops(sameCenters)
            )/nops(sameCenters)
        )),zero);		
		for j in sameCenters do #replace the corresponding centers with this updated one and remover them from consideration
			triangleToVert[ j ]:=c;
			indexSet := indexSet minus{ j };
		od;		
	od;
    ########construct edges
	E:={};	
    #This could be more efficient..
    
    for i from 1 to size do #for each site not a boundary!		
        si := sitelist[i];
		triangles := siteToTriangles[ si ]; #get all triangles associated to it
        for k in triangles do   #for each triangle associated to si          
			Tk := triangleSet[ k ]; 
            pt1 := triangleToVert[ k ];	  #get center 
            
	 		neighborsInTk := { op( Tk ) } minus { i };  #now get other sites of triangle   
           	for j  in neighborsInTk do #for each of the nieghbors                
                sj := sitelist[j];                
                l := op( ( siteToTriangles[ sj ] intersect triangles) minus {k}); #get the triangle that shares the edge si,sj
                pt2 := triangleToVert[ l ];                 
				if nops({ pt1, pt2 })>1 then #if the edge is not a point		
				    E := E union { { pt1, pt2} }; #then add it to the edges
				    edgeToSiteData[  { pt1, pt2 }  ] := { si, sj } intersect { op( sitelist[ 1..size ]) }; #and record the sites (if nonboundary) asociated to it	
       			    for s in edgeToSiteData[  { pt1, pt2 }  ] do	#and for each of the sites record the edge	 
					    if assigned( siteToCellData[ s ] ) then 
					    	siteToCellData[ s ] := siteToCellData[ s ] union {  {pt1, pt2}  };
				   	    else
       			            siteToCellData[ s ] := {  {pt1, pt2}  } ;
				   	    fi;          
         		    od; 
                    #record edges associated to the vert
				    if assigned( vertexToEdgeData[ pt1 ] ) then 
         		    	vertexToEdgeData[ pt1 ] :=  vertexToEdgeData[ pt1 ] union { { pt1, pt2 } };
				    else
       			        vertexToEdgeData[ pt1 ] :=  { { pt1, pt2 } };
    			    fi;
    			    if assigned( vertexToEdgeData[ pt2 ] ) then 
      			        vertexToEdgeData[ pt2 ] :=  vertexToEdgeData[ pt2 ] union { { pt1, pt2 } };
    		        else
    			        vertexToEdgeData[ pt2 ] :=  { { pt1, pt2 } };
    		        fi;
                fi;	 
			od;
		od;             
    od;	
	E, vertexToEdgeData,siteToCellData,edgeToSiteData;
end proc:

ProduceHurwitz := proc(canonicalVertex, sites,size,vertexToEdgeData,siteToCellData,edgeToSiteData,E) ; 
    local i,j, si, siData, crtPt, pstPt, nxtPt, p1, p2, p3, pt, neighbor, crtNeighbors, possibleNxtPts, directionVector,
    sitesOfEdge, s1, s2, s, s1Index, s2Index, sIndex, s1Vector, sVector,mm,  P:=[], 
    isToRight := proc( v1, v2 )
        if member('bug',{args}) then
            print(evalf( v1[1]*v2[2]- v1[2]*v2[1]));
            print(evalf[Digits]( v1[1]*v2[2]- v1[2]*v2[1]));
        fi;        
        return evalb(evalf[Digits]( v1[1]*v2[2]- v1[2]*v2[1])< 0);
    end proc;

    for i from 1 to  size  do #for each site in order
        si := sites[i];
        siData :=  siteToCellData[ si ];  
        crtPt := canonicalVertex;   #start from canonnical vertex
        pstPt := canonicalVertex;       
            
           
        p1 := []; # we need to construct path p1p2p1^{-1} to cell, but we already have path on cell
        p2 := Array(1..nops(siData[2]),siData[2]); 
        
        while not ( crtPt in op~(siData[1]) ) do  #while we travel and are not on the cell
            crtNeighbors :=  op~( vertexToEdgeData[crtPt] ) minus { crtPt, pstPt }; #at crt location, get neighbors
            if member('bug',{args}) then
                print(crtNeighbors);
            fi;
            
            possibleNxtPts := crtNeighbors;               
            for neighbor in crtNeighbors do #choose a neighbor
               
                directionVector := < simplify(fnormal(Re(neighbor - crtPt)),zero), simplify(fnormal(Im(neighbor - crtPt)),zero), 0 >; #get direction vector
                if member('bug',{args}) then
                    print("considering neighbor");
                    print(directionVector);
                fi;
               
                sitesOfEdge := edgeToSiteData[ { crtPt, neighbor } ]; # get the sites associated to the edge in question
                #use these sites to get an understanding of path
                if nops( sitesOfEdge ) = 2 then  
                    s1 := sitesOfEdge[1];
                    s2 := sitesOfEdge[2];	
                    s1Vector := < simplify(fnormal(Re(s1 - crtPt)),zero),  simplify(fnormal(Im(s1 - crtPt)),zero), 0 >;
                    s1Index := searchList( s1, sites );
                    s2Index := searchList( s2, sites );
                     if member('bug',{args}) then
                    print(s1Index);
                    print(s2Index);
                    fi;
                    
                    bool1 := isToRight(directionVector,s1Vector);
                    if s1Index < i and s2Index < i then #both sites old
                        possibleNxtPts := possibleNxtPts minus {neighbor}; #old edge so we do not travel on it
                    elif s1Index > i and s2Index < i then #s1 is new, s2 old:s1 should be on the left
                        if  bool1 then #s1 new and to the right so cannot be the neighbor and actually we have a problem
							error "failed to resolve points: found a new site on the wrong side of the edge";
                        fi;
                    elif s1Index < i and s2Index > i then #s1 is old, s1 is new
                        if not bool1  then # s1 should be on the right
                            error "failed to resolve points: found a new site on the wrong side of the edge";
                        fi;
                    else #s1,s2 are new so not neighbor            
                        possibleNxtPts := possibleNxtPts minus {neighbor}; #future edge so remove
                    fi;
                else #we are on the boundary
                    s := sitesOfEdge[1];	
                    sVector := <simplify(fnormal(Re(s - crtPt)),zero),  simplify(fnormal(Im(s - crtPt)),zero), 0>;	
                    sIndex := searchList( s, sites );
                    bool1 := isToRight(directionVector,sVector);
                    if  bool1 then #site is to the right and new 
                        if sIndex > i then # have |new site, this is legal and may be used later but not now
                            possibleNxtPts := possibleNxtPts minus {neighbor}; 
                        fi; # have |old site, this must be the path of interest
                    else #site is to the left
                        if sIndex < i then #old| is legal and are already used
                            possibleNxtPts := possibleNxtPts minus {neighbor}; #old edge
                        fi; #new| is legal and must be the path of interest
                    fi;
                fi;                  
            od;
            #Removed all nieghbors that do not belong, theoretical there should only be one.NEEDS PROOF 
            if nops(possibleNxtPts)=1 then

                nxtPt :=possibleNxtPts[1];          
            else
               print("error: somehow we ended with more than one possible pt");
               nxtPt :=possibleNxtPts[1];                   
               print("deterimining neighbor");
               print(possibleNxtPts);
                mm:=abs(nxtPt-si);
                for pt in (possibleNxtPts minus {nxtPt}) do
                    if evalf(mm-abs(pt-si)) >0. then
                        nxtPt:=pt;
                        mm:=abs(nxtPt-si);
                    fi
                od;       
            fi;         
            p1 := [ op(p1) , [ crtPt , nxtPt]];
            pstPt := crtPt;
            crtPt := nxtPt;
        od; #at this point we have all of the paths needed. We need to glue them together
            if member('bug',{args}) then
            print("permuting p2");
            
            print(p2_permute);
            print(crtPt);
            fi;
            
        while not ( p2[1][1] = crtPt and  p2[-1][2] = crtPt ) do
            mm := p2[1];
            for j from 1 to ArrayNumElems(p2)-1 do  
                p2[j]:=p2[j+1];
            od;
            p2[-1]:=mm;
        od;
        p2:= convert(p2, list);
        p3 := reverseList(map(reverseList,p1));
          if member('bug',{args}) then
           print("path to cell");
           print(p1); 
           print("path on cell");
            print(p2);
           print("path back home");
            print(p3);
            print( plots[display]( 
               plot( map(a->[[Re(a[1]),Im(a[1])],[Re(a[2]),Im(a[2])]],E ), color ="yellow" ), 		
               plot(map(a->[[Re(a[1]),Im(a[1])],[Re(a[2]),Im(a[2])]],[op(p1), op(p2)]), color="red")) );
            fi;   
            
           
        P := [ op(P), [ op(p1), op(p2), op(p3)] ];  
                  
    od;         
    P;
end proc;  
ProduceCycle := proc( cycleSet::set, site, vertexToEdgeData )
    local undirectedCycle := cycleSet, directedCycle:=[], crtPt,nxtPt,possibleEdges, natPt,prevPt;
    #print(undirectedCycle);
	natPt := sortComplex(  [op(op~(undirectedCycle))],'normal' )[1]; #lex ordering is important for our algorithm	
	#this produces a canonical point of the cell
	possibleEdges:=undirectedCycle intersect vertexToEdgeData[natPt]; #all edges that are in this cell that have the canonnical point
	nxtPt:= sortComplex( 
        [op(op~(possibleEdges)minus {natPt})],
        'Arg',natPt)[1]; #there is only two and the one that we want to move to must match orientation. Due to lex, this must be one that does not go left
    #and it must have smallest argument out of two possible relative to cannonical pt.
	directedCycle:=[ op(directedCycle), [natPt,nxtPt] ];#add this edge
	prevPt:=natPt; #update location since we have moved
	crtPt:=natPt;
	while nxtPt <>  natPt do #just take next edge untill we reach the next edge
		prevPt:=crtPt;
		crtPt:=nxtPt;
		nxtPt:=	op( ((vertexToEdgeData[crtPt] intersect undirectedCycle) minus {{crtPt,prevPt}})[1] minus {nxtPt});		
		directedCycle:=[ op(directedCycle), [crtPt,nxtPt] ];
	od;
end proc;

GenerateMonodromy := proc(curve,E,sites,canonicalVertex,preimages, size,P)
    local infBp,edgeToEdgeData, CL:= round( sqrt(max(1, Digits-8))), edge,
	 dy, i, j, k,l, pt1,pt2,fibre, integers, yi,yj,yk,cycle,sigma_inf,sigma,
	  siteToGenerator:=Array([]),siteToMonoData:=Array([]),aux:=[],mm,val,
	  sitesPlus,pathNumber,fibre1,fibre2,fibre2perm;   
    
    for edge in E do
        pt1 := edge[1]; 
        pt2 := edge[2];       
        fibre := [ fsolve(subs(x=pt1 , curve), y,complex)];    
        edgeToEdgeData[[pt1,pt2]]:= [[0,fibre], Continue(CL, curve,x,y,pt1*(1-t)+t*pt2,t,0,1, fibre)]; #this generates vectors that analtically continue the curve at points
            
        edgeToEdgeData[[pt2,pt1]]:= [ seq(
            [
            1- edgeToEdgeData[ [pt1,pt2] ][i][1], 
            edgeToEdgeData[ [pt1,pt2] ][i][2] 
            ], i=1..nops(edgeToEdgeData[[pt1,pt2]]))
            ];

        edgeToEdgeData[[pt2,pt1]] := reverseList(edgeToEdgeData[[pt2,pt1]]);
		#print(edgeToEdgeData[[pt2,pt1]]);    
		#print(edgeToEdgeData[[pt1,pt2]]);     
    od; 
    #analtic continuation against edges is completed, now I need to take the paths, and these analytic continuation data to determine perm.
    for pathNumber from 1 to size do  #for each path
        #print("path1");
        #print(pathNumber);
        #print("of");
        #print(N);
        pt1 :=  P[pathNumber][1][1];  #get location and the preimages associated.
     
        fibre1 := preimages;
        for edge from 1 to nops(P[pathNumber]) do             
            pt2 := P[pathNumber][edge][2]; #move to next point   
            #print("these fibs should=");
            #print(edgeToEdgeData[[pt1,pt2]][1][2] );
            #print(fibre1);
            #print("yeah");
            fibre2:=edgeToEdgeData[[pt1,pt2]][-1][2]; #get the fiber, obtained by traveling to next point
            #print(fibre2);
            #print(edgeToEdgeData[[pt1,pt2]][-1]);
            fibre2perm:=[];
            #now compare the where you were to were you are now fiberwise
            for i from 1 to nops(fibre1) do #i is sent to
                yi:=fibre1[i];
                mm:=infinity;
                for j from 1 to nops(fibre1) do
                    yj:=edgeToEdgeData[[pt1,pt2]][1][2][j];
                    if evalf(abs(yi-yj)-mm)<0 then                
                        mm:=abs(yi-yj);                
                        k:=j;
                    fi
                od; #k
                fibre2perm:=[op(fibre2perm), fibre2[k]]; #[fib_k1,fib_k2,..,fib_kn]  but send starting fiber to this spot via perm   
                          
            od;
            fibre1 := fibre2perm;
            #print(nops(fibre1));
            pt1:=pt2; #now we know and we can move to pt2               
        od;
        #we now have the permuted fibre  
        ArrayTools[Append](siteToGenerator, Perm([]));
        integers:= { seq(i, i=1..nops(fibre1)) };
        i:= min( integers );      #get all the integers and choose the smallest
        while integers <> {} do  
            integers := integers minus {i}; #remeove current integer
            yi := preimages[ i ];    #need to find where the fiber of this int goes to in the permutation        
            mm := infinity;
            for j from 1 to nops(fibre1) do 
                yj:=fibre1[j];             
                if evalf(abs(yi-yj)-mm)<0 then                
                    mm:=abs(yi-yj);                
                    k:=j;
                fi
            od; 
            #i->k where k is the closed fiber element in the permuted ones #add i->k to running permuation. 
            #Either k is known or not, if known the cycle is complete. If not then move to next portion of cycle.
            #code can be cleaned up!
            if k in integers then 
                siteToGenerator[pathNumber] := siteToGenerator[pathNumber] . Perm([[i,k]]); 
                i := k;                           
            else
                i := min(integers);  
            fi;             
        od; 
        #if path has id perm then point is not a branch point, fix this!!!!!
        if siteToGenerator[pathNumber] = Perm([]) then
            #siteToMonoData :=  [op(siteToMonoData), 'NB' ];
        else 
        	ArrayTools[Append](siteToMonoData,[sites[pathNumber],  [] ] );  
            #aux := [op(aux), [sites_ordered[pathNumber],  [] ] ];
            #print(siteToMonoData);
            integers:= {seq(j, j =1..nops(fibre1))};  
            #print(integers);     
            i:=min(integers);
            cycle:=[];
            while not integers={} do    
                integers := integers minus {i};
                cycle:=[op(cycle), i];            
                j:=siteToGenerator[pathNumber][i]; 
                  #i->j         
                if j in integers then 
                    i:=j;                           
                else
                    i:=min(integers);
                    if nops(cycle)>1 then
                        #Sprint(cycle);
                        #print(siteToMonoData[pathNumber]);
                        siteToMonoData[-1][2]:=[op(siteToMonoData[-1][2]), cycle]; 
                       #aux[-1][2]:=[op(aux[-1][2]), cycle];
                    fi;
                    cycle:=[];
                fi;            
            od;
            #siteToMonoData := [op(siteToMonoData), [sites_ordered[i], siteToMonoData[i][2] ] ];
        fi;   
    od;
    #print(siteToMonoData);
    sigma_inf:=Perm([]);
	sitesPlus:=sites;
    for sigma in siteToGenerator do
        sigma_inf:=sigma_inf . sigma;
    od;
    sigma_inf := GroupTheory[PermInverse](sigma_inf);

	 #checking at x=infty
    if sigma_inf <>Perm([]) then
		infBp:=true;
        sitesPlus:=[op(sitesPlus),infinity];
   		ArrayTools[Append](siteToMonoData, [infinity,  [] ] ); 
        #siteToMonoData := [op(siteToMonoData), [infinity,  [] ] ];   
        integers:= {seq(j, j =1..nops(fibre))};       
        i:=min(integers);
        cycle:=[];
        while not integers={} do    
            integers := integers minus {i};
            cycle:=[op(cycle), i];            
            j:=sigma_inf[i];            
            if j in integers then 
                i:=j;                           
            else
                i:=min(integers);
                if nops(cycle)>1 then
                    siteToMonoData[-1][2]:=[op(siteToMonoData[-1][2]), cycle]; 
                fi;
                cycle:=[];
            fi;            
        od;
		check_cstruct(siteToMonoData[-1][2] ,algcurves[puiseux](curve,x=infinity,y,`cycle structure`));
    else infBp:=false; fi;
	#print(infBp);
    infBp, sitesPlus, siteToMonoData, edgeToEdgeData;       
end proc:
check_cstruct:=proc(perm,s)
	local a;
	#print(s);
	#print(perm);
	if perm<>[] and not type(perm[1],list) then
		procname(convert_disjcyc(perm),s)
	else
		a:=map(nops,perm);
		if sort([1$(convert(s,`+`)-convert(a,`+`)),op(a)])<>sort(s) then
			error "Found wrong cycle structure"
		fi
	fi;
	perm
end:

getLine:=proc(tri,slist)
    local si, sj, sk, m1,m2, x01, x02, y01,y02;
    si := [ Re(slist[tri[1]]), Im(slist[tri[1]]) ];
    sj := [ Re(slist[tri[2]]), Im(slist[tri[2]]) ];
    sk := [ Re(slist[tri[3]]), Im(slist[tri[3]]) ];
    if si[2] = sj[2] then
        m1:= simplify(fnormal(evalf(
                -(sk[1]-si[1])/(sk[2]-si[2])
            )),zero);
        y01 := simplify(fnormal(evalf(
                ( sk[2] + si[2] )/2
        )),zero);
        x01 := simplify(fnormal(evalf(
                ( sk[1] + si[1])/2
        ) ),zero);       
        m2 := simplify(fnormal(evalf(
            -(sk[1]-sj[1])/(sk[2]-sj[2])
        )),zero);
        y02 := simplify(fnormal(evalf(
            (sk[2]+sj[2])/2
        )),zero);
        x02:=simplify(fnormal(evalf(
                (sk[1]+sj[1])/2
        )),zero);
    elif si[2]=sk[2] then
        m1:= simplify(fnormal(evalf(
            -(sj[1]-si[1])/(sj[2]-si[2])
        )),zero);
        y01 := simplify(fnormal(evalf(
            ( sj[2] + si[2] )/2
        )),zero);
        x01 := simplify(fnormal(evalf(
            ( sj[1] + si[1])/2
        ) ),zero);
        m2 := simplify(fnormal(evalf(
            -(sj[1]-sk[1])/(sj[2]-sk[2])
        )),zero);
        y02 := simplify(fnormal(evalf(
            (sj[2]+sk[2])/2
        )),zero);
        x02:=simplify(fnormal(evalf(
            (sj[1]+sk[1])/2
        )),zero);
    else
        m1 := simplify(fnormal(evalf(
            -(si[1]-sj[1])/(si[2]-sj[2])
        )),zero);
        y01 := simplify(fnormal(evalf(
            ( si[2] + sj[2] )/2
        )),zero);
        x01 := simplify(fnormal(evalf(
            ( si[1] + sj[1])/2
        ) ),zero);
           
        m2 := simplify(fnormal(evalf(
            -(si[1]-sk[1])/(si[2]-sk[2])
        )),zero);
        y02 := simplify(fnormal(evalf(
            (si[2]+sk[2])/2
        )),zero);
        x02:=simplify(fnormal(evalf(
            (si[1]+sk[1])/2
        )),zero);
    fi;
    return m1*(x-x01)+ y01,m2*(x-x02)+ y02;
end proc;