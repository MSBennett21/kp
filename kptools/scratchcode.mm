########process centers of triangles
    #if member('bug', _Env_GlobalOptions ) or  member('bugT', _Env_GlobalOptions )  then
    #    print(plots[display](plotComponent,  
    #        plots[pointplot](complexToPt~({seq(triangleToVert[i][1],i=1..num)}))));          
    #    print(cat(" iteration ", 0, " of the resolution of vertricies"));
    #    print(plots[pointplot](  
    #        complexToPt~({ seq(triangleToVert[i][1],i=1..num)}),symbol=solidcircle,
    #        color = "black" )); 
    #fi;      
    #i:=1; 
    #m :=0;
    #mm1, sameCenters := processVerts( triangleToVert,num );
    #mm := mindist(complexLToPtL([seq(mm1[i][1],i=1..nops(convert(mm1,list)))]))[1];    
    #while  evalf(abs( m-mm ))  <> 0 do
    #    m:=mm;
    #    mm1, sameCenters := processVerts( triangleToVert,num );
    #    mm := mindist(complexLToPtL([seq(mm1[i][1],i=1..nops(convert(mm1,list)))]))[1];   
    #    if member('bug', _Env_GlobalOptions ) or  member('bugT', _Env_GlobalOptions )  then
    #        print(plots[display](plotComponent,  
    #            plots[pointplot](complexToPt~({seq(triangleToVert[i][1],i=1..num)}))));          
    #        print(cat(" iteration ", i, " of the resolution of vertricies"));
    #        print(plots[pointplot](  
    #            complexToPt~({ seq(triangleToVert[i][1],i=1..num)}),symbol=solidcircle,
    #            color = "black" )); 
    #    fi;
    #    i:=i+1;
    #od;
    #infDistVerts := mm;    
    #minDigits := max( minDigits, round(-log10(infDistVerts))+1 ); 
    #if member('bug', _Env_GlobalOptions ) then
    #   print(cat("filtered vetricies ",convert(i, string)," times"));
    #   print(cat("the closest that two points in the set of sites and verts can be to an edge is ",
    #   convert(min(infDistVerts,infDist), string))); 
    #   print(cat("and so the accuracy needed to resovle these points is",convert(minDigits , string)));
    #fi;
processVerts := proc(verts, size)
local i, indexSet := { seq(i, i = 1..size) }, sameCenters:=table(), j,k,m, ci,cj, distinctVerts:=[], center;	
    
    while indexSet <> {} do #while we still have triangles
		i := min(indexSet); #choose the smallest
       	
        ci,ri := verts[i][1],verts[i][2]; # get the center of the triangle: the vert for the edge  
        sameCenters[i]:=[i];
        if member('bugTHARD', _Env_GlobalOptions ) then       
            print(cat("triangle index is ", convert(i,string)," and has center ", convert(ci,string)," and radius ",
            convert(ri, string)));
        fi;     
        indexSet := indexSet minus {i};		
        for j in indexSet do #go through the remaining edges and log the one center that is closest
			cj,rj := verts[j][1],verts[j][2];  
			m :=  abs( ci-cj ); 
            
            if evalf(m) <= min(infDist/2,.1)  then 
                    sameCenters[i]:=[op(sameCenters[i]),j];                    
                    if member('bugTHARD',_Env_GlobalOptions ) then                    
                        print(cat("vertex ", convert(j, string), " is within ", convert( m, string), 
                        " to current vertex. The radius is off by", convert( abs(rj-ri), string)));
                        print(cj);   print(rj);                                     
                    fi;                 
                
            fi;                        
        od;        
        indexSet := indexSet minus { op(sameCenters[i]) } ; 
        if member('bugTHARD',_Env_GlobalOptions ) then 
            print("These triangles have the same center");
            print(sameCenters[i]);            
        fi;
        #center := scrub(Digits add( verts[ sameCenters[i][k] ][1], k=1..nops(sameCenters[i]))/nops(sameCenters[i]));
        #radius := scrub( add( verts[ sameCenters[i][k] ][2], k=1..nops(sameCenters[i]))/nops(sameCenters[i]));
        distinctVerts :=[op(distinctVerts),[center, radius]]; 
        vertexData[ center ]:=table();
        vertexData[ center ]["radius"]:=radius;
        vertexData[ center ]["neighbors"]:={};           
        for k in sameCenters[i] do
            verts[k] := [center,radius];               
        od;         
	od; 
        return distinctVerts, sameCenters;
    end proc;
local i,ti,xi,xs,yi,ys,ysApprox;
		if s<-1 or s>1 then
			error "undefined";
		else
			i:=1;
			ti:=acData[i][1];
			while ti<s do
				i:=i+1;
				ti:=acData[i][1];
			od;
			try
				ti:= ifelse( ti-s < abs(acData[i-1][1]-s), ti,acData[i-1][1] );
				i:= ifelse( ti-s < abs(acData[i-1][1]-s), i,i-1 );		
			end try;
			yi:=acData[i][2];
			xs:=subs(t=s,gamma);
			xi:=subs(t=ti,gamma);
			ysApprox :=[seq(subs(x=xi,y=yi[j], yi[j] - F[4]/F[5]*(xs-xi)),j=1..numSheets)];
			ys := DK( subs(x=xs,F[1]),y, ysApprox, 1/(2*numSheets) );
			if ys = FAIL then
				print([fsolve(subs(x=xs,evalf(F[1]))]);
				error "testrun";
				ys := permuteList([fsolve(subs(x=xs,evalf(F[1]))],ysApprox);
			fi;
	return tp(Matrix([seq(scrubToFloat(Digits,gammaDot*subs(x=xs,y=ys[j],differentialBasis)),j=1..numSheets)]));
	end;