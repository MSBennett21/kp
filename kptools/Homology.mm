




Homology:=proc(curve,x::name,y::name)
	description "Written by B. Deconinck with small modifications by M. van Hoeij";
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	Canonicalbasis(Digits,args)
end:

# Given a plane representation of a Riemann surface, 
# this procedure gives a canonical basis for the homology 
# of the Riemann surface. 

Canonicalbasis:=proc(DI,curve,x,y)
	local g,hs2,i,j,ttable,tbasis,tmatrix,
	      tlist,alpha,mon,c;
	global C,p,q;
	options remember, 
	    `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	unassign('C','p','q');
	if indets(curve,'name') minus {x,y}<>{} then
		error "Only 2 variables allowed"
	elif indets(curve,float)<>{} then
		error "No floating point coefficients allowed"
	elif args[nargs]<>`give paths` then
		return procname(args,`give paths`)[1]
	fi;
	g:=algcurves['genus'](curve,x,y);
	if g=0 then return `The curve is rational. All cycles are homologous to zero` fi;
	# Compute the monodromy of the curve  = Hurwitz system.
	mon:=monodromy(args[2..nargs]);
	hs2:=nops(mon[2]);
	#if hs2<>nops(GroupTheory['Orbit']('permgroup'(hs2,{seq(i[2],i=mon[3])}),1)) then
	#	error "monodromy not transitive, so curve is reducible"
	#fi;
	hs2:=[hs2,seq([i[1],
		# Add fixed points to the disjoint cycle decomposition.
		sort([op(i[2]),seq(`if`(has(i[2],j),NULL,[j]),j=1..hs2)],
		(l1,l2)->evalb(op(1,l1)<op(1,l2)) )
	],i=mon[3])];
	ttable:=Tretkofftable(hs2);
	tbasis:=Homologybasis(ttable);
	tlist:=Tretkofflist(ttable);
	c:=nops(tlist)/2;
	tmatrix:=Intersectionmatrix(tlist,[seq(p[i],i=1..c),seq(q[i],i=1..c)]);
	if rank(tmatrix)/2 <> g then
		error "Found inconsistent genus"
	fi;
	alpha:=Frobeniustransform(tmatrix,g);

	# Check computation.
	tmatrix := prod(alpha,tmatrix,inplace=false, outputoptions=[]);
	tmatrix := prod(tmatrix,tp(alpha,inplace=false,outputoptions=[]),
		inplace=false,outputoptions=[]);
	for i to c do for j to c do if
	  tmatrix[i,j]<> `if`(j=i+g and i<=g,1,`if`(i=j+g and j<=g,-1,0))
		then error "Frobeniustransform failed"
	fi od od;

	# Place results in a table.
	c:=table();
	c['basepoint']:=mon[1];
	c['sheets']:=mon[2];
	c['genus']:=g;
	c['cycles']:=tbasis;
	c['linearcombination']:=`if`(rowdim(alpha)>2*g,
		delrows(alpha,2*g+1..rowdim(alpha),outputoptions=[]),alpha);
	c['canonicalcycles']:=Canbasis(tbasis,c['linearcombination']);
	[op(c),op(mon[4..-1])]
end:


# This procedure brings any intersection matrix a
# to its canonical form by a transformation 
# alpha * a * transpose(alpha)=b. If 2g=rank(a) and
# d is the size of the square matrix a, then b has
# d-2g null rows and d-2g null columns. These are 
# moved to the lower right corner. On its diagonal,
# b has 2 gxg null blocks. Above the diagonal is a
# gxg identity block. Below the diagonal is a gxg
# -identity block. The output of the procedure is 
# the transformation matrix alpha. 

Frobeniustransform:=proc(a,g)
	local alpha,b,dim,i,j,k,counter,pivot;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	b:=Matrix(a);
	dim:=rowdim(b);
	# The rank of an antisymmetric matrix is always even, it's 2*g.
	alpha := Matrix([1$dim],scan=[diagonal]);
	for i from 1 to g
	# Create the block below the diagonal
	# Make zeros everywhere else in the first
	# g columns.
	do
		# Make sure column i has a suitable pivot.
		counter:=dim;
		while convert(col(delrows(b,1..g+i-1,outputoptions=[]),i,outputoptions=[]),list)=[seq(0,k=1..dim-i-g+1)] do
			rowop(alpha,[i,counter],
				LinearAlgebra:-NoUserValue,inplace=true,outputoptions=[]);
			colop(b,[i,counter],
				LinearAlgebra:-NoUserValue,inplace=true,outputoptions=[]);
			rowop(b,[i,counter],
				LinearAlgebra:-NoUserValue,inplace=true,outputoptions=[]);
			counter:=counter-1;
		od;
		if b[i+g,i]=0 
		 then 
		   # If the pivot element is zero, 
		   # change rows to make it nonzero.
		   k:=i+g+1;
		   while b[i+g,i]=0
		   do
			if b[k,i]<>0 then
			   pivot:=-1/b[k,i];
			   rowop(alpha,k,pivot,
			      inplace=true,outputoptions=[]);
			   rowop(alpha,[k,i+g],
			      LinearAlgebra:-NoUserValue,inplace=true,outputoptions=[]);
			   rowop(b,k,pivot,inplace=true,outputoptions=[]);
			   rowop(b,[k,i+g],LinearAlgebra:-NoUserValue,inplace=true,outputoptions=[]);
			   colop(b,k,pivot,inplace=true,outputoptions=[]);
			   colop(b,[k,i+g],LinearAlgebra:-NoUserValue,inplace=true,outputoptions=[]);
			fi;
			k:=k+1;
		   od;
		 else
		   # Make the pivot element -1.
		   pivot:=-1/b[i+g,i];
		   rowop(alpha,i+g,pivot,inplace=true,outputoptions=[]);
		   rowop(b,i+g,pivot,inplace=true,outputoptions=[]);
		   colop(b,i+g,pivot,inplace=true,outputoptions=[]);
		   fi;
		for j from i to i+g-1 
		# Use the pivot to create zeros in the rows above it.
	    	do
			pivot:=-b[j,i]/b[i+g,i];
			rowop(alpha,[j,i+g],pivot,inplace=true,outputoptions=[]);
			rowop(b,[j,i+g],pivot,inplace=true,outputoptions=[]);
			colop(b,[j,i+g],pivot,inplace=true,outputoptions=[]);
		od;
		for j from i+g+1 by 1 to dim
		# Use the pivot to create zeros in the rows below it.
		do
			pivot:=-b[j,i]/b[i+g,i];
			rowop(alpha,[j,i+g],pivot,inplace=true,outputoptions=[]);
			rowop(b,[j,i+g],pivot,inplace=true,outputoptions=[]);
			colop(b,[j,i+g],pivot,inplace=true,outputoptions=[]);
		od;
	od;
	for i from 1 to g
	# The block above the diagonal is already there
	# Use it to create zeros everywhere else in the
	# second block of g columns. Automatically all
	# other columns are then zero, because the rank
	# of b is only 2g. 
	do
		for j from i+g+1 by 1 to dim
		do
			pivot:=-b[j,i+g];
			rowop(alpha,[j,i],pivot,inplace=true,outputoptions=[]);
			rowop(b,[j,i],pivot,inplace=true,outputoptions=[]);
			colop(b,[j,i],pivot,inplace=true,outputoptions=[]);
		od;
	od;
	# return the transformation matrix alpha.
	alpha
end:

# This procedure generates the intersection matrix b
# starting from a list of p_i's and q_i's. The ordering
# of the p_i's and q_i's determines the entry in the 
# intersection matrix, which is one of -1,0,1. 

Intersectionmatrix:=proc(lijst, elements)
	local a,b,i,j,qi,qj,pj,length,dimen;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	a:=lijst;
	length:=nops(a);
	dimen:=length/2;
	b := Matrix(dimen,dimen,shape=skewsymmetric);
	# The intersection matrix is antisymmetric.
	for i from 1 to dimen-1 do
		a:=Putfirst(a,elements[i]);
		qi:=Whereis(a,elements[i+dimen]);
		for j from i+1 by 1 to dimen do
			pj:=Whereis(a,elements[j]);
			qj:=Whereis(a,elements[j+dimen]);  
			if pj<qi and qi<qj then b[i,j]:=1;
			elif qj<qi and qi<pj then b[i,j]:=-1;
			else b[i,j]:=0
			fi;
		od;
	od;
	b
end:

# This procedure takes a  cyclic list of elements and
# writes it as a linear list with a certain element 
# in the first position.

Putfirst:=proc(L::list,element)
	local pos,i;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	pos:=Whereis(L,`if`(nargs=1,Smallest(L),element));
	[seq(L[i],i=pos..nops(L)),seq(L[i],i=1..pos-1)]
end:

# Given a linear list, determine the position of 
# a given element.

Whereis:=proc(list,element)
	local i;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	for i from 1 to nops(list)+1 do if list[i]=element then
		break
	fi od;
	i
end:

# This procedure in effect encodes from a given Hurwitz 
# system of a Riemann surface a graph. A spanning tree of
# this graph gives a finite set of cycles which contains a
# basis for the homology of the given Riemann surface. The 
# graph is encoded as such a spanning tree. The particular
# way in which this is done allows us to determine the 
# intersection numbers of the cycles. The input of this 
# procedure is a Hurwitz system, the output is a table, 
# referred to as a Tretkoff table. Various other procedures
# are used to get the cycles and their intersection numbers
# from this table. 

Tretkofftable:=proc(hs)
	local tretkofftable,c,disjointcycle,
	      coveringnr,t,i,j,k,branchpoints,mon, 
	      sourcelist,startingindex, startseq,
	      usedbptlabels,usedsheetlabels,
	      finished,level,entry,previous,
	      newlist,a1,a2,a3,a4,b1,b2,b3,b4,label,
	      finaledges,qcounter,pcounter,edge,qedges;
	global C, p, q;

	# Start with an empty table
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	tretkofftable:=table();
	# The number of sheets of the Riemann surface
	coveringnr:=op(1,hs);
	# The number of branchpoints
	t:=nops(hs)-1;
	# A list of the branchpoints
	branchpoints:=[seq(hs[i+1][1],i=1..t)];
	# A list of monodromy permutations associated 
	# with each branchpoint. In disjoint cycle form.
	mon:=[seq(hs[i+1][2],i=1..t)];
	# The branchpoints, together with their permutations
	c:=Array(1..t,1..coveringnr);
	# The pre-images of the branchpoints on the Riemann
	# surface are indexed by the cycles in the disjoint
	# cycle decomposition of the permutation associated to
	# that branchpoint.
	for i from 1 to t do
		for j from 1 to coveringnr do
			disjointcycle:=Findcycle(mon[i],j);
			c[i,j]:=[branchpoints[i],disjointcycle];
		od;
	od;
	finished:=false; 
	startseq:=[];
	# from sheet one, go to all branchpoints on sheet one, 
	# which are not fixed points.
	for i from 1 to t do
		if nops(c[i,1][2])<>1 then
			startseq:=[op(startseq),c[i,1]];
		else startseq:=[op(startseq),[`stop`,c[i,1]]];
		fi;
	od;
	tretkofftable[C[0]]:=[[1,[],startseq,[]]];
	# usedsheetlabels keeps track of which sheets we
	# have already visited.
	usedsheetlabels:={1};
	# usedbptlabels keeps track of which pre-images 
	# of branchpoints have been used already. This included 
	# one-cycles, which are terminal points.
	usedbptlabels:=seq(c[i,1],i=1..t);
	for i from 1 to t do
		for j from 1 to coveringnr do
			if nops(c[i,j][2])=1 then 
				usedbptlabels:=usedbptlabels,c[i,j]
			fi;
		od;
	od;
	usedbptlabels:={usedbptlabels};
	# level denotes how many levels from the root (C[0])
	# we are in the graph. The total number of levels is 
	# finite, but can be arbitrarily large.
	level:=1;
	finaledges:=[];
	qedges:=[];
	# qcounter keeps track of the final points.
	qcounter:=1;
	while not finished do
		finished:=true;
		# previous=number of branches from the 
		# previous level to this level.
		previous:=nops(tretkofftable[C[level-1]]);
		# odd levels correspond to sheets
		if type(level, odd) then 
			# number of vertices at this level equals the
			# number of branches originating at the 
			# previous level and pointing to this level
			tretkofftable[C[level]]:=[seq([],i=1..previous)];
			for i from 1 to previous do
				# it's possible that a vertex at the previous level
				# doesn't point to this level, i.e. it is a final vertex. 
				# In that case, nothing corresponds to that vertex 
				# at this level. 
				if tretkofftable[C[level-1]][i][3]<>[] then
					finished:=false;
					# sourcelist=vertex at the previous level
					sourcelist:=tretkofftable[C[level-1]][i];
					entry:=[seq([],j=1..nops(op(3,sourcelist)))];
					# Construct a new vertex at this level for
					# every branch leaving the vertex at the
					# previous level. 
					for j from 1 to nops(op(3,sourcelist)) do
						# a1=new vertex
						a1:=op(3,sourcelist)[j];
						# a2=originating vertex
						a2:=op(1,sourcelist);
						# branches to future vertices
						a3:=[];
						label:=op(4,sourcelist);
						# a4=label in the graph. Which path to 
						# follow starting from the root to get here.
						a4:=[op(label),j];
						
						if a1[1]<>`stop` 
						then  
						  
						  newlist:=Putfirst(
							    op(3,sourcelist)[j][2],
							    op(1,sourcelist)
								 );
						  for k from 2 by 1 to nops(newlist) do
						
						  if usedsheetlabels intersect {newlist[k]}={}
						  # the sheet has not been visited yet. 
						  then 
							  a3:=[op(a3),newlist[k]];
							  usedsheetlabels:=usedsheetlabels union {newlist[k]};
						  else
						  # the sheet has already been visited.
							  a3:=[op(a3),[`stop`,newlist[k]]];
							  # edge denotes a branch pointing at an endpoint
							  # of the graph. 
							  edge:=[newlist[k],a1];
							  if {op(finaledges)} intersect {edge}={}
							  then
							  # It's the first time this edge occurs: add it
							  # to the list of final edges and mark it by a 
							  # q-endpoint.
							  	finaledges:=[op(finaledges),edge];
							  	qedges:=[op(qedges),edge];
								tretkofftable[q[qcounter]]:=[
									[`stop`,newlist[k]],
									a1,
									[],
									[op(a4),k-1]
											    ];
								qcounter:=qcounter+1;
							  else
							  # This branch has occured before. Find out
							  # where and give it the corresponding p-endpoint.
							  # If is possible, it only occured because of a 
							  # final edge which does not lead to a cycle. 	
							  	if {op(qedges)} intersect {edge}<>{} then
  							  	  pcounter:=Whereis(qedges,edge);
								  tretkofftable[p[pcounter]]:=[
									a1,
									a2,
									a3,
									[op(a4),k-1]
											      ];
								fi;
							  fi;
						  fi;
						
						  od;
						  
						fi;
						
						# create the new vertex.
						entry:=subsop(j=[a1,a2,a3,a4],entry);
					od;
					# add the vertex to the graph. 
					tretkofftable[C[level]]:=subsop(i=entry,tretkofftable[C[level]]);
				fi;	
			od;
		fi;
		# even levels correspond to pre-images 
		# of branchpoints
		if type(level,even) then
			# number of vertices at this level equals the
			# number of branches originating at the 
			# previous level and pointing to this level
			tretkofftable[C[level]]:=[seq([],i=1..previous)];
			for i from 1 to previous do
				# it's possible that a vertex at the previous level
				# doesn't point to this level, i.e. it is a final vertex. 
				# In that case, nothing corresponds to that vertex 
				# at this level. 
				if tretkofftable[C[level-1]][i][3]<>[] then
					finished:=false;
					# sourcelist=vertex at the previous level
					sourcelist:=tretkofftable[C[level-1]][i];
					# Construct a new vertex at this level for
					# every branch leaving the vertex at the
					# previous level. 
					entry:=[seq([],j=1..nops(op(3,sourcelist)))];
					for j from 1 to nops(op(3,sourcelist)) do
						# b1=new vertex
						b1:=op(3,sourcelist)[j];
						# b2=originating vertex
						b2:=op(1,sourcelist);
						# branches to future vertices
						b3:=[];
						# b4=label in the graph. Which path to 
						# follow starting from the root to get here.
						label:=op(4,sourcelist);
						b4:=[op(label),j];
						
						# Note: the order in which the sheets are visited is
						# important. It is obviously given by the monodromy
						# permutation related to each branch point. As a 
						# consequence, the following is split into two parts,
						# which need to be done in order: first the sheets that 
						# are next in the permutation, then the sheets in the
						# permutation preceding the current one.
						
						if b1[1]<>`stop`
						then
						
						  startingindex:=Whereis(branchpoints,b2[1]);
						  for k from startingindex+1 by 1 to t do
						    if usedbptlabels intersect {c[k,b1]}={} 
						    then
						    # the pre-image of the branchpoint
						    # has not been used  
						      if nops(c[k,b1][2])<>1 then
						      	b3:=[op(b3),c[k,b1]];
							usedbptlabels:=usedbptlabels union {c[k,b1]};
						      fi;
						    else
						    # The pre-image of the branchpoint
						    # has been used.
						      b3:=[op(b3),[`stop`,c[k,b1]]];
						      # edge denotes a branch pointing at an endpoint
						      # of the graph. 
						      edge:=[b1,c[k,b1]];
						      if {op(finaledges)} intersect {edge}={}
							  then
							  # It's the first time this edge occurs: add it
							  # to the list of final edges and mark it by a 
							  # q-endpoint.
							  	finaledges:=[op(finaledges),edge];
								if nops(c[k,b1][2])<>1 then 
								  tretkofftable[q[qcounter]]:=[
									[`stop`,c[k,b1]],
									b1,
									[],
									[op(b4),k-startingindex]
									        	      ];
								  qedges:=[op(qedges),edge];
								  qcounter:=qcounter+1;
								fi;
							  else
							  # This branch has occured before. Find out
							  # where and give it the corresponding p-endpoint.
							  # If is possible, it only occured because of a 
							  # final edge which does not lead to a cycle. 	
							  	if {op(qedges)} intersect {edge}<>{} then
							  	  pcounter:=Whereis(qedges,edge);
								  tretkofftable[p[pcounter]]:=[
									b1,
									b2,
									b3,
									[op(b4),k-startingindex]
											      ];
								fi;
							  fi;
						    fi;
						  od; 
						  
						  for k from 1 to startingindex-1 do
						    if usedbptlabels intersect {c[k,b1]}={} 
						    then
						    # the pre-image of the branchpoint
						    # has not been used  
						      if nops(c[k,b1][2])<>1 then
						      	b3:=[op(b3),c[k,b1]];
						      	usedbptlabels:=usedbptlabels union {c[k,b1]};
						      fi;
						    else
						    # The pre-image of the branchpoint
						    # has been used.
						      b3:=[op(b3),[`stop`,c[k,b1]]];
						      # edge denotes a branch pointing at an endpoint
						      # of the graph. 
						      edge:=[b1,c[k,b1]];
						      if {op(finaledges)} intersect {edge}={}
							  then
							  # It's the first time this edge occurs: add it
							  # to the list of final edges and mark it by a 
							  # q-endpoint.
							  	finaledges:=[op(finaledges),edge];
								if nops(c[k,b1][2])<>1 then 
								  tretkofftable[q[qcounter]]:=[
									[`stop`,c[k,b1]],
									b1,
									[],
									[op(b4),k+t-startingindex]
											      ];
								  qedges:=[op(qedges),edge];
								  qcounter:=qcounter+1;
								fi;
							  else
							  # This branch has occured before. Find out
							  # where and give it the corresponding p-endpoint.
							  # If is possible, it only occured because of a 
							  # final edge which does not lead to a cycle. 	
							  	if {op(qedges)} intersect {edge}<>{} then
							  	  pcounter:=Whereis(qedges,edge);
								  tretkofftable[p[pcounter]]:=[
									b1,
									b2,
									b3,
									[op(b4),k+t-startingindex]
											      ];
								fi;
							  fi;
						    fi;
						  od; 
						  
						fi;
						
						# create a new vertex.
						entry:=subsop(j=[b1,b2,b3,b4],entry);
					od;
					# add the new vertex to the graph. 
					tretkofftable[C[level]]:=subsop(i=entry,tretkofftable[C[level]]);			
				fi;
			od;
		fi;
		# Don't bunch the vertices together according to their origin. 
		# All vertices are treated equal. 
		tretkofftable[C[level]]:=[seq(op(tretkofftable[C[level]][i]),i=1..previous)];
		level:=level+1;
	od;
	# How many levels with new information are there?
	# This excludes the last level, which contains only 
	# final points.
	tretkofftable['depth']:=level-2;
	# How many cycles are generated by the spanning tree?
	tretkofftable['numberofcycles']:=qcounter-1;
	op(tretkofftable)
end:

# This procedure returns the cycle which 
# contains a given element of a permutation.

Findcycle:=proc(mon,k)
	local i,j;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	for i in mon do if member(k,i) then return i fi od
end:

# This procedure compares lists. List1 is less
# than List2 if the words formed from these lists,
# word1 is before word2 in the alphabet with 
# letters given by the integers. 

Listordered:=proc(list1,list2)
	local i,j1,j2;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	j1:=nops(list1);
	j2:=nops(list2);
	if j1<=j2 
	then 
		i:=1;
		while i<=j1 do
		  if list1[i]<>list2[i] then 
		  	return evalb( list1[i]<list2[i] )
		  fi;
		  i:=i+1;
		od;
		true
	else
		not procname(list2,list1)
	fi;
end:

# Tretkofflist determines a sequence of pi, 
# qi symbols which are used later to determine
# the intersection indices of the cycles of the
# homology.

Tretkofflist:=proc(ttable)
	local i,j,n,lijst,result;
	global p, q;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	n:=ttable['numberofcycles'];
	result:=NULL;
	lijst:=[seq(ttable[p[i]][4],i=1..n),seq(ttable[q[i]][4],i=1..n)];
	for i in sort(lijst,Listordered) do
		j:=Whereis(lijst,i);
		result:=result,`if`(j>n,q[j-n],p[j])
	od;
	[result]
end:

# This procedure does not really determine
# a basis for the homology. It determines a 
# finite set containing a basis. Some elements
# in the set may be dependent in the homology, 
# however.
# The cycle is found by following the path that 
# leads to qi from the root. Then we follow the 
# path from the root to pi. These paths are pasted 
# together and their overlap around the root is
# removed. 

Homologybasis:=proc(ttable)
	local ppart,qpart,c,i,j,k,m,qi,pi,qpath,ppath,
		vertex,loc1,loc2,part,Z;
	global C, p, q;

	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	c:=table();
	for i from 1 to ttable['numberofcycles'] do
		pi:=ttable[p[i]];
		ppath:=[seq(pi[4][k],k=1..nops(pi[4])-1)];
		qi:=ttable[q[i]];
		qpath:=qi[4];
	    for Z to 2 do
		part:=1;
		k:=0;
		vertex:=ttable[C[0]][1];
		for j in `if`(Z=1,ppath,qpath) do
			part:=part,vertex[3][j];
			loc1:=Whereis(ttable[C[k]],vertex);
			loc2:=0;
			if loc1>1 
			then for m from 1 to loc1-1 do
				loc2:=loc2+nops(ttable[C[k]][m][3])
			     od;
			fi;
			loc2:=loc2+j;
			k:=k+1;
			vertex:=ttable[C[k]][loc2];
		od;
		if Z=1 then
			ppart:=[part]
		else
			qpart:=[part]
		fi
	    od;
		qpart:=subsop(nops(qpart)=qpart[-1][2],qpart);
		c[i]:=Makecycle(ppart,qpart);
	od;	
	c
end:

# This procedure removes the common parts
# of two lists before putting them together
# to create a cycle. 

Makecycle:=proc(a,b)
	local B,i,H;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	H:=0;
	while a[H+1]=b[H+1] do H:=H+1 od;
	B:=[seq(b[i],i=H+1..nops(b)-1)];
	Putfirst([seq(a[i],i=H..nops(a)),seq(B[-i],i=1..nops(B))])
end:


# The cycles of the homology are written with
# their smallest sheetnumber first.

Smallest:=proc(list1)
	local a,i;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	a:=list1;
	if not type(a[1],integer) 
	then 
		a:=[seq(a[i],i=2..nops(a)),a[1]];
	fi;
	min(seq(a[2*i-1],i=1..nops(a)/2))
end:

# This procedure uses the basic cycles obtained 
# in makecycle and the linear combinations which
# are the output of Frobeniustransform to contruct
# the canonical a and b cycles. 

Canbasis:=proc(ttable,alpha)
	local i,g,can ,a,b;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	g:=rowdim(alpha)/2;
	can:=table();
	for i from 1 to 2*g do
		can[`if`(i>g,b[i-g],a[i])]:=Compresscycle(Simplifycycle(
		 ttable,convert(row(alpha,i,outputoptions=[]),list)))
	od;
	op(can)
end:

# This procedure simplifies a cycle which
# is given as a linear combination of other
# cycles. The input is the set of cycles and
# the coefficients of the linear combination.
# The output is a list of cycles. The sum of 
# the cycles in the list constitutes the new
# cycle. 

Simplifycycle:=proc(ttable,lincom)
	local i,j,k,r,cycle,tempcycle,
	      n,index1,sheets,inter,
	      sht,tempsheets;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	r:=nops(lincom);
	i:=1;
	# initialization of the canonical cycle
	while i<>0 do
		if lincom[i]<>0 
		then
			index1:=i;
			if lincom[i]<0 
			then
				# reverse the cycle
				tempcycle:=ttable[i];
				n:=nops(tempcycle);
				tempcycle:=[tempcycle[1],seq(tempcycle[n-j],j=0..n-2)];
				# take as many copies of the cycle as needed
				cycle:=[seq(op(tempcycle),j=1..-lincom[i])];
			else
				tempcycle:=ttable[i];
				cycle:=[seq(op(tempcycle),j=1..lincom[i])];
			fi;
			i:=0;
		else
			i:=i+1;
		fi;
	od;
	# compose with the other cycles that appear
	# in the linear combination.
	sheets:={};
	for i from 1 to nops(cycle)/2 do
		sheets:=sheets union {cycle[2*i-1]};
	od;
	cycle:=[cycle];
	sheets:=[sheets];
	for i from index1+1 to r do
		if lincom[i]<>0 
		then
			if lincom[i]<0 
			then
				tempcycle:=ttable[i];
				n:=nops(tempcycle);
				# reverse the cycle
				tempcycle:=[tempcycle[1],seq(tempcycle[n-j],j=0..n-2)];
				# find a common sheet
				tempsheets:={};
				for j from 1 to nops(tempcycle)/2 do
					tempsheets:=tempsheets union {tempcycle[2*j-1]};
				od;
				# take as many copies of tempcycle as needed
				tempcycle:=[seq(op(tempcycle),j=1..-lincom[i])];
			else
				tempcycle:=ttable[i];
				n:=nops(tempcycle);
				# find a common sheet
				tempsheets:={};
				for j from 1 to nops(tempcycle)/2 do
					tempsheets:=tempsheets union {tempcycle[2*j-1]};
				od;
				# take as many copies of tempcycle as needed
				tempcycle:=[seq(op(tempcycle),j=1..lincom[i])];
			fi;
			k:=1;
			while k<>0 and k<=nops(cycle) do
				inter:=sheets[k] intersect tempsheets;
				if inter<>{} 
				then
					sht:=op(1,inter);
					cycle:=subsop(k=Putfirst(cycle[k],sht),cycle);
					tempcycle:=Putfirst(tempcycle,sht);
					cycle:=subsop(k=[op(cycle[k]),op(tempcycle)],cycle);
					k:=0;
				else k:=k+1;
				     sheets:=[op(sheets),tempsheets];
				fi;
			od;
			if k=nops(cycle)+1
				then cycle:=[op(cycle),tempcycle];
			fi;
		fi;
	od;
	cycle
end:

# This procedure takes a list, 
# which represents a cycle. It compresses 
# the list, which is cyclic. If the i+2'th 
# element in the list is the same as the i'th 
# element, than the i+1'th and i+2'th element 
# can be removed from the list. Finally, the 
# cycle is made to start from the sheet with
# the smallest sheet number. 

Compresscycle:=proc(llist)
	local i,j,k,n,c,cycle;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	n:=nops(llist);
	cycle:=NULL;
	for i from 1 to n do
		c:=llist[i];
		j:=1;
		while j <= nops(c) do
			if j=nops(c) then
				if op(j,c)=op(2,c) then
					c:=[seq(op(k,c),k=3..nops(c))];
					j:=1;
				else j:=j+1;
				fi;
			elif j=nops(c)-1 then
				if op(j,c)=op(1,c) then
					c:=[seq(op(k,c),k=1..nops(c)-2)];
					j:=1;
				else j:=j+1;
				fi;
			else
				if op(j,c)=op(j+2,c) then
					c:=[seq(op(k,c),k=1..j),seq(op(k,c),k=j+3..nops(c))];
					j:=1;
				else j:=j+1;
				fi;
			fi;
		od;
		cycle:=cycle,Reformcycle(Putfirst(c))
	od;
	[cycle]
end:

# This procedure rewrites a cycle in a specific form.
# The output is a list. The odd entries in the list
# are sheet numbers. The even entries consist of lists
# with two elements: The first element is the location 
# of the branchpoint (actually, its location in the
# complex plane), the second element indicates how many 
# times one needs to go around the branchpoint to go to
# the next sheet.

Reformcycle:=proc(cycle)
	local lijst,n,a,b,c,pos1,pos2,around,i,mini;
	option `Copyright (c) 1999 Waterloo Maple Inc. All rights reserved.`;
	n:=nops(cycle);
	lijst:=cycle;
	for i from 1 to n/2 do
		a:=lijst[2*i-1];
		b:=lijst[2*i];
		if 2*i=n then c:=lijst[1];
			 else c:=lijst[2*i+1];
		fi;
		pos1:=Whereis(b[2],a);
		pos2:=Whereis(b[2],c);
		mini:=min(abs(pos2-pos1),pos2-pos1+nops(b[2]),abs(pos2-pos1-nops(b[2])));
		if abs(pos2-pos1)=mini then around:=pos2-pos1;
		elif pos2-pos1+nops(b[2])=mini then around:=pos2-pos1+nops(b[2]); 
		else around:=pos2-pos1-nops(b[2]);
		fi;
		b:=subsop(2=around,b);
		lijst:=subsop(2*i=b,lijst);
	od;
	lijst
end:

