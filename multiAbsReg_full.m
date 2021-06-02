function [mfpt,Pnode,Qnode,Pedge,Qedge] = multiAbsReg_full(xvals,ka,D,startunif,domrefbound,...
																							getTimetoSoma,vel,options)
% find mean first passage time to absorbance in a region
% reflecting boundary at 0
% reflecting or absorbing boundary at xvals(end)
% default: particle starts at 0, absorbing boundary at xvals(end)
% starting state can be specified as an optional input
% xvals (N-dim vector) divides the region into N edges
% ka (N-dim vector) provides the absorbance within each edge
opt = struct();
% initial particle distribution
opt.startstate = 0;

if(exist('options','var')); opt = copyStruct(options,opt); end
startstate = opt.startstate;


nstate = length(xvals);
if(size(xvals,1)~=1); xvals = xvals'; end
nodelens = diff([0,xvals]); % N-dim vector

if(~any(startstate))
	if(startunif)
		startstate = nodelens./sum(nodelens);
% 		startstate = 1./nstate*ones(1,nstate);
	else
		startstate = zeros(1,nstate);
		startstate(1) = 1;
	end
end

%% set up matrix of transition probabilities
Qedge = zeros(nstate,1);
Pedge = zeros(nstate,nstate);
Qnode = zeros(nstate,1);
Pnode = zeros(nstate,nstate);
mfpt = nan;


% matrices for starting on nodes
for n = 1:nstate
	pq_startunif = false;
	if(n==1) % first node
		L1 = 0;
		L2 = nodelens(1);
		refbounds = [true,(nstate==n)*domrefbound];
		kavals = [0,ka(1)];
		[Q,~,Pp] = partAbsPQ_full(kavals,D,L1,L2,refbounds,pq_startunif);
		if(nstate>1);	Pnode(n,n+1) = Pp; end
		Qnode(1) = Q;
	else
		if(n==2) % node just after first reflecting node
			L1 = nodelens(n-1);
			L2 = nodelens(n);
			refbounds = [true,(nstate==n)*domrefbound];
			kavals = [ka(n-1),ka(n)];
			[Q,~,Pp,Pm] = partAbsPQ_full(kavals,D,L1,L2,refbounds,pq_startunif);
		else % nstate>2
			if(nstate==n && domrefbound)
					L1 = nodelens(n);
					L2 = nodelens(n-1);
					kavals = [ka(n),ka(n-1)];
					refbounds = [true,false];
					[Q,~,Pm,Pp] = partAbsPQ_full(kavals,D,L1,L2,refbounds,pq_startunif);
			else
				L1 = nodelens(n-1);
				L2 = nodelens(n);
				kavals = [ka(n-1),ka(n)];
				refbounds = [false,false];
				[Q,~,Pp,Pm] = partAbsPQ_full(kavals,D,L1,L2,refbounds,pq_startunif);
			end
		end
		if(n<nstate); Pnode(n,n+1) = Pp; end
		Pnode(n,n-1) = Pm;
		Qnode(n) = Q;
	end
end

% time to reach soma
if(getTimetoSoma)
	% get centers of regions
	binedges = xvals(end)-[0,xvals];
	bincenters = (binedges(1:end-1)+binedges(2:end))/2;
	
	% get time to reach soma from bin center
	movetime = bincenters'/vel;
	movetime(end) = 0;
	
	% probabilities of being absorbed at node
	absprob = 1-sum(Pnode,2);
	
	% additional factor for active transport
	Qnode = Qnode+absprob.*movetime;
end

% matrices for starting uniformly
if(startunif)
	pq_startunif = true;
	for n = 1:nstate
		L1 = nodelens(n);
		L2 = 0;
		
		if(n==nstate)
			refbounds = [n==1,domrefbound];
		else
			refbounds = [n==1,false];
		end
		
		kavals = [ka(n),0];
		[Q,~,Pp,Pm] = partAbsPQ_full(kavals,D,L1,L2,refbounds,pq_startunif);
		if(n<nstate); Pedge(n,n+1) = Pp; end
		if(n>1); Pedge(n,n) = Pm; end
		Qedge(n) = Q;
		
	end
	
	% get mfpt starting from an edge
	IPmat = (eye(nstate)-Pnode);
	% starting positon
	IPQ = IPmat\Qnode;
	
% 	mfpt = startstate*(Qedge+Pedge*IPQ);
	mfpt = (startstate*Qedge)+(startstate*Pedge)*IPQ;
	
else
	% get mfpt starting from a node
	IPmat = (eye(nstate)-Pnode);
	% starting positon
	IPQ = IPmat\Qnode;
	
	mfpt = startstate*IPQ;
	
end


end