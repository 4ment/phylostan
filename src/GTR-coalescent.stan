functions{
	real coalescent_log(real []x, real pop){
		real p;
		real n;
		real l;
		int i;
		real temp;
		n = size(x);
		p = 0;
		l = n + 1.0;
		i = 1;

		while(i <= floor(n)){
			temp = ((l*(l-1.0))/2.0)/pop;
			p = p + log(temp) - (temp * x[i]);
			l = l - 1.0;
			i = i + 1;
		}
		return p;
	}
	
	real oneOnX_log(real x){
		return -log(x);
	}
}

data {
    int <lower=0> L;               // alignment length
    int <lower=0> S;               // number of tips
    real<lower=0,upper=1> tipdata[S,L,4];
    int <lower=0,upper=2*S> peel[S-1,3];   // list of nodes for peeling
    vector<lower=0>[4] frequencies_alpha;
    vector<lower=0>[6] rates_alpha;
    int map[S-1,S];
}

transformed data {
    int bcount; // number of branches
    int tauCount;
    bcount = 2*S-2;
    tauCount = S-1;
}

parameters {
	real <lower=0,upper=100> taus[tauCount]; // intercoalescent times
	simplex[4] freqs;
	simplex[6] rates;
	real <lower=0,upper=100>theta;
}


model {
	vector[4] partials[2*S,L];   // partial probabilities for the S tips and S-1 internal nodes
    matrix[4,4] fttm[bcount]; // finite-time transition matrices for each branch
	matrix[4,4] R; // symmetric rate matrix
	matrix[4,4] Q; // rate matrix
	matrix[4,4] P2; // diagonal sqrt frequencies
	matrix[4,4] P2inv; // diagonal inverse sqrt frequencies
	matrix[4,4] A; // another symmetric matrix
	vector[4] eigenvalues;
	matrix[4,4] eigenvectors;
	matrix[4,4] m1;
	matrix[4,4] m2;
	real s;
	vector [bcount] blens; // branch lengths
	
    // set some branch length priors
    theta ~ oneOnX();
    
    taus ~ coalescent(theta);
    
    freqs ~ dirichlet(frequencies_alpha);
	
	rates ~ dirichlet(rates_alpha);
	//rates ~ gamma(0.05,10);
	

	R[2,1] = rates[1];
	R[1,2] = R[2,1];
	R[3,1] = rates[2];
	R[1,3] = R[3,1];
	R[4,1] = rates[3];
	R[1,4] = R[4,1];
	R[3,2] = rates[4];
	R[2,3] = R[3,2];
	R[4,2] = rates[5];
	R[2,4] = R[4,2];
	R[4,3] = rates[6];
	R[3,4] = R[4,3];
	
	P2 = rep_matrix(0,4,4);
	P2inv = rep_matrix(0,4,4);
	
	for (i in 1:4) {
		R[i,i] = 0.0;
		P2[i,i] = sqrt(freqs[i]);
		P2inv[i,i] = sqrt(1.0/freqs[i]);
	}
	
	Q = R * diag_matrix(freqs);
	s = 0;
	for (i in 1:4) {
		Q[i,i] = 0.0;
		Q[i,i] = -sum(Q[i,1:4]);
		s = s + Q[i,i] * freqs[i];
	}
	Q = -Q/s;
	
	A = P2 * Q * P2inv;

	eigenvalues = eigenvalues_sym(A);
	eigenvectors = eigenvectors_sym(A);
    
    m1 = P2inv * eigenvectors;
    m2 = eigenvectors' * P2;
    
    
    for( b in 1:bcount ){
   		blens[b] = 0;
    }
    
    for( i in 1:tauCount ){
    	for( j in 1:S ){
    		if(map[i,j] > 0){
    			blens[map[i,j]] = blens[map[i,j]] + taus[i];
    		}
    	}
    }
    
    for( b in 1:bcount ){
		fttm[b] = m1 * diag_matrix(exp(eigenvalues*blens[b])) * m2;
    }

    // copy tip data into node probability vector
    for( n in 1:S ) {
        for( i in 1:L ) {
            for( a in 1:4 ) {
                partials[n,i][a] = tipdata[n,i,a];
            }
        }
    }
    
    // calculate tree likelihood for the topology encoded in peel
	for( i in 1:L ) {
		for( n in 1:(S-1) ) {
            //partials[peel[n,3],i] = (fttm[n*2-1]*partials[peel[n,1],i]) .* (fttm[n*2]*partials[peel[n,2],i]);
            partials[peel[n,3],i] = (fttm[peel[n,1]]*partials[peel[n,1],i]) .* (fttm[peel[n,2]]*partials[peel[n,2],i]);
        }
        
        // multiply by background nt freqs (assuming uniform here)
        for(j in 1:4){
        	partials[2*S,i][j] = partials[peel[S-1,3],i][j] * freqs[j];
        }
        // add the site log likelihood
        target += log(sum(partials[2*S,i]));
    }
}
