
data {
    int <lower=0> L;               // alignment length
    int <lower=0> S;               // number of tips
    real<lower=0,upper=1> tipdata[S,L,4];
    int <lower=0,upper=2*S> peel[S-1,3];   // list of nodes for peeling
    real weights[L];
    vector<lower=0>[4] frequencies_alpha;
    vector<lower=0>[6] rates_alpha;
    vector<lower=0>[2] thetas_alpha;
}

transformed data {
    int bcount; // number of branches
    int root;
    bcount = 2*S-2;
    root = peel[S-1,3]; // root is the final peeling
}

parameters {
	vector <lower=0,upper=10>[bcount-1] blens; // branch lengths
	// For GTR model
	simplex[6] rates;
	simplex[4] freqs;
	simplex[2] theta;
}


model {
    vector[4] node[2,2*S,L];   // partial probabilities for the S tips and S-1 internal nodes
    matrix[4,4] fttm[2,bcount]; // finite-time transition matrices for each branch, plus zero-length root branch
    real ps[2];

    // For GTR
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

    // set some branch length priors
    blens ~ exponential(20);
    rates ~ dirichlet(rates_alpha);
    freqs ~ dirichlet(frequencies_alpha);
    theta ~ dirichlet(thetas_alpha);

    // calculate finite time transition matrices for each branch
    // under the Jukes-Cantor model
    for( b in 1:bcount-1 ) {
	    for( i in 1:4 ) {
        	for( j in 1:4 ) {
                fttm[1,b][i,j] = 0.25 - 0.25*exp(-4*blens[b]/3);
            }
            fttm[1,b][i,i] = 0.25 + 0.75*exp(-4*blens[b]/3);
        }
    }
    // GTR model
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

    for( b in 1:bcount-1 ){
		fttm[2,b] = m1 * diag_matrix(exp(eigenvalues*blens[b])) * m2;
    }

    // zero-length root branch identity matrix
    fttm[1,bcount] = rep_matrix(0,4,4);
    fttm[2,bcount] = rep_matrix(0,4,4);
    for( i in 1:4 ) {
        fttm[1,bcount][i,i] = 1;
        fttm[2,bcount][i,i] = 1;
    }

    // copy tip data into node probability vector
    for( n in 1:S ) {
        for( i in 1:L ) {
            for( a in 1:4 ) {
                node[1,n,i][a] = tipdata[n,i,a];
                node[2,n,i][a] = tipdata[n,i,a];
            }
        }
    }

    // calculate tree likelihood for the topology encoded in peel
	for( i in 1:L ) {

        for( n in 1:(S-1) ) {
            for( m in 1:2){
                node[m,peel[n,3],i] = (fttm[m,n*2-1]*node[m,peel[n,1],i]) .* (fttm[m,n*2]*node[m,peel[n,2],i]);
            }
        }

        // multiply by background nt freqs (assuming uniform here)
        node[1,2*S,i] = node[1,root,i] * 0.25;

        for(j in 1:4){
        	node[2,2*S,i][j] = node[2,peel[S-1,3],i][j] * freqs[j];
        }

        // add the site log likelihood
        for( m in 1:2){
           ps[m] = log(theta[m]) + log(sum(node[m,2*S,i]));
        }

        target += log_sum_exp(ps)*weights[i];
    }

}
