
data {
    int <lower=0> L;               // alignment length
    int <lower=0> S;               // number of tips
    real<lower=0,upper=1> tipdata[S,L,4];
    int <lower=0,upper=2*S> peel[S-1,3];   // list of nodes for peeling
    real weights[L];
    vector<lower=0>[4] frequencies_alpha;
    vector<lower=0>[6] rates_alpha;
    vector<lower=0>[3] thetas_alpha;
}

transformed data {
    int bcount; // number of branches
    int root;
    int M; // number of substitution models
    bcount = 2*S-2;
    root = peel[S-1,3]; // root is the final peeling
    M = 3;
}

parameters {
	vector <lower=0,upper=10>[bcount-1] blens; // branch lengths
	// For GTR model
	simplex[6] rates;
	simplex[4] freqs;
	// HKY
	//simplex[4] freqs_hky;
	real <lower=0,upper=10> kappa;

	simplex[M] theta;
}


model {
    vector[4] node[M,2*S,L];   // partial probabilities for the S tips and S-1 internal nodes
    matrix[4,4] fttm[M,bcount]; // finite-time transition matrices for each branch, plus zero-length root branch
    real ps[M];

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

	// For HKY
	real RR;
	real Y;
	real k1;
	real k2;
	real r;
	real exp1;
	real exp22;
	real exp21;

    // set some branch length priors
    blens ~ exponential(20);
    rates ~ dirichlet(rates_alpha);
    freqs ~ dirichlet(frequencies_alpha);
    //freqs_hky ~ dirichlet(frequencies_alpha);
    kappa ~ uniform(0,100);
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

    // HKY model
    RR  = freqs[1] + freqs[3];
	Y  = freqs[4] + freqs[2];
	k1 = kappa * Y + RR;
	k2 = kappa * RR + Y;

	r = 1.0 / (2.0 * (freqs[1] * freqs[2] + freqs[2] * freqs[3] + freqs[1] * freqs[4] + freqs[3] * freqs[4] + kappa * (freqs[2] * freqs[4] + freqs[1] * freqs[3])));

    for( b in 1:bcount-1 ){

		exp1  = exp(-blens[b]*r);
		exp22 = exp(-k2 * blens[b] * r);
		exp21 = exp(-k1 * blens[b] * r);

		//A
		fttm[3,b][1,1]  = freqs[1] * (1. + (Y/RR) * exp1) + (freqs[3]/RR) * exp22; //A
		fttm[3,b][1,2]  = freqs[2] * (1. -         exp1);                        //C
		fttm[3,b][1,3]  = freqs[3] * (1. + (Y/RR) * exp1) - (freqs[3]/RR) * exp22; //G
		fttm[3,b][1,4]  = freqs[4] * (1. -         exp1);                        //T, U

		//C
		fttm[3,b][2,1]  = freqs[1] * (1. -         exp1);                        //A
		fttm[3,b][2,2]  = freqs[2] * (1. + (RR/Y) * exp1) + (freqs[4]/Y) * exp21; //C
		fttm[3,b][2,3]  = freqs[3] * (1. -         exp1);                        //G
		fttm[3,b][2,4]  = freqs[4] * (1. + (RR/Y) * exp1) - (freqs[4]/Y) * exp21; //T, U

		//G
		fttm[3,b][3,1] = freqs[1] * (1. + (Y/RR) * exp1) - (freqs[1]/RR) * exp22; //A
		fttm[3,b][3,2] = freqs[2] * (1. -         exp1);                        //C
		fttm[3,b][3,3] = freqs[3] * (1. + (Y/RR) * exp1) + (freqs[1]/RR) * exp22; //G
		fttm[3,b][3,4] = freqs[4] * (1. -         exp1);                        //T, U

		//T, U
		fttm[3,b][4,1] = freqs[1] * (1. -         exp1);                        //A
		fttm[3,b][4,2] = freqs[2] * (1. + (RR/Y) * exp1) - (freqs[2]/Y) * exp21; //C
		fttm[3,b][4,3] = freqs[3] * (1. -         exp1);                        //G
		fttm[3,b][4,4] = freqs[4] * (1. + (RR/Y) * exp1) + (freqs[2]/Y) * exp21; //T, U
    }

    // zero-length root branch identity matrix
    for( m in 1:M ) {
        fttm[m,bcount] = rep_matrix(0,4,4);
    }
    for( i in 1:4 ) {
        for( m in 1:M ) {
            fttm[m,bcount][i,i] = 1;
        }
    }

    // copy tip data into node probability vector
    for( n in 1:S ) {
        for( i in 1:L ) {
            for( a in 1:4 ) {
                for( m in 1:M ) {
                    node[m,n,i][a] = tipdata[n,i,a];
                }
            }
        }
    }

    // calculate tree likelihood for the topology encoded in peel
	for( i in 1:L ) {

        for( n in 1:(S-1) ) {
            for( m in 1:M){
                node[m,peel[n,3],i] = (fttm[m,n*2-1]*node[m,peel[n,1],i]) .* (fttm[m,n*2]*node[m,peel[n,2],i]);
            }
        }

        // JC69
        node[1,2*S,i] = node[1,root,i] * 0.25;

        // GTR and HKY
        for(j in 1:4){
        	node[2,2*S,i][j] = node[2,peel[S-1,3],i][j] * freqs[j];
        	node[3,2*S,i][j] = node[3,peel[S-1,3],i][j] * freqs[j];
        }

        // add the site log likelihood
        for( m in 1:M){
           ps[m] = log(theta[m]) + log(sum(node[m,2*S,i]));
        }

        target += log_sum_exp(ps)*weights[i];
    }

}
