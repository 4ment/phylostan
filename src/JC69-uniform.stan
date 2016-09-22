
data {
    int <lower=0> L;               // alignment length
    int <lower=0> S;               // number of tips
    real<lower=0,upper=1> tipdata[S,L,4];
    int <lower=0,upper=2*S> peel[S-1,3];   // list of nodes for peeling
}

transformed data {
    int bcount; // number of branches
    int alphalen;
    int root;
    bcount = 2*S-2;
    alphalen=4; // assume DNA
    root = peel[S-1,3]; // root is the final peeling
}

parameters {
	vector <lower=0,upper=10>[bcount] blens; // branch lengths
}


model {
    vector[alphalen] node[2*S,L];   // partial probabilities for the S tips and S-1 internal nodes
    matrix[alphalen,alphalen] fttm[bcount+1]; // finite-time transition matrices for each branch, plus zero-length root branch

    // set some branch length priors
    blens ~ uniform(0,10);

    // calculate finite time transition matrices for each branch
    // under the Jukes-Cantor model
    for( b in 1:bcount ) {
	    for( i in 1:alphalen ) {
        	for( j in 1:alphalen ) {
                fttm[b][i,j] = 0.25 - 0.25*exp(-4*blens[b]/3);
            }
            fttm[b][i,i] = 0.25 + 0.75*exp(-4*blens[b]/3);
        }
    }
    // zero-length root branch identity matrix
    fttm[bcount+1] = rep_matrix(0,alphalen,alphalen);
    for( i in 1:alphalen ) {
        fttm[bcount+1][i,i] = 1;
    }

    // copy tip data into node probability vector
    for( n in 1:S ) {
        for( i in 1:L ) {
            for( a in 1:alphalen ) {
                node[n,i][a] = tipdata[n,i,a];
            }
        }
    }
    
    // calculate tree likelihood for the topology encoded in peel
	for( i in 1:L ) {
        for( n in 1:(S-1) ) {
            node[peel[n,3],i] = (fttm[n*2-1]*node[peel[n,1],i]) .* (fttm[n*2]*node[peel[n,2],i]);
        }
        // multiply by background nt freqs (assuming uniform here)
        node[2*S,i] = node[root,i] / alphalen;

        // add the site log likelihood
        target += log(sum(node[2*S,i]));
    }

}
