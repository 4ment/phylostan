#!/usr/bin/env python


def autocorrelated_prior(heterochronous):
	str = '''
	real logn_autocorrelated_log(real[] rates, real[] heights, int[,] map, real nu{0}){{
		int S = size(heights) + 1;
		int nodeCount = S + size(heights);
		real logP = 0.0;
		// no rate at root and rate of first child is exponentialy distributed
		for(i in 3:nodeCount){{
			if(map[i,2] == nodeCount){{
				if(map[i,1] > S){{
					logP += normal_lpdf(log(rates[map[i,1]]) | log(rates[map[2,1]]), nu*(2.0*heights[map[i,2]-S] - heights[map[i,1]-S] - heights[map[2,1]-S])/2.0);
				}}
				else{{
					logP += normal_lpdf(log(rates[map[i,1]]) | log(rates[map[2,1]]), nu*(2.0*heights[map[i,2]-S]{1} - heights[map[2,1]-S])/2.0);
				}}
			}}
			else if(map[i,1] > S){{
				logP += normal_lpdf(log(rates[map[i,1]]) | log(rates[map[i,2]]), nu*(heights[map[i,2]-S] - heights[map[i,1]-S])/2.0);
			}}
			else{{
				logP += normal_lpdf(log(rates[map[i,1]]) | log(rates[map[i,2]]), nu*(heights[map[i,2]-S]{1})/2.0);
			}}
		}}
		return logP;
	}}
'''
	if heterochronous:
		return str.format(', real[] lowers', ' - lowers[map[i,1]]')
	else:
		return str.format('', '')


def get_weibull(invariant=False):
	weibull_pinv_site_rates = """
		{
			real m = 0;
			real cat = C - 1;
			real pvar = 1.0 - pinv;
			rs[1] = 0.0;
			ps[1] = pinv;
			for(i in 2:C){
				rs[i] = pow(-log(1.0 - (2.0*(i-2)+1.0)/(2.0*cat)), 1.0/wshape); // weibull inverse cdf with lambda=1
				ps[i] = pvar/cat;
			}
			m = sum(rs)*pvar/cat;
			for(i in 2:C){
				rs[i] /= m;		
			}
		}
"""
	weibull_site_rates = """
		{
			real m = 0;
			for(i in 1:C){
				rs[i] = pow(-log(1.0 - (2.0*(i-1)+1.0)/(2.0*C)), 1.0/wshape);
			}
			m = sum(rs)/C;
			for(i in 1:C){
				rs[i] /= m;		
			}
		}
"""
	if invariant:
		return weibull_pinv_site_rates
	else:
		return weibull_site_rates


def constant_coalescent(heterochronous=False):
	constant_coalescent_str = """
	real constant_coalescent_log(real[] heights, real popSize, int[,] map{0}){{
		int S = size(heights)+1; // number of leaves from the number of internal nodes
		int nodeCount = size(heights) + S;

		real logP = 0.0;
		real lineageCount = 0.0; // first 2 intervals are sampling events

		int indices[nodeCount];
		int childCounts[nodeCount];
		real times[nodeCount];

		real start;
		real finish;
		real interval;
		real logPopSize = log(popSize);

		for( i in 1:nodeCount ){{
			// internal node: transform
			if(map[i,1] > S){{
				times[map[i,1]] = heights[map[i,1]-S];
				childCounts[map[i,1]] = 2;
			}}
			else{{
				times[map[i,1]] = {1};
				childCounts[map[i,1]] = 0;
			}}
		}}

		// calculate intervals
		indices = sort_indices_asc(times);

		// first tip
		start = times[indices[1]];

		for (i in 1:nodeCount) {{
			finish = times[indices[i]];
			
			interval = finish - start;
			if(interval != 0.0){{
				logP -= interval*((lineageCount*(lineageCount-1.0))/2.0)/popSize;
			}}
			
			// sampling event
			if (childCounts[indices[i]] == 0) {{
				lineageCount += 1.0;
			}}
			// coalescent event
			else {{
				lineageCount -= 1.0;
				logP -= logPopSize;
			}}
			
			start = finish;
		}}

		return logP;
	}}
	"""

	if heterochronous:
		return constant_coalescent_str.format(', real[] lowers', 'lowers[map[i,1]]')
	else:
		return constant_coalescent_str.format('', '0')


def skyride_coalescent(heterochronous):
	skyride_coalescent_str = """
	real skyride_coalescent_log(real[] heights, vector pop, int[,] map{0}){{
		int S = size(heights)+1; // number of leaves from the number of internal nodes
		int nodeCount = size(heights) + S;

		real logP = 0.0;
		int index = 1;
		real lineageCount = 0.0;

		int indices[nodeCount];
		int childCounts[nodeCount];
		real times[nodeCount];
		real start;
		real finish;
		real interval;

		for( i in 1:nodeCount ){{
			// internal node: transform
			if(map[i,1] > S){{
				times[i] = heights[map[i,1]-S];
				childCounts[i] = 2;
			}}
			else{{
				times[i] = {1};
				childCounts[i] = 0;
			}}
		}}

		// calculate intervals
		indices = sort_indices_asc(times);

		// first tip
		start = times[indices[1]];

		for (i in 1:nodeCount) {{
			finish = times[indices[i]];

			interval = finish - start;
			// consecutive sampling events
			if(interval != 0.0){{
				logP -= interval*((lineageCount*(lineageCount-1.0))/2.0)/exp(pop[index]);
				if (childCounts[indices[i]] != 0) {{
					logP -= pop[index];
					index += 1;
				}}
			}}

			// sampling event
			if (childCounts[indices[i]] == 0) {{
				lineageCount += 1.0;
			}}
			// coalescent event
			else {{
				lineageCount -= 1.0;
			}}

			start = finish;
		}}

		return logP;
	}}
"""

	if heterochronous:
		return skyride_coalescent_str.format(', real[] lowers', 'lowers[map[i,1]]')
	else:
		return skyride_coalescent_str.format('', '0')


def skygrid_coalescent(heterochronous):
	skygrid_coalescent_str = """
	real skygrid_coalescent_log(real[] heights, vector pop, int[,] map, vector grid{0}){{
		int G = rows(grid);
		int S = size(heights)+1; // number of leaves from the number of internal nodes
		int nodeCount = size(heights) + S;

		real logP = 0.0;
		int index = 1;
		real lineageCount = 0.0;

		int indices[nodeCount];
		int childCounts[nodeCount];
		real times[nodeCount];
		real start;
		real finish;
		real end;

		real popSize;
		real logPopSize;

		for( i in 1:nodeCount ){{
			// internal node: transform
			if(map[i,1] > S){{
				times[map[i,1]] = heights[map[i,1]-S];
				childCounts[map[i,1]] = 2;
			}}
			else{{
				times[map[i,1]] = {1};
				childCounts[map[i,1]] = 0;
			}}
		}}

		// calculate intervals
		indices = sort_indices_asc(times);

		// first tip
		start = times[indices[1]];
		logPopSize = pop[index];
		popSize = exp(logPopSize);

		for (i in 1:nodeCount) {{
			finish = times[indices[i]];

			while(index < G && finish > grid[index]){{
				end = fmin(grid[index], finish);
				logP -= (end - start)*((lineageCount*(lineageCount-1.0))/2.0)/popSize;
				start = end;

				if(index < G){{
					index += 1;
					logPopSize = pop[index];
					popSize = exp(logPopSize);
				}}
			}}
			logP -= (finish - start)*((lineageCount*(lineageCount-1.0))/2.0)/popSize;
			if (childCounts[indices[i]] != 0) {{
				logP -= logPopSize;
			}}

			// sampling event
			if (childCounts[indices[i]] == 0) {{
				lineageCount += 1.0;
			}}
			// coalescent event
			else {{
				lineageCount -= 1.0;
			}}

			start = finish;
		}}
		return logP;
	}}
"""

	if heterochronous:
		return skygrid_coalescent_str.format(', real[] lowers', 'lowers[map[i,1]]')
	else:
		return skygrid_coalescent_str.format('', '0')


def GMRF():
	gmrf_logP = """
	// Not time aware
	real gmrf_log(vector logPopSize, real precision){
		int N = rows(logPopSize);
		real realN = N;
		real s = 0;
		for (i in 2:N){
			s += pow(logPopSize[i-1]-logPopSize[i], 2.0);
		}
		return log(precision)*(realN - 1.0)/2.0 - s*precision/2.0 - (realN - 1.0)/2.0 * log(2.0*pi());
	}
"""
	return gmrf_logP


def heights_to_blens(heterochronous=False, strict=True):
	model_heights_to_blens = """
	// populate blens from heights array in preorder
	for( j in 2:(bcount+1) ){{
		// internal node
		if(map[j,1] > S){{
			blens[map[j,1]] = rate*(heights[map[j,2]-S] - heights[map[j,1]-S]);
		}}
		else{{
			blens[map[j,1]] = rate*(heights[map[j,2]-S]{});
		}}
	}}
"""
	if heterochronous:
		model_heights_to_blens = model_heights_to_blens.format(' - lowers[map[j,1]]')
	else:
		model_heights_to_blens = model_heights_to_blens.format('')
	if not strict:
		return model_heights_to_blens.replace('rate', 'substrates[map[j,1]]')
	return model_heights_to_blens
	
def transform_heights(heterochronous=False):
	transform_str = """
	// transform node heights to proportion, except for the root
	real[] transform(real[] p, real rootHeight, int[,] map{0}){{
		int S = size(p)+2;
		int nodeCount = S*2-1;
		
		real heights[S-1];
		int j = 1;
		
		heights[map[1,1]-S] = rootHeight;
		for( i in 2:nodeCount ){{
			// internal node: transform
			if(map[i,1] > S){{
				heights[map[i,1]-S] = {1}(heights[map[i,2]-S]{2})*p[j];
				j += 1;
			}}
		}}
		return heights;
	}}
	"""

	if heterochronous:
		return transform_str.format(', real[] lowers', 'lowers[map[i,1]] + ', ' - lowers[map[i,1]]')
	else:
		return transform_str.format('', '', '')


def jacobian(heterochronous=False):
	log_det_jacobian = """
	// add log det jacobian
	for( i in 2:(bcount+1) ){{
		// skip leaves
		if(map[i,1] > S ){{
			target += log(heights[map[i,2]-S]{});
		}}
	}}
"""
	if heterochronous:
		return log_det_jacobian.format(' - lowers[map[i,1]]')
	else:
		return log_det_jacobian.format('')


def JC69(C=1, invariant=False):
	jc69_function_str = '''
	matrix[] calculate_jc69_p_matrices(vector blens{0}){{
		{1}
		int bcount = rows(blens);
		matrix[4,4] pmats[bcount{2}]; // probability matrices
		int index = 1;
		real d;
		{3}
			for( b in 1:bcount ) {{
				pmats[index] = rep_matrix(0.25 - 0.25*exp(-blens[b]{4}/0.75), 4, 4);
				d = 0.25 + 0.75*exp(-blens[b]{4}/0.75);
				for( i in 1:4 ) {{
					pmats[index][i,i] = d;
				}}
				index += 1;
			}}
		{5}
		return pmats;
	}}
	'''

	if C > 1 or invariant:
		return jc69_function_str.format(', vector rs', 'int C = rows(rs);', '*C', 'for(c in 1:C){', '*rs[c]', '}')
	else:
		return jc69_function_str.format(*['']*6)


def HKY(C=1, invariant=False):
	hky_function_str = '''
	matrix[] calculate_hky_p_matrices(vector freqs, real kappa, vector blens{0}){{
		{1}
		int bcount = rows(blens);
		matrix[4,4] pmats[bcount{2}]; // probability matrices

		matrix[4,4] Q; // rate matrix
		matrix[4,4] P2 = diag_matrix(sqrt(freqs));        // diagonal sqrt frequencies
		matrix[4,4] P2inv = diag_matrix(1.0 ./ sqrt(freqs)); // diagonal inverse sqrt frequencies
		matrix[4,4] A; // another symmetric matrix
		vector[4] eigenvalues;
		matrix[4,4] eigenvectors;
		matrix[4,4] m1;
		matrix[4,4] m2;
		// symmetric rate matrix
		matrix[4,4] R = [[0.0, 1.0, kappa, 1.0],
						 [1.0, 0.0, 1.0, kappa],
						 [kappa, 1.0, 0.0, 1.0],
						 [1.0, kappa, 1.0, 0.0]];
		real s = 0.0;
		int index = 1;

		Q = R * diag_matrix(freqs);
		for (i in 1:4) {{
			Q[i,i] = 0.0;
			Q[i,i] = -sum(Q[i,1:4]);
			s -= Q[i,i] * freqs[i];
		}}
		Q /= s;

		A = P2 * Q * P2inv;

		eigenvalues = eigenvalues_sym(A);
		eigenvectors = eigenvectors_sym(A);

		m1 = P2inv * eigenvectors;
		m2 = eigenvectors' * P2;

		{3}
			for( b in 1:bcount ){{
				pmats[index] = m1 * diag_matrix(exp(eigenvalues*blens[b]{4})) * m2;
				index += 1;
			}}
		{5}

		return pmats;
	}}
'''

	if C > 1 or invariant:
		return hky_function_str.format(', vector rs', 'int C = rows(rs);', '*C', 'for(c in 1:C){', '*rs[c]', '}')
	else:
		return hky_function_str.format('', '', '', '', '', '')


def GTR(C=1, invariant=False):
	gtr_function_str = '''
	matrix[] calculate_gtr_p_matrices(vector freqs, vector rates, vector blens{0}){{
		{1}
		int bcount = rows(blens);
		matrix[4,4] pmats[bcount{2}]; // probability matrices
		
		matrix[4,4] Q; // rate matrix
		matrix[4,4] P2 = diag_matrix(sqrt(freqs));        // diagonal sqrt frequencies
		matrix[4,4] P2inv = diag_matrix(1.0 ./ sqrt(freqs)); // diagonal inverse sqrt frequencies
		matrix[4,4] A; // another symmetric matrix
		vector[4] eigenvalues;
		matrix[4,4] eigenvectors;
		matrix[4,4] m1;
		matrix[4,4] m2;
		// symmetric rate matrix
		matrix[4,4] R = [[0.0, rates[1], rates[2], rates[3]],
						 [rates[1], 0.0, rates[4], rates[5]],
						 [rates[2], rates[4], 0.0, rates[6]],
						 [rates[3], rates[5], rates[6], 0.0]];
		real s = 0.0;
		int index = 1;

		Q = R * diag_matrix(freqs);
		for (i in 1:4) {{
			Q[i,i] = 0.0;
			Q[i,i] = -sum(Q[i,1:4]);
			s -= Q[i,i] * freqs[i];
		}}
		Q /= s;

		A = P2 * Q * P2inv;

		eigenvalues = eigenvalues_sym(A);
		eigenvectors = eigenvectors_sym(A);

		m1 = P2inv * eigenvectors;
		m2 = eigenvectors' * P2;

		{3}
			for( b in 1:bcount ){{
				pmats[index] = m1 * diag_matrix(exp(eigenvalues*blens[b]{4})) * m2;
				index += 1;
			}}
		{5}

		return pmats;
	}}
	'''

	if C > 1 or invariant:
		return gtr_function_str.format(', vector rs', 'int C = rows(rs);', '*C', 'for(c in 1:C){', '*rs[c]', '}')
	else:
		return gtr_function_str.format('', '', '', '', '', '')


one_on_X = """
	real oneOnX_log(real x){
		return -log(x);
	}
"""

def likelihood(mixture, clock=True):
	init_tip_partials = """
	// copy tip data into node probability vector
	for( n in 1:S ) {
		for( i in 1:L ) {
			for( a in 1:4 ) {
				partials[n,i][a] = tipdata[n,i,a];
			}
		}
	}
"""
	init_tip_partials_mixture="""
	// copy tip data into node probability vector
	for( n in 1:S ) {
		for( i in 1:L ) {
			for( a in 1:4 ) {
				for(c in 1:C){
					partials[c,n,i][a] = tipdata[n,i,a];
				}
			}
		}
	}
"""
	model_calculate_logP = """
	// calculate tree likelihood
	for( i in 1:L ) {
		for( n in 1:(S-1) ) {
			partials[peel[n,3],i] = (pmats[peel[n,1]]*partials[peel[n,1],i]) .* (pmats[peel[n,2]]*partials[peel[n,2],i]);
		}

		for(j in 1:4){
			partials[2*S,i][j] = partials[peel[S-1,3],i][j] * freqs[j];
		}
		// add the site log likelihood
		target += log(sum(partials[2*S,i]))*weights[i];
	}
"""
	model_calculate_mixture_logP = """
	// calculate tree likelihood
	for( i in 1:L ) {
		for( n in 1:(S-1) ) {
			for(c in 1:C){
				partials[c,peel[n,3],i] = (pmats[peel[n,1]+(c-1)*bcount]*partials[c,peel[n,1],i]) .* (pmats[peel[n,2]+(c-1)*bcount]*partials[c,peel[n,2],i]);
			}
		}
		for(c in 1:C){
			probs[c] = ps[c] * sum(partials[c,peel[S-1,3],i] .* freqs);
		}
		// add the site log likelihood
		target += log(sum(probs))*weights[i];
	}
"""
	model_calculate_unconstrained_logP = """
	// calculate tree likelihood
	for( i in 1:L ) {
		for( n in 1:(S-2) ) {
			partials[peel[n,3],i] = (pmats[peel[n,1]]*partials[peel[n,1],i]) .* (pmats[peel[n,2]]*partials[peel[n,2],i]);
		}
		partials[peel[S-1,3],i] = (pmats[peel[S-1,1]]*partials[peel[S-1,1],i]) .* partials[peel[S-1,2],i];

		// add the site log likelihood
		target += log(sum(partials[peel[S-1,3],i] .* freqs))*weights[i];
	}
"""
	model_calculate_unconstrained_mixture_logP = """
	// calculate tree likelihood
	for( i in 1:L ) {
		for( n in 1:(S-2) ) {
			for(c in 1:C){
				partials[c,peel[n,3],i] = (pmats[peel[n,1]+(c-1)*bcount]*partials[c,peel[n,1],i]) .* (pmats[peel[n,2]+(c-1)*bcount]*partials[c,peel[n,2],i]);
			}
		}
		for(c in 1:C){
			partials[c,peel[S-1,3],i] = (pmats[peel[S-1,1]+(c-1)*bcount]*partials[c,peel[S-1,1],i]) .* partials[c,peel[S-1,2],i];
			probs[c] = ps[c] * sum(partials[c,peel[S-1,3],i] .* freqs);
		}
		// add the site log likelihood
		target += log(sum(probs))*weights[i];
	}
"""

	if not mixture:
		model = init_tip_partials
		if clock:
			model += '\n' + model_calculate_logP
		else:
			model += '\n' + model_calculate_unconstrained_logP
	else:
		model = init_tip_partials_mixture
		if clock:
			model += '\n' + model_calculate_mixture_logP
		else:
			model += '\n' + model_calculate_unconstrained_mixture_logP

	return model


def get_model(params):
	functions_block = []

	data_block = []
	transformed_data_declarations = []
	transformed_data_block = []

	parameters_block = []
	transformed_parameters_declarations = []
	transformed_parameters_block = []

	model_block_declarations = []
	model_priors = []
	model_block = []

	data_block.append('int <lower=0> L;                      // alignment length')
	data_block.append('int <lower=0> S;                      // number of tips')
	data_block.append('real<lower=0,upper=1> tipdata[S,L,4]; // alignment as partials')
	data_block.append('int <lower=0,upper=2*S> peel[S-1,3];  // list of nodes for peeling')
	data_block.append('real weights[L];')

	if params.clock is None:
		transformed_data_declarations.append('int bcount = 2*S-3; // number of branches')
	else:
		data_block.append('int map[2*S-1,2];                     // list of node in preorder [node,parent]')
		transformed_data_declarations.append('int bcount = 2*S-2; // number of branches')

	# Site model
	if params.invariant or params.categories > 1:
		data_block.append('int C;')
		model_block_declarations.append('real probs[C];')
		model_block_declarations.append('vector[4] partials[C,2*S,L];   // partial probabilities for the S tips and S-1 internal nodes')
		model_block_declarations.append('matrix[4,4] pmats[bcount*C]; // finite-time transition matrices for each branch')
	else:
		model_block_declarations.append('vector[4] partials[2*S,L];   // partial probabilities for the S tips and S-1 internal nodes')
		model_block_declarations.append('matrix[4,4] pmats[bcount]; // finite-time transition matrices for each branch')

	if params.categories > 1 and params.heterogeneity == 'weibull':
		transformed_parameters_declarations.append('vector[C] ps = rep_vector(1.0/C, C);')

		parameters_block.append('real<lower=0.1> wshape;')
		if params.invariant:
			parameters_block.append('real<lower=0.0, upper=1.0> pinv;')
			model_priors.append('pinv ~ uniform(0.0,1.0);')

		transformed_parameters_declarations.append('vector[C] rs;')
		transformed_parameters_block.append(get_weibull(params.invariant))

		model_priors.append('wshape ~ exponential(1.0);')
	elif params.categories > 1 and not params.invariant:
		parameters_block.append('simplex[C]  ps;')
		parameters_block.append('simplex[C] rate_unscaled;')

		transformed_parameters_declarations.append('vector[C] rs;')
		transformed_parameters_declarations.append('simplex[C] constraint;')

		transformed_parameters_block.append('constraint = ps .* rate_unscaled; // not actually a simplex yet')
		transformed_parameters_block.append('rs = rate_unscaled / sum(constraint);')
		transformed_parameters_block.append('constraint /= sum(constraint); // is now a simplex that equals p .* x ')
	elif params.invariant and params.categories == 1:
		transformed_data_declarations.append('int C = 2;')
		
		parameters_block.append('real<lower=0.0, upper=1.0> pinv;')

		transformed_parameters_declarations.append('vector[2] ps;')
		transformed_parameters_declarations.append('vector[2] rs;')

		transformed_parameters_block.append('ps[1] = pinv;')
		transformed_parameters_block.append('ps[2] = 1.0 - pinv;')
		transformed_parameters_block.append('rs[1] = 0.0;')
		transformed_parameters_block.append('rs[2] = 1.0/(1.0 - pinv);')
		
		model_priors.append('pinv ~ uniform(0.0,1.0);')
	elif params.invariant and params.categories > 1:
		raise ValueError('Cannot use proportion of invariant and discrete rate heterogeneity yet.')

	# Clock model
	if params.clock is not None:
		model_block_declarations.append('vector [bcount] blens; // branch lengths')

		transformed_data_declarations.append('int pCount = S-2; // number of proportions')

		parameters_block.append('real <lower=0,upper=1> props[pCount]; // proportions')
		transformed_parameters_declarations.append('real <lower=0> heights[S-1];')

		if params.estimate_rate:
			if params.clock == 'strict':
				parameters_block.append('real <lower=0> rate;')
				model_priors.append('rate ~ exponential(1000);')
			elif params.clock == 'autocorrelated':
				parameters_block.append('real <lower=0> nu;')
				parameters_block.append('real <lower=0.001, upper=0.01> substrates[2*S-2];')
				functions_block.append(autocorrelated_prior(params.heterochronous))
				if params.heterochronous:
					model_priors.append('substrates ~ logn_autocorrelated(heights, map, nu, lowers);')
				else:
					model_priors.append('substrates ~ logn_autocorrelated(heights, map, nu);')
				model_priors.append('substrates[map[2,1]] ~ exponential(1000);')
			elif params.clock == 'uncorrelated':
				parameters_block.append('real <lower=0> substrates[2*S-2];')
				model_priors.append('substrates ~ exponential(1000);')

		else:
			data_block.append('real <lower=0> rate;')

		functions_block.append(transform_heights(params.heterochronous))
		data_block.append('real lower_root;')
		if params.heterochronous:
			data_block.append('real lowers[2*S-1]; // list of lower bounds for each internal node (for reparametrization)')
			parameters_block.append('real <lower=lower_root> height; // root height')
			transformed_parameters_block.append('heights = transform(props, height, map, lowers);')
		else:
			parameters_block.append('real<lower=lower_root> height; // root height')
			transformed_parameters_block.append('heights = transform(props, height, map);')

		model_block.append(heights_to_blens(params.heterochronous, params.clock == 'strict' or not params.estimate_rate))

		# Coalescent
		if params.coalescent == 'constant':
			functions_block.append(one_on_X)
			parameters_block.append('real <lower=0> theta;')
			model_priors.append('theta ~ oneOnX();')
			functions_block.append(constant_coalescent(params.heterochronous))
			if params.heterochronous:
				model_priors.append('heights ~ constant_coalescent(theta, map, lowers);')
			else:
				model_priors.append('heights ~ constant_coalescent(theta, map);')
		elif params.coalescent == 'skyride':
			functions_block.append(skyride_coalescent(params.heterochronous))
			functions_block.append(GMRF())

			data_block.append('int I; // number of intervals')

			parameters_block.append('vector[I] thetas; // log space')
			parameters_block.append('real<lower=0> tau;')

			if params.heterochronous:
				model_priors.append('heights ~ skyride_coalescent(thetas, map, lowers);')
			else:
				model_priors.append('heights ~ skyride_coalescent(thetas, map);')

			model_priors.append('thetas ~ gmrf(tau);')
			model_priors.append('tau ~ gamma(0.001, 0.001);')
		elif params.coalescent == 'skygrid':
			functions_block.append(skygrid_coalescent(params.heterochronous))
			functions_block.append(GMRF())

			data_block.append('int G; // number of grid interval')
			data_block.append('vector<lower=0>[G] grid;')

			parameters_block.append('vector[G] thetas; // log space')
			parameters_block.append('real<lower=0> tau;')

			if params.heterochronous:
				model_priors.append('heights ~ skygrid_coalescent(thetas, map, grid, lowers);')
			else:
				model_priors.append('heights ~ skygrid_coalescent(thetas, map, grid);')
			model_priors.append('thetas ~ gmrf(tau);')
			model_priors.append('tau ~ gamma(0.001, 0.001);')
	else:
		parameters_block.append('vector<lower=0> [bcount] blens; // branch lengths')
		model_priors.append('blens ~ exponential(10);')

	# Substitution model
	if params.model == 'GTR':
		data_block.append('vector<lower=0>[4] frequencies_alpha; // parameters of the prior on frequencies')
		data_block.append('vector<lower=0>[6] rates_alpha;       // parameters of the prior on rates')

		parameters_block.append('simplex[6] rates;')
		parameters_block.append('simplex[4] freqs;')

		model_priors.append('rates ~ dirichlet(rates_alpha);')
		model_priors.append('freqs ~ dirichlet(frequencies_alpha);')

		functions_block.append(GTR(params.categories, params.invariant))
		if params.invariant or params.categories > 1:
			model_block.append('pmats = calculate_gtr_p_matrices(freqs, rates, blens, rs);')
		else:
			model_block.append('pmats = calculate_gtr_p_matrices(freqs, rates, blens);')
	elif params.model == 'HKY':
		data_block.append('vector<lower=0>[4] frequencies_alpha; // parameters of the prior on frequencies')

		parameters_block.append('real<lower=0> kappa;')
		parameters_block.append('simplex[4] freqs;')

		model_priors.append('kappa ~ lognormal(1.0,1.25);')
		model_priors.append('freqs ~ dirichlet(frequencies_alpha);')

		functions_block.append(HKY(params.categories, params.invariant))
		if params.invariant or params.categories > 1:
			model_block.append('pmats = calculate_hky_p_matrices(freqs, kappa, blens, rs);')
		else:
			model_block.append('pmats = calculate_hky_p_matrices(freqs, kappa, blens);')
	elif params.model == 'JC69':
		transformed_data_declarations.append('vector[4] freqs = rep_vector(0.25,4);')
		functions_block.append(JC69(params.categories, params.invariant))
		if params.invariant or params.categories > 1:
			model_block.append('pmats = calculate_jc69_p_matrices(blens, rs);')
		else:
			model_block.append('pmats = calculate_jc69_p_matrices(blens);')
	else:
		raise ValueError('Supports JC69, HKY and GTR only.')

	# Tree likelihood
	model_block.append(likelihood(params.categories > 1 or params.invariant, params.clock is not None))

	if params.clock is not None:
		model_block.append(jacobian(params.heterochronous))

	script = ''

	if len(functions_block) > 0:
		script += 'functions{' + '\n'.join(functions_block) + '\n}\n\n'

	script += 'data{\n' + '\t' + '\n\t'.join(data_block) + '\n}\n\n'
	
	if len(transformed_data_declarations) != 0:
		script += 'transformed data{\n'
		script += '\t' + '\n\t'.join(transformed_data_declarations) + '\n\n'
		script += '\t' + '\n\t'.join(transformed_data_block) + '\n}\n\n'

	script += 'parameters{\n' + '\t' + '\n\t'.join(parameters_block) + '\n}\n\n'
	
	if len(transformed_parameters_declarations) != 0:
		script += 'transformed parameters{\n'
		script += '\t' + '\n\t'.join(transformed_parameters_declarations) + '\n\n'
		script += '\t' + '\n\t'.join(transformed_parameters_block) + '\n}\n\n'

	script += 'model{\n'
	script += '\t' + '\n\t'.join(model_block_declarations) + '\n\n'
	script += '\t' + '\n\t'.join(model_priors) + '\n\n'
	script += '\t' + '\n\t'.join(model_block) + '\n}\n\n'

	return script


def main():
	print(get_model())


if __name__ == "__main__":
	main()
