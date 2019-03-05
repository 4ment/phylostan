#!/usr/bin/env python


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

		times[map[1,1]] = heights[map[1,1]-S];
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


def heights_to_blens(heterochronous=False):
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
		return model_heights_to_blens.format(' - lowers[map[j,1]]')
	else:
		return model_heights_to_blens.format('')

	
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


def joacobian(heterochronous=False):
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

		for(j in 1:4){
			partials[2*S,i][j] = partials[peel[S-1,3],i][j] * freqs[j];
		}
		// add the site log likelihood
		target += log(sum(partials[2*S,i]))*weights[i];
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


def get_model(substitution='GTR', coalescent=None, heterochronous=True, estimate_rate=False, categories=1, clock=None, invariant=False, **kwargs):
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

	if clock is None:
		transformed_data_declarations.append('int bcount = 2*S-3; // number of branches')
	else:
		data_block.append('int map[2*S-1,2];                     // list of node in preorder [node,parent]')
		transformed_data_declarations.append('int bcount = 2*S-2; // number of branches')

	# Site model
	if invariant or categories > 1:
		model_block_declarations.append('real probs[C];')
		model_block_declarations.append('vector[4] partials[C,2*S,L];   // partial probabilities for the S tips and S-1 internal nodes')
		model_block_declarations.append('matrix[4,4] pmats[bcount*C]; // finite-time transition matrices for each branch')
	else:
		model_block_declarations.append('vector[4] partials[2*S,L];   // partial probabilities for the S tips and S-1 internal nodes')
		model_block_declarations.append('matrix[4,4] pmats[bcount]; // finite-time transition matrices for each branch')
		
	if categories > 1 and not invariant:
		data_block.append('int C;')
		
		parameters_block.append('simplex[C]  ps;')
		parameters_block.append('simplex[C] rate_unscaled;')

		transformed_parameters_declarations.append('vector[C] rs;')
		transformed_parameters_declarations.append('simplex[C] constraint;')

		transformed_parameters_block.append('constraint = ps .* rate_unscaled; // not actually a simplex yet')
		transformed_parameters_block.append('rs = rate_unscaled / sum(constraint);')
		transformed_parameters_block.append('constraint /= sum(constraint); // is now a simplex that equals p .* x ')
	elif invariant and categories == 1:
		transformed_data_declarations.append('int C = 2;')
		
		parameters_block.append('real<lower=0.0, upper=1.0> pinv;')

		transformed_parameters_declarations.append('vector[2] ps;')
		transformed_parameters_declarations.append('vector[2] rs;')

		transformed_parameters_block.append('ps[1] = pinv;')
		transformed_parameters_block.append('ps[2] = 1.0 - pinv;')
		transformed_parameters_block.append('rs[1] = 0.0;')
		transformed_parameters_block.append('rs[2] = 1.0/(1.0 - pinv);')
		
		model_priors.append('pinv ~ uniform(0.0,1.0);')
	elif invariant and categories > 1:
		raise ValueError('Cannot use proportion of invariant and discrete rate heterogeneity yet.')

	# Clock model
	if clock is not None:
		model_block_declarations.append('vector [bcount] blens; // branch lengths')

		transformed_data_declarations.append('int pCount; // number of proportions')
		transformed_data_block.append('pCount = S-2;')

		parameters_block.append('real <lower=0,upper=1> props[pCount]; // proportions')
		transformed_parameters_declarations.append('real <lower=0> heights[S-1];')

		if estimate_rate:
			if clock == 'strict':
				parameters_block.append('real <lower=0> rate;')
				#data_block.append('real <lower=0> meanRate;')
				model_priors.append('rate ~ exponential(1.0/1000);')
		else:
			data_block.append('real <lower=0> rate;')

		functions_block.append(transform_heights(heterochronous))
		if heterochronous:
			data_block.append('real lowers[2*S-1]; // list of lower bounds for each internal node (for reparametrization)')
			data_block.append('real lower_root;')
			parameters_block.append('real <lower=lower_root> height; // root height')
			transformed_parameters_block.append('heights = transform(props, height, map, lowers);')
		else:
			parameters_block.append('real height; // root height')
			transformed_parameters_block.append('heights = transform(props, height, map);')
		
		model_block.append(heights_to_blens(heterochronous))

		# Coalescent
		if coalescent == 'constant':
			functions_block.append(one_on_X)
			parameters_block.append('real <lower=0> theta;')
			model_priors.append('theta ~ oneOnX();')
			functions_block.append(constant_coalescent(heterochronous))
			if heterochronous:
				model_priors.append('heights ~ constant_coalescent(theta, map, lowers);')
			else:
				model_priors.append('heights ~ constant_coalescent(theta, map);')
	else:
		parameters_block.append('vector<lower=0,upper=10> [bcount] blens; // branch lengths')
		model_priors.append('blens ~ exponential(10);')

	# Substitution model
	if substitution == 'GTR':
		data_block.append('vector<lower=0>[4] frequencies_alpha; // parameters of the prior on frequencies')
		data_block.append('vector<lower=0>[6] rates_alpha;       // parameters of the prior on rates')

		parameters_block.append('simplex[6] rates;')
		parameters_block.append('simplex[4] freqs;')

		model_priors.append('rates ~ dirichlet(rates_alpha);')
		model_priors.append('freqs ~ dirichlet(frequencies_alpha);')

		functions_block.append(GTR(categories, invariant))
		if invariant or categories > 1:
			model_block.append('pmats = calculate_gtr_p_matrices(freqs, rates, blens, rs);')
		else:
			model_block.append('pmats = calculate_gtr_p_matrices(freqs, rates, blens);')
	elif substitution == 'JC69':
		model_block_declarations.append('vector[4] freqs = rep_vector(0.25,4);')
		functions_block.append(JC69(categories, invariant))
		if invariant or categories > 1:
			model_block.append('pmats = calculate_jc69_p_matrices(blens, rs);')
		else:
			model_block.append('pmats = calculate_jc69_p_matrices(blens);')
	else:
		raise ValueError('Supports JC69 and GTR only.')

	# Tree likelihood
	model_block.append(likelihood(categories > 1 or invariant, clock is not None))

	if clock is not None:
		model_block.append(joacobian(heterochronous))

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
