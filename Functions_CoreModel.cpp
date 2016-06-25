//
// Created by Gleb Kichaev on 10/21/15.
//


#include <unordered_map>
#include "Functions_CoreModel.h"
#include "Functions_Optimize.h"
#include "Functions_IO.h"

void BuildCausalVector(VectorXd& vec2Build , VectorXd& index){
    vec2Build.setZero();
    for(int i = 0; i <index.size(); i++){
        if(index[i] >0){
            vec2Build[index[i]-1] = 1;
        }

    }
}

void BuildCausalVector(VectorXd& vec2Build , VectorXd& index, vector<int>& index_as_vector){
    vec2Build.setZero();
    for(int i = 0; i <index.size(); i++){
        if(index[i] >0){
            vec2Build[index[i]-1] = 1;
            index_as_vector.push_back(index[i]-1);
        }
    }
}

inline double LogSum(double val1, double val2){
    double logsum = 0;
    if(val1 > val2){
        logsum = log(1 + exp(val2-val1)) + val1;
    }
    else{
        logsum = log(1 + exp(val1-val2)) + val2;
    }

    return(logsum);
}

int NextCombo(VectorXd& c, int k, int n) {
    for (int i= k; --i >= 0;) {
        if (++c[i] <= n-(k-i)) {
            while (++i < k)
                c[i]= c[i-1]+1;
            return 1;
        }
    }
    return 0;
}

double MVN_Density_ZeroMean(VectorXd& x, MatrixXd& sigma){
    MatrixXd b;
    b.setIdentity(sigma.rows(),sigma.cols());
    MatrixXd sigma_inv = sigma.ldlt().solve(b);
    double exponential_term = x.transpose()*(sigma_inv*x);
    double log_det = -.5*log(sigma.determinant());
    return log_det + -.5*exponential_term;
}

int Sample_Int(mt19937 & generator ,vector<double>& weights){
    /* Function that is analagous to sample.int function in R */
    discrete_distribution<int> distribution(weights.begin(), weights.end());
    int sampled_int = distribution(generator);
    return sampled_int;
}

double Get_Rand_Unif( mt19937& generator){
    uniform_real_distribution<double> distribution(0.0,1.0);
    double rand_uniform = distribution(generator);
    return rand_uniform;
}

int Get_Rand_Poiss( mt19937& generator, double poisson_rate){
    poisson_distribution<int> distribution(poisson_rate);
    return distribution(generator);
}

int Next_Move(int set_size, double move_probability, mt19937& generator){
    // #return 0 if swap,  1 if add, 2 if delete
    double rand_uniform= Get_Rand_Unif(generator);
    int next_move;
    int power = set_size-1;
    if(set_size ==2){
        if(rand_uniform < move_probability){
            next_move=1;
        }
        else {
            next_move = 0;
        }
    }
    else{
        if(rand_uniform < pow(move_probability,power)){
            next_move = 1;
        }
        //     else if(rand_uniform < pow((1-move_probability),power)){
        //       next_move = 2;
        //  }
        else{
            next_move = 2;
            // next_move = 2;
        }
    }
    return(next_move);
}
int binomialCoeff(int n, int k)
{
    int C[k+1];
    memset(C, 0, sizeof(C));
 
    C[0] = 1;  // nC0 is 1
 
    for (int i = 1; i <= n; i++)
    {
        // Compute next row of pascal triangle using
        // the previous row
        for (int j = min(i, k); j > 0; j--)
            C[j] = C[j] + C[j-1];
    }
    return C[k];
}

double number_causals(int num_snps, int num_causals){
	int number_of_causals = 0;
	for(int i =1; i <= num_causals; i++){
		number_of_causals += binomialCoeff(num_snps, i);
	}
	return number_of_causals;
} 


int Next_Deletion(VectorXd & zscore_squared, vector<int>& current_set, mt19937& generator){
    //drop a snp from the causal set proportional to 1/zscore^2
    vector<double> probs(current_set.size());
    for(unsigned int i=0; i < current_set.size(); i++){
        probs[i] = 1/zscore_squared[current_set[i]];
    }
    int drop_from_current_set = Sample_Int(generator,probs);
    return drop_from_current_set;
}

int Next_Addition(VectorXd & zscore_squared, MatrixXd& LDsquared_complement, vector<int>& current_set, mt19937& generator){
    //get the next addition based on ld and z_score
    unsigned int N = zscore_squared.size();
    VectorXd ld_weight(N);
    ld_weight.setZero();
    //get max ld to the causal set
    for(unsigned int i =0; i < N; i++){
        for(unsigned int j=0 ; j<current_set.size(); j++){
            if(ld_weight[i]< LDsquared_complement(current_set[j], i)){
                ld_weight[i] = LDsquared_complement(current_set[j], i);
            }
        }
    }
    for(unsigned int j=0 ; j<current_set.size(); j++){
        ld_weight[current_set[j]]=0;
    }
    vector<double> probs(N);
    for(int i=0; i < N; i++){
        probs[i]=zscore_squared[i]*ld_weight[i];
    }
    int add_to_current_set =  Sample_Int(generator, probs);
    return add_to_current_set;
}
double Calculate_LogBayesFactor2(VectorXd& zscore, MatrixXd& ld_matrix, vector<int>& causal_config, double prior_variance){
    unsigned int set_size = causal_config.size();
    MatrixXd ld_c(set_size, set_size);
    VectorXd zscores_c(set_size);
    for(unsigned int i =0; i < set_size; i++){
        zscores_c(i)=zscore(causal_config[i]);
        //cout << "Z s  :" << zscores_c << endl;
        for(unsigned int j = 0; j < set_size; j++){
            ld_c(i,j)= ld_matrix(causal_config[i], causal_config[j]);
        }
    }
    MatrixXd offset = ld_c+prior_variance*(ld_c*ld_c);
    double bf_num = MVN_Density_ZeroMean(zscores_c, offset);
    double bf_denom = MVN_Density_ZeroMean(zscores_c, ld_c);
    /*
       if(isnan(bf_num-bf_denom)){
    //  cout << "The offset matrix is:" << endl;
    //  cout <<offset << endl;
    cout << "num, denom" << endl;
    cout << bf_num << " " << bf_denom << endl;
    cout << "numerator" << endl;
    MVN_Density_ZeroMean_print(zscores_c, offset);
    cout << "denom" << endl;
    MVN_Density_ZeroMean_print(zscores_c, ld_c);
    }
    */  
    return bf_denom-bf_num;
}

double Calculate_LogBayesFactor(VectorXd& zscore, MatrixXd& ld_matrix, vector<int>& causal_config, double prior_variance){
    unsigned int set_size = causal_config.size();
    MatrixXd ld_c(set_size, set_size);
    VectorXd zscores_c(set_size);
    for(unsigned int i =0; i < set_size; i++){
        zscores_c(i)=zscore(causal_config[i]);
        //cout << "Z s  :" << zscores_c << endl;
        for(unsigned int j = 0; j < set_size; j++){
            ld_c(i,j)= ld_matrix(causal_config[i], causal_config[j]);
        }
    }
    MatrixXd offset = ld_c+prior_variance*(ld_c*ld_c);
    double bf_num = MVN_Density_ZeroMean(zscores_c, offset);
    double bf_denom = MVN_Density_ZeroMean(zscores_c, ld_c);
    /*
       if(isnan(bf_num-bf_denom)){
    //  cout << "The offset matrix is:" << endl;
    //  cout <<offset << endl;
    cout << "num, denom" << endl;
    cout << bf_num << " " << bf_denom << endl;
    cout << "numerator" << endl;
    MVN_Density_ZeroMean_print(zscores_c, offset);
    cout << "denom" << endl;
    MVN_Density_ZeroMean_print(zscores_c, ld_c);
    }
    */  
    return bf_num-bf_denom;
}

double Calculate_LogPrior(VectorXd& per_snp_prior,VectorXd& gammas, vector<int> causal_set){
    int N = per_snp_prior.size();
    VectorXd causal_indicator(N);
    causal_indicator.setZero();
    for(unsigned int i=0; i < causal_set.size(); i++){
        causal_indicator[causal_set[i]] = 1;
    }
    double probs = 0;
    for(unsigned int i = 0; i < per_snp_prior.size(); i++){
        if(causal_indicator[i] == 1){
            probs += log(per_snp_prior[i]);
        }
        else{
            probs += log(1-per_snp_prior[i]);
        }
    }
    return(probs);
}

inline double Prior_Snp_Probabilty(VectorXd& beta , VectorXd& aij){
    double dotprod = beta.dot(aij);
    double prob = 1/(1+exp(dotprod));
    return(prob);
}

void Compute_SNP_Priors(MatrixXd& annotation_matrix,VectorXd& gammas, VectorXd& per_snp_prior){
    for(int i= 0; i < annotation_matrix.rows(); i ++){
        VectorXd snp_annotation = annotation_matrix.row(i);
        per_snp_prior[i] = Prior_Snp_Probabilty(gammas, snp_annotation);
    }
}

double Prior_CausalSet_Probabilty(VectorXd& priorJ,VectorXd& beta, VectorXd& C_vector){
    double probs = 0;
    for(int i = 0; i < C_vector.size(); i++){
        if(C_vector[i] == 1){
            probs += log(priorJ[i]);
        }
        else{
            probs += log(1-priorJ[i]);
        }
    }
    return(probs);
}

void CasualSet_To_String(string& causal_string, vector<int>& causal_set){
    causal_string = "";
    for(unsigned int i =0; i < causal_set.size(); i++){
        causal_string += to_string(causal_set[i]);
        causal_string += ",";
    }
    causal_string.pop_back();
}

void Locus_Sampler(VectorXd& marginal, VectorXd& zscores, VectorXd& gammas, MatrixXd& annotations,  MatrixXd& ld_matrix, double& fullLikeli, double prior_variance, double  move_probability, int num_samples, int sampling_seed){
    unsigned int num_snps = zscores.size();
    double log_runsum = -1e150;
    VectorXd per_snp_priors(num_snps);
    random_device rd;
    mt19937 generator(rd());
    generator.seed(sampling_seed);

    Compute_SNP_Priors(annotations, gammas, per_snp_priors);

    unordered_map <string, double> causal_set_posteriors;
    causal_set_posteriors.reserve(num_samples);
    vector<int> causal_set;
    causal_set.push_back(0);
    double log_bayes_factor;
    double log_prior;
    //compute single causal sets
    for(unsigned int i = 0; i < num_snps; i++){
        causal_set[0] = i;
        log_bayes_factor = Calculate_LogBayesFactor(zscores, ld_matrix, causal_set, prior_variance);
        log_prior = Calculate_LogPrior(per_snp_priors, gammas, causal_set);
        marginal[i] = log_bayes_factor+log_prior;
        log_runsum = LogSum(log_runsum, log_bayes_factor+log_prior);
    }
    causal_set[0] = 0;
    causal_set.push_back(num_snps/2);
    MatrixXd ld_sq_complement(num_snps,num_snps);
    VectorXd Zscore_sq(num_snps);
    ld_sq_complement = 1-ld_matrix.array().square();
    Zscore_sq = zscores.array().square();
    int counter = 0;
    //sample multiple causal sets
    int move;
    int index_add;
    int index_delete;
    string causal_string;
    unordered_map<string,double>::const_iterator found;
    while(counter <= num_samples){
        counter++;
        //sample next causal set
        move = Next_Move(causal_set.size(), move_probability,generator);
        if(move == 0){
            index_delete = Next_Deletion(Zscore_sq,causal_set,generator);
            causal_set.erase(causal_set.begin() + index_delete);
            index_add = Next_Addition(Zscore_sq,ld_sq_complement, causal_set,generator);
            causal_set.push_back(index_add);
        }
        else if(move ==1){
            index_add = Next_Addition(Zscore_sq,ld_sq_complement, causal_set,generator);
            causal_set.push_back(index_add);
        }
        else{
            index_delete = Next_Deletion(Zscore_sq,causal_set,generator);
            causal_set.erase(causal_set.begin()+index_delete);
        }
        //check if causal set has already been evaluated and store in hash table if it has not been
        CasualSet_To_String(causal_string, causal_set);
        found = causal_set_posteriors.find (causal_string);
        if ( found == causal_set_posteriors.end()) {
            log_bayes_factor = Calculate_LogBayesFactor(zscores, ld_matrix, causal_set, prior_variance);
            log_prior = Calculate_LogPrior(per_snp_priors, gammas, causal_set);
            causal_set_posteriors[causal_string] = log_bayes_factor+log_prior;
        }
    }
    Marginalize_Sets(causal_set_posteriors, marginal, log_runsum);

    for(int i = 0; i < num_snps; i ++){
        marginal[i] = marginal[i]-log_runsum;
    }

    fullLikeli = fullLikeli+log_runsum;
}

void Generate_Proposal_chi_sq(vector<VectorXd>& zscores, VectorXd& proposal, int num_snps, int num_sets){
    VectorXd divisor(num_snps);
    divisor.setZero();
    divisor = divisor.array()+(num_sets);
    proposal.setZero();
    //get average chi-squares across traits/pops
    for(unsigned int j=0; j < num_sets; j++) {
        proposal = proposal.array() + zscores[j].array().square();
    }
    proposal = proposal.array()/divisor.array();
    double normalizer = proposal.sum();
    //normalize to be a probability
    for(unsigned int i=0; i < num_snps; i++){
        proposal[i] = proposal[i]/normalizer;
    }
}

void Generate_Proposal_Maller(vector<VectorXd>& zscores, VectorXd& proposal, int num_samples){
    VectorXd maller_probs = Zscores2Post(zscores[0]);
    for(unsigned int i = 0; i <zscores.size(); i++){
        proposal[i] = max(maller_probs[i], 3.0/num_samples);
    }

}



void Locus_Importance_Sampler(VectorXd& marginal, vector<VectorXd>& zscores, VectorXd& gammas, MatrixXd& annotations, vector<MatrixXd>& ld_matrix, double& fullLikeli, double prior_variance, int num_samples, int sampling_seed){
    unsigned int num_snps = zscores[0].size();
    unsigned int num_sets = zscores.size();

    random_device rd;
    mt19937 generator(rd());
    generator.seed(sampling_seed);

    VectorXd per_snp_priors(num_snps);
    Compute_SNP_Priors(annotations, gammas, per_snp_priors); //pre-compute priors for each SNP

    vector<vector<int>> sampled_causal_sets;
    vector<double> log_posterior_weight;

    VectorXd proposal(num_snps);
    Generate_Proposal_chi_sq(zscores, proposal, num_snps, num_sets);

    int counter=0;
    vector<int> causal_set;
    double log_bayes_factor;
    double log_prior;
    double log_weight;
    double log_bayes_factor_pop;
    double log_importance_normalizer=-150;
    double log_result;
    double log_import_weight;

    while(counter <= num_samples){
        counter++;
        causal_set.clear();
        Sample_Causal_Vector(proposal,generator, causal_set);
        sampled_causal_sets.push_back(causal_set);
        if ((causal_set.size()>0)) {
            log_bayes_factor=0;
            for (int i = 0;  i<num_sets ; i++) {
                log_bayes_factor_pop = Calculate_LogBayesFactor(zscores[i], ld_matrix[i], causal_set, prior_variance);
                if (std::isnan(log_bayes_factor_pop)) {
                    log_bayes_factor_pop = -1e150;
                    cout << "unstable set" << endl;
                }
                log_bayes_factor += log_bayes_factor_pop;
            }
            log_import_weight = Calc_Importance_logWeight(per_snp_priors, proposal, causal_set);
            log_result = log_bayes_factor + log_import_weight;
            log_importance_normalizer = LogSum(log_importance_normalizer, log_result);

        }
        else{
            log_bayes_factor =-1e150;
            log_import_weight = Calc_Importance_logWeight(per_snp_priors, proposal, causal_set);
            log_result = log_bayes_factor + log_import_weight;
            log_importance_normalizer = LogSum(log_importance_normalizer, log_result);
        }
        log_posterior_weight.push_back(log_result);
    }
    marginal.setZero();
    marginal = marginal.array()+ -1e150;
    Marginalize_Importance_Sets(sampled_causal_sets, log_posterior_weight, log_importance_normalizer, marginal);
    fullLikeli += log_importance_normalizer;

}

void Marginalize_Importance_Sets(vector<vector<int>>& sampled_causal_sets,  vector<double>& log_posterior_weight,  double log_normalizer, VectorXd& marginal_probs){
    vector<int> curr_set;
    int index;
    for(unsigned int i=0; i < sampled_causal_sets.size(); i++){
        curr_set = sampled_causal_sets[i];
        if(curr_set.size()>0){
            for(int j=0; j<curr_set.size(); j++){
                index = curr_set[j];
                marginal_probs[index] = LogSum(marginal_probs[index], log_posterior_weight[i]-log_normalizer);
            }
        }
    }
}


void Sample_Causal_Vector(VectorXd& causal_prior, mt19937& generator, vector<int>& causal_set){
    for(int i =0; i < causal_prior.size(); i++){
        bernoulli_distribution distribution(causal_prior[i]);
        if(distribution(generator)) {
            causal_set.push_back(i);
        }
    }
}


double Calc_Importance_logWeight(VectorXd& per_snp_prior, VectorXd& proposal_prior, vector<int>& causal_set){
    double log_prior = 0;
    double log_proposal = 0;
    for(int i = 0; i < per_snp_prior.size(); i++){
        if(find(causal_set.begin(), causal_set.end(), i) != causal_set.end()){
            log_prior += log(per_snp_prior[i]);
            log_proposal += log(proposal_prior[i]);
        }
        else{
            log_prior += log(1-per_snp_prior[i]);
            log_proposal += log(1-proposal_prior[i]);
        }
    }
    return(log_prior-log_proposal);
}


void string_to_int(string& input_string, vector<int>& output_int){
    vector<string> string_split = split(input_string, ',');
    for(int i = 0; i < string_split.size(); i++){
        output_int.push_back(stoi(string_split[i]));
    }
}

void Marginalize_Sets(unordered_map<string,double>& causal_posteriors, VectorXd& marginals, double& running_sum){
    string string_set;
    double posterior;
    int counter =0;
    vector<int> num_k_sampled;
    for(int k = 0; k < 6; k++){
        num_k_sampled.push_back(0);
    }
    for( auto it = causal_posteriors.begin(); it != causal_posteriors.end(); ++it){
        counter++;
        string_set =  it->first;
        posterior =  it->second;
        vector<int> causal_set;
        string_to_int(string_set, causal_set);

        int causal_set_size = causal_set.size();
        for(int j =0; j<causal_set_size; j++){
            marginals[causal_set[j]] = LogSum(marginals[causal_set[j]], posterior);
        }
        if(causal_set_size> 5){
            num_k_sampled[5]++;
        }
        else{
            num_k_sampled[causal_set_size]++;
        }
        running_sum = LogSum(running_sum, posterior);
    }

    for(int k = 2; k < 5; k++){
        cout << "Sampled " << num_k_sampled[k] << " unique multi-causal sets of size: " << k << endl;
    }
    cout << "Sampled " << num_k_sampled[5] << " unique multi-causal sets of size 5 or greater " << endl;
    cout << "Sampled " << counter << " unique multi-causal sets in total" << endl << endl;
}

void Enumerate_Posterior2(VectorXd& Marginal, vector<VectorXd>& zscores, VectorXd& beta, MatrixXd& annotations,  vector<MatrixXd>& ld_matrix, int NC, double& fullLikeli, double prior_variance, VectorXd& Marginal1, VectorXd& Marginal2, VectorXd& Marginal3, VectorXd &gamma2){
    int numsnps = zscores[0].size();
    int num_pops = zscores.size();
    double runsum = -1e150;
    for(int i =0 ; i < Marginal.size(); i++){
        Marginal(i) = -1e150;
        Marginal1(i) = -1e150;
        Marginal2(i) = -1e150;
        Marginal3(i) = -1e150;
    }
    double runsumone = -1e150;
    double runsumtwo = -1e150;
    double runsumthree = -1e150;
    double c_prob = 0;
    double sum = 0;

    VectorXd mark_marginal(NC);

    VectorXd causConfig(numsnps);
    causConfig.setZero();
    VectorXd causal_index(NC);
    causal_index.setZero();
    VectorXd per_snp_prior(numsnps);
    //cout << "BETA = ";
    for(int i=0; i < beta.size(); i++){
        cout << beta[i] << endl;
    }
    VectorXd m1(numsnps);
    VectorXd m2(numsnps);
    VectorXd m3(numsnps); 
    double log_bayes_factor = 0;
    double log_bayes_factor_pop;
    double log_bayse_factor_pop_one;
    double log_bayse_factor_pop_two;
    double sumone;
    double sumtwo;
    double sumthree;
    double directional_bayes_one;
    //runsum = sum;
    int counter = 1;
    //cout << "Enumaret again" << endl;
    VectorXd zscores1 = zscores[0];
    VectorXd zscores2 = zscores[1];
    VectorXd all_snps(numsnps);
    VectorXd running_sums(numsnps);
    cout << "What gammas are we running with " << endl;
    cout << gamma2 << endl;
    double c_prob1, c_prob2, c_prob3;
    VectorXd test(1);
    test[0] = log(1/(gamma2[0])- 1); 
    Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
    c_prob1 = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
    test[0] = log(1/(gamma2[1])- 1); 
    Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
    c_prob2= Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
    //cout << "log_bayse_factor_pop_one " << endl;
    test[0] = log(1/(gamma2[2])- 1); 
    Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
    c_prob3 = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
    sumtwo = log_bayse_factor_pop_two;
	int number_of_causals_total = number_causals(numsnps, NC);
	cout << "Number of causals total " << number_of_causals_total << endl;
    MatrixXd causal_sets(numsnps, number_of_causals_total); 
	cout << causal_sets.rows()  << " " << causal_sets.cols() << endl;
    VectorXd sumone_vec(number_of_causals_total);
    VectorXd sumtwo_vec(number_of_causals_total);
    while(NextCombo(causal_index, NC, numsnps+1) == 1){
        VectorXd test(1);
        test[0] = gamma2[0];
        Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
        c_prob1 = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
        test[0] = log(1/(gamma2[1])- 1); 
        Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
        c_prob2= Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
        //cout << "log_bayse_factor_pop_one " << endl;
        test[0] = log(1/(gamma2[2])- 1); 
        Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
        c_prob3 = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
        //out << " SNP number " << causal_index << endl;
        vector<int> causal_index_as_vector;
        BuildCausalVector(causConfig, causal_index,causal_index_as_vector);
        //cout << "Beta= " << beta << endl;
        test[0] = log(1/(gamma2[0])- 1); 
        //cout << "TEST " << test[0] << endl;
        Compute_SNP_Priors(annotations, test, per_snp_prior);
        c_prob = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
        //cout << c_prob << endl;
        //cout << c_prob << endl;
        log_bayes_factor = 0;
        log_bayse_factor_pop_one = 0;
        log_bayse_factor_pop_two = 0;
        // Let's try and do som region calculations. 
        for (int i = 0;  i<num_pops ; i++) {
            //cout << causal_index_as_vector[0] << endl;
            log_bayes_factor_pop=Calculate_LogBayesFactor(zscores[i], ld_matrix[i], causal_index_as_vector, prior_variance);
            //cout << "Log bayes " << log_bayes_factor_pop << " " << log_bayes_factor_not_causal << endl;
            if(i==0){
                log_bayse_factor_pop_one = log_bayes_factor_pop;
            }else{
                log_bayse_factor_pop_two = log_bayes_factor_pop;
            }
            //  cout << "Colocalisation " << log_bayes_factor << endl;
            if(std::isnan(log_bayes_factor_pop)){
                cout << "Nan" << endl;
                log_bayes_factor_pop = -1e150;
            }
            log_bayes_factor += log_bayes_factor_pop;
        }
        //cout << "Causal index " << causal_index.size() << " " << causal_index[0] << " " << causal_index[1] << " " << causal_index[2] << endl;
        VectorXd zscores1 = zscores[0];
        VectorXd zscores2 = zscores[1];
        int causal_idx = causal_index[causal_index.size()-1];
        double tmp_c_prob3 =0.0;
        //cout << "Sum " << endl;
        //cout << sum << endl; 
        //cout << "Gamma 2 " << gamma2[1];
        test[0] = log(1/(gamma2[1])- 1);
        Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
        //cout << per_snp_prior << endl;
        c_prob = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
        sumone = log_bayse_factor_pop_one;//
         //   + c_prob; 
        //cout << "log_bayse_factor_pop_one " << endl;
        test[0] = log(1/(gamma2[2])- 1); 
        Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
        c_prob = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
        sumtwo = log_bayse_factor_pop_two;
        //cout << "Causal IDX " << causal_idx << endl; 
        //cout << "Gamma 2 " << gamma2[1] << endl;
        //cout << "log bayes " << endl;
        //cout << (c_prob + c_prob) << endl;
        sum = log_bayes_factor; 
        sumone_vec[counter-1] = sumone; 
        sumtwo_vec[counter-1] = sumtwo; 
        for(int i=0; i < causConfig.size(); i++){
            causal_sets(i, counter -1 ) = causConfig[i]; 
        }
//        for (int i=0; i < causal_index.size(); i++){
//            causal_idx = causal_index[i];
//            if(causal_index[i] > 0){
//                if((zscores1(causal_idx-1) >  0 && zscores2(causal_idx-1)>0)|| (zscores1(causal_idx-1)< 0 && zscores2(causal_idx-1)< 0)){
//                    sumone_vec[i]= c_prob+log_bayes_factor;
//                    //cout << c_prob << " " << c_prob+log_bayes_factor << endl;
//                    sumtwo_vec[i] = 0;
//                    mark_marginal[i] = 1;
//                }else{
//                    sumone_vec[i]= 0;
//                    sumtwo_vec[i] = c_prob+log_bayes_factor;
//                    mark_marginal[i] = 0;
//                }
//            }
//        }
        //cout << "Zscores" << endl;
        //cout << zscores2(causal_idx-1) << endl;
        //cout << zscores1(causal_idx-1) << endl;
        //cout << "TEST" << endl;
            //+c_prob;
        
        test[0] = gamma2[2];
        Compute_SNP_Priors(annotations, test, per_snp_prior); //pre-compute priors for each SNP
        c_prob = Prior_CausalSet_Probabilty(per_snp_prior, test, causConfig);
        //cout << "C _prob " << c_prob << endl;
        ////cout << "Sum " << exp(sum) << endl;
        ////cout << "Causal INDEX " << causal_index << endl; 
        //cout << exp(c_prob) << endl;
        //cout << "Prior odds * log bayse " << (sum) << endl;
        //cout << "Prior odds * log bayse pop one " << (sumone) << endl;
        //cout << "Prior odds * log bayse pop two" << (sumtwo) << endl;
        ////cout << "Running summary" << endl;
        //cout << exp(c_prob)*exp(log_bayse_factor_pop_two) << endl;
        runsum = LogSum(runsum, sumone);
        runsum = LogSum(runsum, sumtwo);
        runsumone = LogSum(runsumone, sumone); 
        runsumtwo = LogSum(runsumtwo, sumtwo);
      //   for(int j = 0; j < causal_index.size(); j++){
      //    //   //cout << "Causal index bleg " << causal_index[j] -1 << endl;
      //           if(causal_index[j] >0){
      //           Marginal[causal_index[j]-1] = LogSum(Marginal[causal_index[j]-1], sum); 
      //           Marginal1[causal_index[j]-1] = LogSum(Marginal1[causal_index[j]-1], sumone);
      //           Marginal2[causal_index[j]-1] = LogSum(Marginal2[causal_index[j]-1], sumtwo);
      //           }
      //       }
       // //cout << "Causal probability " << c_prob << endl;
     //   cout << sum << endl;
        counter ++;
    }
    double four_probs; 
    double tmp_c_prob2;
    double llk = -1e-150;
    for(int i = 0; i < causal_sets.cols(); i++){
        cout << "I " << i << endl;
        for(int j =0; j < causal_sets.cols(); j++){
            double tmp_c_prob = 0.0;
            for(int k = 0; k < causal_sets.rows(); k++){
                if(causal_sets(k,j) == 1 && causal_sets(k,i) ==1){
                   tmp_c_prob += log(gamma2[0]);
                }else if(causal_sets(k, j) == 1 && causal_sets(k, i) == 0){
                   tmp_c_prob += log(gamma2[1]);
                }else if(causal_sets(k, j) == 0 && causal_sets(k, i) == 1){
                    tmp_c_prob += log(gamma2[2]);
                }else{
                    tmp_c_prob += log(1-gamma2[0] - gamma2[1] - gamma2[2]);
                }

            }
            four_probs =  tmp_c_prob + sumone_vec[i] + sumtwo_vec[j]; 
            for(int k = 0; k < causal_sets.rows(); k++){
                 if(causal_sets(k,j) == 1 && causal_sets(k,i) ==1){ 
                        Marginal[k] = LogSum(Marginal[k], four_probs); 
                  }else if(causal_sets(k,j) == 1 && causal_sets(k,i) == 0){
                        Marginal1[k] = LogSum(Marginal1[k], four_probs);
                  }else if(causal_sets(k,j) == 0 && causal_sets(k,i) == 1){
                        Marginal2[k] = LogSum(Marginal2[k], four_probs);
                   }
            }
            //cout << "C1= " << i << " " << " C2= " << j << " " <<(tmp_c_prob2) << " " << marginals(j,i) << " " << marginals(i,j) << endl;
            //cout << "C3 " << (four_probs) << endl; 
            llk = LogSum(llk, four_probs);
        }
    }
    for(int f = 0 ; f < Marginal.size(); f++){
        Marginal[f] = Marginal[f] - llk; 
        Marginal1[f] = Marginal1[f] - llk; 
        Marginal2[f] = Marginal2[f] - llk;
        Marginal3[f] = log(1-exp(Marginal[f]) - exp(Marginal1[f]) - exp(Marginal2[f])); 
    }
    fullLikeli = fullLikeli+llk;
       cout << "Full likeli " << fullLikeli << endl; 
}

void Enumerate_Posterior(VectorXd& Marginal, vector<VectorXd>& zscores, VectorXd& beta, MatrixXd& annotations,  vector<MatrixXd>& ld_matrix, int NC, double& fullLikeli, double prior_variance){
    int numsnps = zscores[0].size();
    int num_pops = zscores.size();
    double runsum = -1e150;
    for(int i =0 ; i < Marginal.size(); i++){
        Marginal(i) = -1e150;
    }
    double c_prob = 0;
    double sum = 0;
    VectorXd causConfig(numsnps);
    causConfig.setZero();
    VectorXd causal_index(NC);
    causal_index.setZero();
    VectorXd per_snp_prior(numsnps);
    Compute_SNP_Priors(annotations, beta, per_snp_prior);

    double log_bayes_factor = 0;
    double log_bayes_factor_pop;
    //runsum = sum;
    int counter = 1;
    //cout << "Enumaret again" << endl;
    while(NextCombo(causal_index, NC, numsnps+1) == 1){
        //cout << "Causal index " << causal_index << endl;
        vector<int> causal_index_as_vector;
        BuildCausalVector(causConfig, causal_index,causal_index_as_vector);
        for (int i = 0 ; i < causal_index_as_vector.size(); i++){
            //cout << "test" << endl;
            //cout << causal_index_as_vector[i] << " ";
        }
        //cout << endl;
        c_prob = Prior_CausalSet_Probabilty(per_snp_prior, beta, causConfig);
        log_bayes_factor = 0;

        // Let's try and do som region calculations. 
        for (int i = 0;  i<num_pops ; i++) {
            log_bayes_factor_pop=Calculate_LogBayesFactor(zscores[i], ld_matrix[i], causal_index_as_vector, prior_variance);
            //cout << "Colocalisation " << log_bayes_factor << endl;
            if(std::isnan(log_bayes_factor_pop)){
                log_bayes_factor_pop = -1e150;
            }
            log_bayes_factor += log_bayes_factor_pop;
        }

        sum = c_prob+log_bayes_factor;
        //cout << "Sum " << exp(sum) << endl;
        //cout << "Prior odds * log bayse " << exp(c_prob)*exp(log_bayes_factor) << endl;
        runsum = LogSum(runsum, sum);
        for(int j = 0; j < causal_index.size(); j++){
            //cout << "Causal index bleg " << causal_index[j] -1 << endl;
            if(causal_index[j] >0){
                Marginal[causal_index[j]-1] = LogSum(Marginal[causal_index[j]-1], sum);
            }
        }
        //cout << "Causal probability " << c_prob << endl;
        //cout << sum << endl;
        counter ++;
    }
    for(int f = 0 ; f < Marginal.size(); f++){
        Marginal[f] = Marginal[f]- runsum;
    }
    
    fullLikeli = fullLikeli+runsum;
}


double Estep2(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance,int  max_causals, VectorXd &gamma_intitial2){
    //cout << "E-STEP" << endl;
    vector<VectorXd> marginal_i;
    vector<VectorXd> marginal_j;
    vector<VectorXd> marginal_k;
    vector<vector <double>> test_data = E_out.prob_vectors;
    cout << test_data[0].size() << endl;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;
    for(unsigned int i = 0; i < Zscores.size(); i ++){
        VectorXd locus_log_marginals(Zscores[i][0].size());
        VectorXd locus_log_marginals1(Zscores[i][0].size());
        VectorXd locus_log_marginals2(Zscores[i][0].size());
        VectorXd locus_log_marginals3(Zscores[i][0].size());
        Enumerate_Posterior2(locus_log_marginals,Zscores[i], betas, Aijs[i], ld_matrices[i], max_causals, fullLikeli, prior_variance, locus_log_marginals1, locus_log_marginals2, locus_log_marginals3, gamma_intitial2);
        exp_temp = locus_log_marginals.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        for(int j = 0; j < exp_temp.size(); j++){
            test_data[0].push_back(stack_temp[j]);
        }
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
        // Add new marginal likelihoods
        exp_temp = locus_log_marginals1.array().exp();
        marginal_j.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        for(int j = 0; j < exp_temp.size(); j++){
            test_data[1].push_back(stack_temp[j]);
        }
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
        
        exp_temp = locus_log_marginals2.array().exp();
        marginal_k.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        for(int j = 0; j < exp_temp.size(); j++){
            test_data[2].push_back(stack_temp[j]);
            cout << test_data[2][j] << endl;
        }
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
        exp_temp = locus_log_marginals3.array().exp();
        marginal_k.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        for(int j = 0; j < exp_temp.size(); j++){
            test_data[3].push_back(stack_temp[j]);
            cout << test_data[2][j] << endl;
        }
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }
    E_out.probs_locs = marginal_i;
    E_out.probs_locs1 = marginal_j;
    E_out.probs_locs2 = marginal_k;
    E_out.probs_stacked = stacker;
    E_out.prob_vectors = test_data;
    return(fullLikeli);
}

//E step with full enumeration
double Estep(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance,int  max_causals){
    vector<VectorXd> marginal_i;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;
    for(unsigned int i = 0; i < Zscores.size(); i ++){
        VectorXd locus_log_marginals(Zscores[i][0].size());
        //Enumerate_Posterior2(locus_log_marginals,Zscores[i], betas, Aijs[i], ld_matrices[i], max_causals, fullLikeli, prior_variance);
        // uSe my new function ???
        //Enumerate_Posterior(locus_log_marginals,Zscores[i], betas, Aijs[i], ld_matrices[i], max_causals, fullLikeli, prior_variance);
        exp_temp = locus_log_marginals.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }
    E_out.probs_locs = marginal_i;
    E_out.probs_stacked = stacker;
    return(fullLikeli);
}

double PAINTOR_Importance_Sampling(vector<vector<VectorXd>> &Zscores, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& ld_matrices, CausalProbs &E_out, double prior_variance, int num_samples, int sampling_seed){
    vector<VectorXd> marginal_i;
    VectorXd temp;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;
    for(unsigned int i = 0; i < Zscores.size(); i ++){
        VectorXd locus_log_marginals(Zscores[i][0].size());
        Locus_Importance_Sampler(locus_log_marginals, Zscores[i], betas, Aijs[i], ld_matrices[i], fullLikeli, prior_variance, num_samples, sampling_seed);
        exp_temp = locus_log_marginals.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }
    E_out.probs_locs = marginal_i;
    E_out.probs_stacked = stacker;

    return(fullLikeli);
}

vector<double> eigen2vec(VectorXd &vec){
    vector<double> outvec;
    for(unsigned i = 0; i < vec.size(); i++){
        outvec.push_back(vec[i]);
    }
    return(outvec);
}

vector<vector<double> > Stack_Matrices(vector<MatrixXd> &mats){
    vector<vector<double> > out_stack;

    for(unsigned i = 0; i < mats.size(); i++){
        MatrixXd tempMat = mats[i];
        for(unsigned j = 0; j < tempMat.rows(); j++){
            VectorXd temp = tempMat.row(j);
            out_stack.push_back(eigen2vec(temp));
        }
    }
    return(out_stack);
}

void Stack_EigenMatrices(vector<MatrixXd> &mats, MatrixXd& stacked_matrices){
    int place_holder = 0;
    for(unsigned int i = 0; i < mats.size(); i++){
        int num_snps = mats[i].rows();
        stacked_matrices.middleRows(place_holder, num_snps) = mats[i];
        place_holder += num_snps;
    }
}

void Copy_CausalProbs(CausalProbs& source, CausalProbs& destination ){
    unsigned int num_loci = source.probs_locs.size();
    if(destination.probs_locs.size()==0){
        for(unsigned int i = 0; i < num_loci ; i++){
            destination.probs_locs.push_back(source.probs_locs[i]);
        }
        destination.probs_stacked= source.probs_stacked;
    }
    else {
        for (unsigned int i = 0; i < num_loci; i++) {
            destination.probs_locs[i] = source.probs_locs[i];
        }
        destination.probs_stacked = source.probs_stacked;
    }
}

/*
   double EM_Run(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &beta_int, vector<MatrixXd> &annotations,vector<vector<MatrixXd>> &ld_matrix , double prior_variance, int num_samples, int sampling_seed){
//Run EM with sampling
vector<double> beta_run = eigen2vec(beta_int);
ObjectiveData Opt_in;
Opt_in.Aijs = Stack_Matrices(annotations);
int iterations = 0;
VectorXd beta_update;
double likeliOld;
double likeli =1;
double gradient_tolerance = 1e-5;
int max_ascents = 100;
int succesful_opt;
CausalProbs curr_probs;
double EM_likelihood;
while(iterations < iter_max) {
likeli = Estep(Zscores, beta_int, annotations, ld_matrix, probabilites, prior_variance, num_samples, sampling_seed);
VectorXd stacked_probabilities = vector2eigen(probabilites.probs_stacked);
MatrixXd stacked_annotations(stacked_probabilities.size(), beta_int.size());
Stack_EigenMatrices(annotations, stacked_annotations);
if(likeliOld <= likeli || iterations==0){
Copy_CausalProbs(probabilites, curr_probs);
VectorXd new_gammas(beta_int.size());
VectorXd return_values(2);
Gradient_Ascent(beta_int,  new_gammas, stacked_probabilities, stacked_annotations, gradient_tolerance, max_ascents, return_values);
likeliOld = likeli;
succesful_opt = return_values[0];
EM_likelihood = return_values[1];
if(succesful_opt == 0){
beta_int = new_gammas;
//   //cout << std::setprecision(9) << "Sum of log Bayes Factors at iteration: " << iterations+1 << ": " << likeli << endl;
// //cout << "Parameter Estimates:" << endl << new_gammas << endl << endl;
E=////cout << "EM Likelihood :" << EM_likelihood << endl;
}
else{
//cout << "Warning: Optimization unstable at iteration: " << iterations+1 << endl;
//cout << "Maximum liklelihood estimate for enrichments potentially not reached" << endl;
//cout << "Program is exiting and returning values corresponding to estimates from iteration: " << iterations << endl;
break;
}
}
else{
//cout << "Last Iteration did not improve model fit. Resampling at prevoius estimates. " << endl;
}
iterations++;
}
Copy_CausalProbs(curr_probs,probabilites);
return(EM_likelihood);
}
*/

double EM_Run(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &beta_int, vector<MatrixXd> &Aijs,vector<vector<MatrixXd>> &ld_matrix , double prior_variance,int max_causals){
    //Run EM with full enumeration
    vector<double> beta_run = eigen2vec(beta_int);
    ObjectiveData Opt_in;
    Opt_in.Aijs = Stack_Matrices(Aijs);
    int iterations = 0;
    VectorXd beta_update;
    double likeliOld;
    double likeli =1;
    double gradient_tolerance = 1e-5;
    int max_ascents = 1000;
    int succesful_opt;
    CausalProbs curr_probs;
    double EM_likelihood;
    //cout << "EM run" << endl;
    for(int i=0; i < beta_int.size(); i++){
        cout << beta_int[i] << endl;
    }
    beta_int(4);
    for(int i=0; i < beta_int.size(); i++){
        beta_int[i] = 0.25; 
        cout << beta_int[i] << endl;
    }
    while(iterations < iter_max) {
        likeli = Estep(Zscores, beta_int, Aijs, ld_matrix, probabilites, prior_variance, max_causals);
        VectorXd stacked_probabilities = vector2eigen(probabilites.probs_stacked);
        MatrixXd stacked_annotations(stacked_probabilities.size(), beta_int.size());
        Stack_EigenMatrices(Aijs, stacked_annotations);
        if(abs(likeliOld - likeli)> 0.01){
            Copy_CausalProbs(probabilites, curr_probs);
            VectorXd new_gammas(beta_int.size());
            VectorXd return_values(2);
            Gradient_Ascent(beta_int,  new_gammas, stacked_probabilities, stacked_annotations, gradient_tolerance, max_ascents, return_values);
            likeliOld = likeli;
            succesful_opt = return_values[0];
            EM_likelihood = return_values[1];
            if(succesful_opt == 0){
                beta_int = new_gammas;
                //cout << std::setprecision(9) << "Sum of log Bayes Factors at iteration: " << iterations+1 << ": " << likeli << endl;
                cout << "Parameter Estimates:" << endl << new_gammas << endl << endl;
                cout << "EM Likelihood :" << EM_likelihood << endl;
            }
            else{
                cout << "Warning: Optimization unstable at iteration: " << iterations+1 << endl;
                cout << "Maximum liklelihood estimate for enrichments potentially not reached" << endl;
                cout << "Program is exiting and returning values corresponding to estimates from iteration: " << iterations << endl;
                break;
            }
        }
        else {
            break;
        }
        iterations++;
    }
    Copy_CausalProbs(curr_probs,probabilites);
    return(EM_likelihood);
}

double PreCompute_Enrichment(int iter_max, vector<vector<VectorXd>> &Zscores,  VectorXd &gamma_intitial, vector<MatrixXd> &annotations,vector<vector<MatrixXd>> &ld_matrix , double prior_variance, int max_causals){
    //Run EM with full enumeration
    CausalProbs probabilites;
    vector<double> beta_run = eigen2vec(gamma_intitial);
    ObjectiveData Opt_in;
    Opt_in.Aijs = Stack_Matrices(annotations);
    int iterations = 0;
    VectorXd beta_update;
    double likeliOld;
    double likeli =1;
    double gradient_tolerance = 1e-5;
    int max_ascents = 1000;
    int succesful_opt;
    CausalProbs curr_probs;
    vector<vector <double>> test_data(3);
    vector <double> rowone;
    test_data.push_back(rowone);
    vector <double> rowtwo;
    test_data.push_back(rowtwo);
    vector <double> rowthree;
    test_data.push_back(rowthree);
    probabilites.prob_vectors = test_data;
    
    double EM_likelihood;
    for(int i=0; i < gamma_intitial.size(); i++){
        cout << "Initial gamma " << gamma_intitial[i] ;
    }
    VectorXd gamma_intitial2(3);
    for(int i=0; i < gamma_intitial2.size(); i++){
            gamma_intitial2[i] = 1/(1+exp(gamma_intitial[0]));
            cout << gamma_intitial2[i]<< endl;
    }

    while(iterations < iter_max) {
        //cout << "Iterations = " << iterations << endl;
        cout << gamma_intitial.size() << endl;
        cout << gamma_intitial2.size() << endl;
        cout << gamma_intitial2 << endl; 
        //for(int i=0; i < gamma_intitial2.size(); i++){
        //        gamma_intitial2[i] = gamma_intitial[0];
        //}
        likeli = Estep2(Zscores, gamma_intitial, annotations, ld_matrix, probabilites, prior_variance, max_causals, gamma_intitial2);
        cout << "Are we getting this far????" << endl;
        cout << "likli hood " << likeli <<endl;
        VectorXd stacked_probabilities = vector2eigen(probabilites.probs_stacked);
        //cout << "lel" << endl;
        MatrixXd stacked_annotations(stacked_probabilities.size()/4,1); 
        MatrixXd stacked_probabilities2(stacked_probabilities.size()/4, 4);
        int counter = 0;
        for(int i=0; i < stacked_probabilities2.cols(); i++){
            for (int j=0; j < stacked_probabilities2.rows(); j++){
                stacked_probabilities2(j,i) = probabilites.prob_vectors[i][j]; 
            }
        }
        ofstream myfile;
        myfile.open("probs.txt");
        for(int i=0; i < stacked_probabilities2.rows(); i++){
            //cout << "SNP NO " << i << " " << endl; 
            for (int j=0; j < stacked_probabilities2.cols(); j++){
                myfile << stacked_probabilities2(i,j) << " ";
                counter +=1;
            }
            myfile << endl;
        }
        myfile.close();
      //  exit(1);
      //  cout << "END OUTPUT" << endl;
        //cout << "Does this work ??" << endl;
        //cout << "Testing 1 2 " << endl;
        //cout << annotations.size()  << " " << stacked_annotations.cols() <<" "  << stacked_annotations.rows() << endl; 
        Stack_EigenMatrices(annotations, stacked_annotations);
        //cout << "Likelihood = " << likeli <<endl;
        if(abs(likeliOld - likeli)> 0.01){
            Copy_CausalProbs(probabilites, curr_probs);
            VectorXd new_gammas(gamma_intitial.size());
            VectorXd return_values(2);
            Gradient_Ascent(gamma_intitial,  new_gammas, stacked_probabilities, stacked_probabilities2, gradient_tolerance, max_ascents, return_values, gamma_intitial2,Zscores);
            //cout << "Do we complete gradient ascent "<< endl;
            likeliOld = likeli;
            succesful_opt = return_values[0];
            for(int i=0; i < gamma_intitial2.size(); i++){
                cout << gamma_intitial2[i] << endl;
            }
            //cout << "WHAT IN THE ACTUALY FUCK" << endl;
            EM_likelihood = return_values[1];
            if(0==0){
                cout << "Are we actually fucking converging WHAT ????? " << endl;
                cout << new_gammas << endl;
            }
            else{
                cout << "Warning: Optimization unstable at iteration: " << iterations+1 << endl;
                cout << "Maximum liklelihood estimate for enrichments potentially not reached" << endl;
                cout << "Program is exiting and returning values corresponding to estimates from iteration: " << iterations << endl;
                break;
            }
        }
        else {
            cout << gamma_intitial2<< endl;
            cout << "Ran twice " << endl; 
            cout << ":)  " << endl;
            exit(1);
            break;
        }
        iterations++;
    }
    cout << "Lilki ___ ummm " << likeli;
    Copy_CausalProbs(curr_probs,probabilites);
    return(likeli);
}

VectorXd Zscores2Post(VectorXd& Zs){
    VectorXd post(Zs.size());
    VectorXd Zsq  = Zs.array().square();
    for(int i = 0; i < Zsq.size(); i ++){
        VectorXd Ztemp = (Zsq.array() - Zsq[i])/2;
        VectorXd Zexp = Ztemp.array().exp();
        post[i] = 1/Zexp.sum();
    }
    return(post);
}

VectorXd vector2eigen(vector<double> &vec){
    VectorXd outVec(vec.size());
    for(unsigned i =0; i <vec.size(); i++){
        outVec[i] = vec[i];
    }
    return(outVec);
}

void Eigen2Vector(vector<double> &vector_in, VectorXd& eigen_out){
    for(unsigned int i =0; i < eigen_out.size(); i++){
        vector_in.push_back(eigen_out[i]);
    }
}

void PrintVector(vector<int>& in_vec){
    for(int i = 0; i < in_vec.size(); i++){
        cout << in_vec[i] << " ";
    }
    cout << endl;
}
