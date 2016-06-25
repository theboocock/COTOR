//
//

#include "Functions_Optimize.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen>
//#include "nlopt.hpp"
#include <iterator>
#include <string>
#include <iomanip>
#include <algorithm>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>


using namespace Eigen;
using namespace std;

double convert_exp(double input){
    return 1/(1+exp(input));
}

double Compute_Objective(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, vector <double> claudia_gammas){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double objective_value=0;
    double value_c1 = 0;
    double value_c0 = 0;
    double value_c2 = 0;
    double value_c3 = 0;
    double m0;
    cout << claudia_gammas[0] << " " << claudia_gammas[1] <<  " " << claudia_gammas[2];
    m0 = log(1 - convert_exp(claudia_gammas[0]) - convert_exp(claudia_gammas[1]) - convert_exp(claudia_gammas[2])); 
    cout << " Cout m0 " << m0 << endl;
    for(int i =0; i < total_snps/3; i++){
        value_c0 = (1- stacked_annotations.row(i)[1]) * (1- stacked_annotations.row(i)[2]) * log(1+exp(m0));
        value_c1 = (1- stacked_annotations.row(i)[1]) * stacked_annotations.row(i)[2] * log(1+exp(claudia_gammas[0]));
        value_c2 = (1- stacked_annotations.row(i)[2]) * stacked_annotations.row(i)[1] * log(1+exp(claudia_gammas[1]));
        value_c3 = (stacked_annotations.row(i)[1]) * stacked_annotations.row(i)[2] * log(1+exp(claudia_gammas[2]));
        cout << "Values " << value_c0 << " " << value_c1 << " " << value_c2 << " " << value_c3 << endl;
        //objective_value = objective_value -value_c0 - value_c1 - value_c2 - value_c3; 
    }
	cout << "Objective " << endl;
	cout << objective_value << endl;    
    return objective_value;
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

double Compute_Objective_perSNP(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, vector <double> claudia_gammas, vector<vector<VectorXd>> &Zscores){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double objective_value=0;
    double value_c1 = 0;
    double value_c0 = 0;
    double param_sum = 0;
    double value_c2 = 0;
    double m0,m1,m2,m3;
    m2 = -1000;
    m3 = -1000;
    m0 = -1000;
    m1 = -1000;
    double objective;
    double double_tmp, double_tmp2, double_tmp3;
    cout << "Are we running this objective function () " << endl; 
    for(int overall_idx =0; overall_idx < total_snps/4; overall_idx++){
            //cout << "SNP #  " << i <<  endl;
            //cout << "Total SNPs" << total_snps << endl;
            //cout << "Stacked annotations " << stacked_annotations.rows() << endl;
            m0  = (1 - (1/(1+exp(claudia_gammas[0]))) - (1/(1+exp(claudia_gammas[1])))- (1/(1+exp(claudia_gammas[2]))));
            m0 = log(1/m0 - 1);
            m0 = (1-stacked_annotations.row(overall_idx)[1] - stacked_annotations.row(overall_idx)[2] - stacked_annotations.row(overall_idx)[0]) * log(1/(1+exp(m0)));
            double_tmp = (stacked_annotations.row(overall_idx)[1]) * log(1+exp(claudia_gammas[1]));
            double_tmp2 = (stacked_annotations.row(overall_idx)[2]) * log(1+exp(claudia_gammas[2]));
            double_tmp3 = (stacked_annotations.row(overall_idx)[0]) * log(1+exp(claudia_gammas[0]));
            objective = objective + m0 - double_tmp - double_tmp2 - double_tmp3;
    }
    cout <<  "Objective so far " << objective << endl;
    return(-objective);
}
double Compute_Objective3(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, vector <double> claudia_gammas, vector<vector<VectorXd>> zscores){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double objective_value=-1000;
    double value_c1 = 0;
    double value_c0 = 0;
    double param_sum = 0;
    double value_c2 = 0;
    double m0,m1,m2,m3;
    m2 = -1000;
    m2 = -1000;
    m3 =-1000;
    for(int j =0; j < zscores.size(); j++){
        int zscore_length = zscores[j][0].size();
        for(int i =0; i < zscore_length; i++){
                //cout << "SNP #  " << i <<  endl;
                //cout << "Total SNPs" << total_snps << endl;
                //cout << "Stacked annotations " << stacked_annotations.rows() << endl;
                param_sum = 0;
                value_c1 = 0;
                //cout << (stacked_annotations.row(i)[0]) * (claudia_gammas[0]) << endl;
                m1= LogSum(m1,(log(stacked_annotations.row(i)[0]) + log(claudia_gammas[0])));
                m2= LogSum(m2,(log(stacked_annotations.row(i)[1]) + log(claudia_gammas[1])));
                m3= LogSum(m3,(log(stacked_annotations.row(i)[0]) + log(stacked_annotations.row(i)[1]) + log(claudia_gammas[2])));
            }
    m0 = log(1- claudia_gammas[0] - claudia_gammas[1] - claudia_gammas[2]);
    double tmp = 1 + exp(m0-m3) + exp(m1 -m3) + exp(m2-m3); 
    cout << "M1 " << m1 << " " << m2 << " " << m3 << endl; 
    tmp = m3 + log(tmp);
    objective_value += tmp;
    }
    return(-objective_value);
}
double Compute_Objective2(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, vector <double> claudia_gammas, vector<vector<VectorXd>> &Zscores){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double objective_value=-1000;
    double value_c1 = 0;
    double value_c0 = 0;
    double param_sum = 0;
    double value_c2 = 0;
    double m0,m1,m2,m3;
    m2 = -1000;
    m3 = -1000;
    m0 = -1000;
    m1 = -1000;
    double overall_idx = 0;
    double objective = -1000;
    double double_tmp, double_tmp2, double_tmp3;
    cout << "Are we running this objective function () " << endl; 
    for(int j =0; j < Zscores.size(); j++){
        int zscore_length = Zscores[j][0].size();
    m2 = -1000;
    m3 = -1000;
    m0 = -1000;
    m1 = -1000;
    for(int i =0; i < zscore_length; i++){
            //cout << "SNP #  " << i <<  endl;
            //cout << "Total SNPs" << total_snps << endl;
            double_tmp = log(stacked_annotations.row(overall_idx)[1]) - log(1+exp(5.304371));
            double_tmp2 = log(stacked_annotations.row(overall_idx)[2]) - log(1+exp(5.304371));
            double_tmp3 = log(stacked_annotations.row(overall_idx)[0]) - log(1+exp(5.304371));
            m1 = LogSum(m1, double_tmp);
            m2 = LogSum(m2, double_tmp2);
            m3 = LogSum(m3, double_tmp3);
            overall_idx += 1;
    }
    m0  = (2 - (2/(1+exp(claudia_gammas[0]))) - (2/(1+exp(claudia_gammas[1])))- (2/(1+exp(claudia_gammas[2]))));
    cout << m0 << endl; 
    m0 = log(m0);
    m1 = m1 + log(2/(1+exp(claudia_gammas[1])));
    m2 = m2 + log(2/(1+exp(claudia_gammas[2])));
    m3 = m3 + log(2/(1+exp(claudia_gammas[0])));
    cout << "M1 " << m1 << " M2 " << m2 << "M3 " << m3 << endl;
    cout << "idx " << overall_idx << endl;
    cout << "Gammas " << stacked_annotations.row(overall_idx- 1 )[0] << " M2 " << stacked_annotations.row(overall_idx - 1)[1]  << "M3 " << stacked_annotations.row(overall_idx- 1 )[2] << endl;
    double tmp = exp(m0) + exp(m1) + exp(m2)+ exp(m3);
    cout << " TMP 1 " << log(tmp) << endl;
    tmp = 1+exp(m0-m3)+ exp(m1-m3)+exp(m2-m3);
    tmp = m3+log(tmp);
    cout << " TMP 2 " << tmp << endl;
    objective = objective  + tmp; 
    cout <<  "Objective so far " << objective << endl;
    }
    return(-objective);
}

struct f_params{
    VectorXd current_gammas;
    VectorXd stacked_probabilites;
    MatrixXd stacked_annotations;
    vector <double> claudia_gammas;
    vector<vector<VectorXd>> Zscores;
};

double GSL_llk(const gsl_vector *x, void *params){
	//first set times
	int na = 3; 
    
	for (int i = 0; i < na; i++){
		((struct f_params *) params)->claudia_gammas[i]= gsl_vector_get(x, i);
	}
    
    VectorXd current_gammas = ((struct f_params *) params)->current_gammas;
    VectorXd stacked_probabilites = ((struct f_params *) params)->stacked_probabilites;
    MatrixXd stacked_annotations = ((struct f_params *) params)->stacked_annotations;
    vector<vector<VectorXd>> Zscores = ((struct f_params *) params)->Zscores;
    double object = Compute_Objective_perSNP(current_gammas, stacked_probabilites, stacked_annotations, ((struct f_params *) params)->claudia_gammas, Zscores);
    return object;
}


void Maximise(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, VectorXd& claudia_gammas,vector<vector<VectorXd>> &Zscores){
        double size;
        int nparam = 3;
        int status; 
        cout <<  "YAY" << endl;
        const gsl_multimin_fminimizer_type *T =
         gsl_multimin_fminimizer_nmsimplex2;
         gsl_multimin_fminimizer *s;
         gsl_vector *x;
         gsl_vector *ss;
         gsl_multimin_function lm;
         lm.n = 3;
         cout << claudia_gammas << endl;
         lm.f = &GSL_llk;
		 vector <double> claudia_gammas2(claudia_gammas.size());
		 for(int i = 0; i < 3; i++){
			claudia_gammas2[i] = log(1/claudia_gammas[i] -1);
            cout << claudia_gammas2[i] << endl; 
         }
         struct f_params f = {current_gammas, stacked_probabilites, stacked_annotations, claudia_gammas2,Zscores};
         lm.params = &f;
         x = gsl_vector_alloc (lm.n);
         for (int i = 0; i < 3; i++) gsl_vector_set(x, i,(claudia_gammas2[i])); 
         
         ss = gsl_vector_alloc(3);
         gsl_vector_set_all(ss, 0.1);
         s = gsl_multimin_fminimizer_alloc (T, nparam);
         cout << " Do we iterate ??? " << endl;
		 gsl_multimin_fminimizer_set (s, &lm, x, ss);
         int iter = 0;
			do
			 {
					 iter++;
					 status = gsl_multimin_fminimizer_iterate (s);

					 if (status){
							 printf ("error: %s\n", gsl_strerror (status));
							 break;
					 }
					 size = gsl_multimin_fminimizer_size(s);
					 status = gsl_multimin_test_size (size, 0.01);
					 //cout << iter << " "<< iter %10 << "\n";
                     claudia_gammas2 =  (f).claudia_gammas;
			         cout <<"iteration: "<< iter << " "<< claudia_gammas2[0]<< " "<< claudia_gammas2[1]<< " "<< claudia_gammas2[2] << endl; 
						 cout << " "<< s->fval << " "<< size <<  "\n";
			 }
			 while (status == GSL_CONTINUE && iter <5000);
			 if (iter > 4999) {
					 cerr << "WARNING: failed to converge. Maybe try -mcmc to check parameter estimates.\n";
					 //exit(1);
		} 
    for(int i = 0; i < 3; i ++){
        claudia_gammas[i]= 1/(1+exp(claudia_gammas2[i])); 
        cout << "Final gammas " << claudia_gammas[i] << endl;
    }
}

void Compute_Gradient(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, VectorXd&gradient, VectorXd& claudia_gammas){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double scaling_factor_c1;
    double scaling_factor_c0;
    double objective_value=0;
    VectorXd scaling_vectorc_1(gradient.size());
    VectorXd scaling_vectorc_2(gradient.size());
    gradient.setZero();
    double value_c1 = 0;
    double value_c0 = 0;
    double param_sum = 0;
    double value_c2 = 0;
    double value_c3 = 0;
    double m0;
    m0 = log(1/(1- claudia_gammas[0] - claudia_gammas[1] - claudia_gammas[2]) -1);
    for(int i =0; i < total_snps/3; i++){
            value_c0 = (1- stacked_annotations.row(i)[1]) * (1- stacked_annotations.row(i)[2]) * log(1+exp(m0));
            value_c1 = (1- stacked_annotations.row(i)[1]) * stacked_annotations.row(i)[2] * log(1+exp(claudia_gammas[0]));
            value_c2 = (1- stacked_annotations.row(i)[2]) * stacked_annotations.row(i)[1] * log(1+exp(claudia_gammas[1]));
            value_c3 = (stacked_annotations.row(i)[1]) * stacked_annotations.row(i)[2] * log(1+exp(claudia_gammas[2]));
        
        }
}

void nelder_mead(){
 	    // call minimization routine
//    Maths::Optimization::nelder(Compute_Objective, 3, gamma_claudia, 1E-8, 500);
 
    // display results
    std::cout << "Minimization points: " << std::endl;
  //  for (int i = 0 ; i <= N; i++) {
  //      for (int j = 0; j < N; j++)
  //          std::cout << "  " << std::setprecision(10) << P[i][j];
  //      std::cout << std::endl;
  //  }
 

}


void Compute_Hessian(VectorXd& current_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, MatrixXd& hessian){
    unsigned int total_snps = stacked_probabilites.size();
    double exponential_term;
    double scaling_factor_c1;
    double scaling_factor_c0;
    hessian.setZero();
    MatrixXd outer_product;
    for(int i = 0; i < total_snps; i++){
        outer_product = stacked_annotations.row(i).transpose()*stacked_annotations.row(i);
        exponential_term = current_gammas.dot(stacked_annotations.row(i));
        scaling_factor_c1 =  stacked_probabilites[i]*(exp(-1*exponential_term))/((1+exp(-1*exponential_term))*(1+exp(-1*exponential_term)));
        scaling_factor_c0 = (1-stacked_probabilites[i])*(exp(exponential_term))/((1+exp(exponential_term))*(1+exp(exponential_term)));
        hessian = hessian - scaling_factor_c1*outer_product - scaling_factor_c0*outer_product;
    }
}

//int Gradient_Ascent(VectorXd& current_gammas, VectorXd& new_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, double gradient_tolerance, int max_iterations){
//    unsigned int num_parameters = current_gammas.size();
//    VectorXd gradient(num_parameters);
//    VectorXd gamma_iterate(num_parameters);
//    int max_line_search = 500;
//    int search_counter;
//    double tuner;
//    double new_objective;
//    double current_objective = Compute_Objective(current_gammas, stacked_probabilites, stacked_annotations);
//    for(int i = 0; i < max_iterations;i++){
//        tuner = 0.9;
//        Compute_Gradient(current_gammas, stacked_probabilites,stacked_annotations, gradient);
//        new_gammas = current_gammas + tuner * gradient;
//        new_objective = Compute_Objective(new_gammas, stacked_probabilites, stacked_annotations);
//        if(gradient.norm()>gradient_tolerance){
//            search_counter=0;
//            while (new_objective - current_objective <= 0){
//                search_counter++;
//                if(search_counter > max_line_search){
//                    return -9;
//                }
//                else {
//                    tuner = tuner * .9;
//                    new_gammas = current_gammas + tuner * gradient;
//                    new_objective = Compute_Objective(new_gammas, stacked_probabilites, stacked_annotations);
//                }
//            }
//            current_objective = new_objective;
//            current_gammas = new_gammas;
//        }
//        else{
//            break;
//        }
//    }
//    return 0;
//}

void Gradient_Ascent(VectorXd& current_gammas, VectorXd& new_gammas, VectorXd& stacked_probabilites, MatrixXd& stacked_annotations, double gradient_tolerance, int max_iterations, VectorXd& return_values, VectorXd &gamma_initial2,vector<vector<VectorXd>> &Zscores){
    // TODO: remove this
    unsigned int num_parameters = gamma_initial2.size(); 
    VectorXd gradient(num_parameters);
    VectorXd gamma_iterate(num_parameters);
    int max_line_search = 500;
    int search_counter;
    double tuner;
    double new_objective;
    VectorXd current_gammas2(gamma_initial2.size());
    Maximise(current_gammas, stacked_probabilites, stacked_annotations, gamma_initial2, Zscores);
    return_values[0] = 0;
}
