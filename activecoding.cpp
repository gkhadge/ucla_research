#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <ctime>
// #include <thread>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <random>
// #include <boost/random/uniform_int.hpp>
#include <fstream>
#include <assert.h>
// #include <chrono>
using namespace std;

static unsigned int seed;

long double binary_entropy(const long double p)
{
	// Hbp = binary_entropy(p)
	//  Returns the binary entropy of probability p
	//   must have 0<=p<=1

	// if (p < 0 || p > 1)
	// {
	//     printf("ERROR: p must be a probability (between 0 and 1)!");
	//     return -1;
	// }

	// calculate the binary entropy
	long double Hbp = (p==0 || p==1)  ? 0 : -p*log2(p) - (1-p)*log2(1-p);
	return Hbp;
}

// function prob_vec_next = update_probs(prob_vec_current, y, x_vec, del)
void update_probs(const long double* prob_vec_current, const char y, const char* x_vec, const long double del, const unsigned long long M, long double* prob_vec_next) 
{
	// function prob_vec_next = update_probs(prob_vec_current, y, x_vec, del)
	//
	//  Update belief vector (posterior probabilities) at time (t+1)
	//   according to eq. in App. A2, based on belief vector at time t and the
	//   transition probabilities
	//
	//   Arguments:
	//   prob_vec_current = 1xM vector of beliefs at time t
	//   y = channel output at time (t+1), a scalar
	//   x_vec = 1xM vector of channel outputs at time (t+1)
	//   del = crossover probability of BSC
	//
	//   Returns:
	//   prob_vec_next = 1xM vector of beliefs at time (t+1)
	//
	//   Note:
	//   This function was initially called update_probs_test.m, and was
	//   determined to be faster than the previous update_probs.m, so the name
	//   was updated.

	//M = length(prob_vec_current);

	// Determine transition probability P(Y|X) for time (t+1), on the BSC
	// dist_vec = (x_vec ~= y);

	// Update belief for each j in {1,...,M}
	//prob_vec_next = zeros(size(prob_vec_current));
	//int prob_vec_next[M];
	// for (int i=0;i<M;++i){
	// 	prob_vec_next[i] = 0;
	// }

	long double prob_vec_next_sum = 0;
	for (int mdx = 0; mdx < M; ++mdx) {
		//printf("\n %d  %d", x_vec[mdx], y);
	    // if (x_vec[mdx] != y)
	    //     prob_vec_next[mdx] = prob_vec_current[mdx] * del;
	    // else
	    //     prob_vec_next[mdx] = prob_vec_current[mdx] * (1-del);

	    prob_vec_next[mdx] = prob_vec_current[mdx] * ((x_vec[mdx] != y) ? del : (1-del));

	    prob_vec_next_sum += prob_vec_next[mdx];
	}

	// normalize
	//prob_vec_next = prob_vec_next ./ sum(prob_vec_next);
	//printf("\n");
	for (unsigned long long i=0; i<M; ++i){
		prob_vec_next[i] = prob_vec_next[i]/prob_vec_next_sum;
		//printf("%4.4f ", prob_vec_next[i]);
	}

}


// //function [set_S0, set_S1] = partition_beliefs(prob_vec)
// void partition_beliefs(const double* prob_vec, const int M, int* set_num)
// {
// 	// Partition beliefs into two sets S0 and S1 that are as close together in
// 	// probability, but S0 is greater than S1
// 	//
// 	//   prob_vec = M-vector of probabilities that add to 1
// 	//
// 	// Implements Algorithm 2 of
// 	//   [4] Naghshvar & Javidi, "Extrinic Jensen-Shannon (EJS) Divergence", arXiv
// 	//   2013
// 	// Complexity is of order O(M^2)

// 	int s0_ind = 0;

// 	// Initialize the probability weights of the sets
// 	double prob_S0 = 1;
// 	double prob_S1 = 0;
// 	//set_S0 = [1:M];
// 	//set_S1 = [];

// 	//int set_num[M]; //s0_ind if in S0, !s0_ind if in S1
// 	for (int i=0; i<M;++i){
// 		set_num[i] = s0_ind;
// 	}

// 	double min_prob = 0;       // smallest probability of element in S0
// 	double delta = 1;          // difference between set probabilities of S0 and S1

// 	// Initialize prob_k and kdx
//     double prob_k = 1; // min element
//     int kdx;
//     double prob_max = -1;
//     int k_max = 0;
//    	for (int i=0; i<M; ++i){
// 		if (prob_vec[i] < prob_k){
// 			prob_k = prob_vec[i];
// 			kdx = i;
// 		}
// 		if (prob_vec[i] > prob_max){
// 			prob_max = prob_vec[i];
// 			k_max = i;
// 		}
//     }
//     if (prob_max >= 0.5) {
//     	for (int i=0; i<M; ++i){
//     		set_num[i] = (i == k_max);
//     	}
//     } else {
// 	    int n = 0;
// 		while (min_prob < delta) {
// 		    // Find the element in S0 w/ smallest probability
// 		    //[prob_k, kdx] = min(prob_vec(set_S0));

// 		    // if isempty(kdx) {
// 		    //     //error('Index of smallest element is empty.');
// 		    //     printf("Index of smallest element is empty.");
// 		    // }

// 		    // Put this element in S1, take it out of S0
// 		    //set_S1 = [set_S1 set_S0(kdx)];
// 		    //set_S0 = [set_S0(1:kdx-1) set_S0(kdx+1:length(set_S0))];
// 		    set_num[kdx] = !s0_ind;
		    
// 		    // if (length(set_S0) + length(set_S1) != M) {
// 		    //     error('Error in switching set elements. Investigate further.');
// 		    // }
// 		    // Update the set probabilities
// 		    prob_S0 = prob_S0 - prob_k;
// 		    prob_S1 = prob_S1 + prob_k;
		    
// 		    // prob_S1 = 0;
// 		    // for (int i=0; i<M; ++i){
// 		    // 	printf("\n\n%d  %4.4f",set_num[i],prob_vec[i]);
// 		    // 	if (set_num[i])
// 		    // 		prob_S1 += prob_vec[i];
// 		    // }
// 		    // prob_S0 = 1 - prob_S1;
// 		    // printf("\n%d",n);
// 		    // n+=1;

// 		    // Check which set is biggest now
// 		    if (prob_S0 < prob_S1) {
// 		        // Swap sets
// 		        s0_ind = !s0_ind;
// 		        // set_temp = set_S0;
// 		        // set_S0 = set_S1;
// 		        // set_S1 = set_temp;
// 		    	// for (int i=0; i<M; ++i){
// 		    	// 	set_num[i] = !(set_num[i]);
// 		    	// }
// 		        // Update the set probabilities
// 		        double prob_temp = prob_S0;
// 		        prob_S0 = prob_S1;
// 		        prob_S1 = prob_temp;
// 		    }
		    
// 		    // Compute the difference between sets (which is positive now)
// 		    delta = prob_S0 - prob_S1;
// 		    // Compute the smallest probability of an element in S0
// 		    //min_prob = min(prob_vec(set_S0));

// 		    prob_k = 2;
// 		   	for (int i=0; i<M; ++i){
// 		    	if(set_num[i] == s0_ind){
// 		    		if (prob_vec[i] < prob_k){
// 		    			prob_k = prob_vec[i];
// 		    			kdx = i;
// 		    		}
// 		    	}
// 		    }
// 		    min_prob = prob_k;
// 		}	
// 	}
// 	// Correct the set output if s0 is inverted to 1 instead of 0.
// 	if (s0_ind == 1) {
//     	for (int i=0; i<M; ++i){
//     		set_num[i] = !(set_num[i]);
//     	}
// 	}
// }



//function [set_S0, set_S1] = partition_beliefs(prob_vec)
void partition_beliefs(const long double* prob_vec, const unsigned long long M, char* set_num)
{
	// Partition beliefs into two sets S0 and S1 that are as close together in
	// probability, but S0 is greater than S1
	//
	//   prob_vec = M-vector of probabilities that add to 1
	//
	// Implements Algorithm 2 of
	//   [4] Naghshvar & Javidi, "Extrinic Jensen-Shannon (EJS) Divergence", arXiv
	//   2013
	// Complexity is of order O(M^2)

    auto compare = [&] (unsigned long long i1, unsigned long long i2) {
	    return prob_vec[i1] > prob_vec[i2];
	};

    auto compare_less = [&] (unsigned long long i1, unsigned long long i2) {
	    return prob_vec[i1] < prob_vec[i2];
	};

	auto insert_sorted = [&] ( std::vector<unsigned long long> & vec, unsigned long long const& item )
	{
	    vec.insert(std::upper_bound(vec.begin(), vec.end(), item, compare), item );
	    return;
	};

	std::vector<unsigned long long> indices(M) ; // vector with 100 ints.
	std::iota (indices.begin(), indices.end(), 0); // Fill with 0, 1, ..., 99.

    unsigned long long k_max = *std::max_element(indices.begin(), indices.end(), compare_less);
	long double prob_max = prob_vec[k_max];
    if (prob_max >= 0.5) {
    	for (unsigned long long i=0; i<M; ++i){
    		set_num[i] = (i == k_max);
    	}
    	return;
    }

    std::sort(indices.begin(), indices.end(), compare);
    
    // for (vector<int>::const_iterator i = indices.begin(); i != indices.end(); ++i)
    //     cout << *i << ": " << prob_vec[*i] << endl;

    // cout << endl;
	long double prob_A = 0;
	long double prob_B = 1;
	// double prob = 0;
	std::vector<unsigned long long>::iterator it = indices.begin();
	while(true) {
	// for(std::vector<int>::iterator it = indices.begin(); prob_vec[*it] < (prob_S0 - prob_S1) ; ++it) {
	    /* std::cout << *it; ... */

	    if (prob_vec[*it] > (prob_B - prob_A)) 
	    	break;
	    // prob += prob_vec[*it];
	    prob_A += prob_vec[*it];
	    prob_B -= prob_vec[*it];
	    // cout << *it << ": " << prob_vec[*it] << endl;
	    
	    ++it;
	}
	// cout << "prob_S0: " << prob_A << endl;
	// cout << "prob_S1: " << prob_B << endl;
	
    std::vector<unsigned long long> setA(indices.begin(), it);
    std::vector<unsigned long long> setB(it, indices.end());
    
    while (true) {
    	// if (prob_A + prob_B != 1) {
	    // 	printf("\nProb A: %4.4f \t Prob B: %4.4f", prob_A, prob_B);
	    // 	printf("\nProb A + Prob B - 1: %e ",  prob_A + prob_B - 1);
	    // 	printf("\nMin Prob A: %e ",  prob_vec[setA.back()]);
	    // 	printf("\nMin Prob B: %e ",  prob_vec[setB.back()]);
    	// }

    	if (prob_A + prob_B - 1 > 1e-200){
	    	// printf("\nProb A: %4.4f \t Prob B: %4.4f", prob_A, prob_B);
	    	// printf("\nProb A + Prob B - 1: %e ",  prob_A + prob_B - 1);
	    	// printf("\nMin Prob A: %e ",  prob_vec[setA.back()]);
	    	// printf("\nMin Prob B: %e ",  prob_vec[setB.back()]);
	    	prob_A = 0;
		    for (vector<unsigned long long>::const_iterator i = setA.begin(); i != setA.end(); ++i)
		    	prob_A += prob_vec[*i];	    

	    	prob_B = 1 - prob_A;
		    // for (vector<int>::const_iterator i = setB.begin(); i != setB.end(); ++i)
		    // 	prob_B += prob_vec[*i];	
    		// printf("\nProb A + Prob B - 1: %e ",  prob_A + prob_B - 1);
			assert (prob_A + prob_B - 1 < 1e-200); 
    	}

        if (prob_A > prob_B){
            long double prob = prob_vec[setA.back()];
            if (prob < (prob_A - prob_B)){
                prob_B += prob;
                prob_A -= prob;
                insert_sorted(setB,setA.back());
                setA.pop_back();
            } else {
                break;
            }
        } else {
            long double prob = prob_vec[setB.back()];
            if (prob < (prob_B - prob_A)){
                prob_A += prob;
                prob_B -= prob;
                insert_sorted(setA,setB.back());
                setB.pop_back();
            } else {
                break;
            }
        }
    }
 //    if (prob_A + prob_B != 1){
 //    	printf("\nProb A: %4.4f \t Prob B: %4.4f", prob_A, prob_B);
 //    	printf("\nProb A + Prob B - 1: %e ",  prob_A + prob_B - 1);
 //    	printf("\nMin Prob A: %e ",  prob_vec[setA.back()]);
 //    	printf("\nMin Prob B: %e ",  prob_vec[setB.back()]);

 //    	prob_A = 0;
 //    	prob_B = 0;

	//     for (vector<int>::const_iterator i = setA.begin(); i != setA.end(); ++i)
	//     	prob_A += prob_vec[*i];	    
	//     for (vector<int>::const_iterator i = setB.begin(); i != setB.end(); ++i)
	//     	prob_B += prob_vec[*i];	    
 //    	printf("\nProb A + Prob B - 1: %e ",  prob_A + prob_B - 1);
		
	// 	assert (prob_A + prob_B  == 1);
	// 	assert (1 == 0);
	// }
    // if (probA + probB != 1)
    // 	throw std::runtime_error("Error probability does not add up");

	char setA_val = prob_B > prob_A;
	char setB_val = !setA_val;
    for (vector<unsigned long long>::const_iterator i = setA.begin(); i != setA.end(); ++i)
    	set_num[*i] = setA_val;	    
    for (vector<unsigned long long>::const_iterator i = setB.begin(); i != setB.end(); ++i)
    	set_num[*i] = setB_val;
}

//function [Etau_sim, num_msg_errors, num_msgs] = active_coding_M(M, del, epsi, num_msgs, coding_scheme)
void active_coding_M(unsigned long long M, long double del, long double epsi, unsigned long num_msgs, int coding_scheme, unsigned long *sum_tau_out, unsigned long *_num_msg_errors, int print_output)
{

	// function [Etau_sim, num_msg_errors, num_msgs] = active_coding_M(M, del, epsi, num_msgs, coding_scheme)
	//
	//   Sequential coding scheme for this particular number of messages M, BSC 
	//   crossover probability del, error probability epsi, and number of
	//   messages to simulate num_msgs.
	//
	//   coding_scheme is the particular coding scheme used
	//
	//   [2] Naghshvar & Javidi, "Optimal reliability over a class of
	//   binary-input channels with feedback", ITW 2012
	//   
	//   Includes an implementation of Algorithm 2 of [4] in partition_beliefs.m
	//   [4] Naghshvar & Javidi, "Extrinic Jensen-Shannon (EJS) Divergence", arXiv
	//   2013

	//Random number generator
	//default_random_engine generator (seed++);
	mt19937 generator (seed);
	// std::uniform_int_distribution<unsigned long long> dist0M(0,M-1);
	std::uniform_int_distribution<> dist0M(0,M-1);
	bernoulli_distribution dist_bern(del);
	srand(seed);
	//// Simulation parameters

	// BSC channel capacity
	long double capacity = 1 - binary_entropy(del); 

	// Debugging switch
	bool DBG_ON = false; 
	// DBG_ON = 1; 

	// Max. number of noise samples per message
	int num_bits_per_msg = ceil(log2(M)/capacity * 20);
	// printf("\ncapacity: %4.4f", capacity);
	// printf("\nnum_msgs: %u", num_msgs);
	// printf("\nnum_bits_per_msg: %u", num_bits_per_msg);

	// Initialize the number of message errors
	unsigned long num_msg_errors = 0;

	// Scheme from [2] & [4] for binary channels, based on partitioning messages
	// into two sets based on posterior probabilities (beliefs)
	int belief_sets_scheme = 1;
	// Generalized Horstein-Burnashev-Zigangirov scheme (deterministic)
	int GHBZ_determ_scheme = 2;
	// Generalized Horstein-Burnashev-Zigangirov scheme (randomized)
	int GHBZ_random_scheme = 3;

	//// Plot color constants

	// HSV_COLORS = [
	//     .8      0       1;  // purple
	//     1       0       0;  // red
	//     1       .6      0;  // orange
	//     .7      .2      0;  // brown
	//     .2      .6      0;  // dark green
	//     0.0833    1.0000    0.5882;
	//     .3      1       .59];   
	// PURPLE_IDX = 1;

	//// Iterate through each of the messages to be sent
	// (so we don't have to allocate such a huge amount of memory for all the
	// arrays)
	       
	if (print_output >= 2) {
		printf("\n");
		printf("\ncapacity: %4.4f",capacity);
		printf("\nnum_bits_per_msg: %d",num_bits_per_msg);
	}
	//// Run the simulation

	// fprintf('\n');

	// Previously mdx was the message index from 1 to num_msgs, but now it is a
	// constant
	unsigned long mdx = 1;

	unsigned long sum_tau = 0;

	if (coding_scheme == belief_sets_scheme) {	    
	    ////
	    // Deterministic message-partitioning scheme based on beliefs
	    //
	    // double rand01;
	    // double prob_vec_next[M], prob_vec_current[M];
	    long double array_one[M];
	    long double array_two[M];
		long double *prob_vec_next = array_one;
		long double *prob_vec_current = array_two;
	    long double *temp;
	    unsigned long long W, W_hat;
	        // Decoded message \hat W
	    unsigned long tau;
	    char x_vec[M];
	    char tx_bit, rx_bit;

	    for (unsigned long msg_idx = 1; msg_idx <= num_msgs; ++msg_idx) {
	        
	        // Modulo-2 additive noise for BSC(del)
	        //noise_vec = (rand(1, num_bits_per_msg) < del);

	  //   	int noise_vec[num_bits_per_msg]; 
	  //   	//printf("\n");
	  //   	for (int i = 0; i < num_bits_per_msg; ++i) {
			// 	double rand01 = rand() / double(RAND_MAX);
			// 	noise_vec[i] = (int) (rand01 < del);
			// 	//printf("%d ",noise_vec[i]);
			// }

	        // Randomly assign the true message W
	        // unsigned long W;
	        if (DBG_ON) {
	            // Dummy message (always 1)
	            W = 1;
	        } else {
	            // Random assignment
	            //W = ceil(rand(1)*M);
	            // rand01 = rand() / (double(RAND_MAX) + 0.0001);
	            // W = floor(rand01*M);

	            W = dist0M(generator);

	            //printf("\nrand01 = %4.4f", rand01);
	            //printf("\nmsg_idx = %u, W = %u\n", msg_idx, W);
	        }
	        // Stopping times
	        tau = 0;
	        
	        
	        // Start with a new message: initialize posterior probabilities at the
	        // receiver
	        //prob_vec_current = ones(1,M)./M;
	        // double prob_vec_current[M];
	        for(unsigned long long i=0;i<M;++i){
	        	prob_vec_current[i] = 1.0/M;
	        }
	        
	        //printf("... simulating message %u\n", msg_idx);
	        // if (print_output >= 2)
		       //  if ( msg_idx % 10000 == 0)
		       //      printf("\n... simulating message %lu", msg_idx);
	        
	        // Start sending symbols across the channel (tdx is the transmission
	        // index)
		    // unsigned long tdx = 0;
	        //for(int tdx = 0; tdx < num_bits_per_msg; ++tdx){
	        while(true) {
	            ++tau;
	            // Partition messages into two sets based on the current beliefs
	            // int x_vec[M];
	            //[set_S0, set_S1] = partition_beliefs(prob_vec_current, M, set_num);
				partition_beliefs(prob_vec_current, M, x_vec);
				// printf("\n");
				// for(int i=0;i<M;++i)
				// {
				// 	printf("%d ", x_vec[i]);
				// }
	            
	            // Create a vector of the bits that would be transmitted in this
	            // time slot for each of the possible messages
	            //x_vec = zeros(1,M);
	            //int x_vec[M];

	            // for (int jdx = 0; jdx < M; ++jdx) {
	            //     // if sum(set_S1 == jdx) {
	            //     //     // Would send a 1 if message j is in set S_1
	            //     //     x_vec[jdx] = 1;
	            //     // } else {    // if sum(set_S0 == jdx)
	            //     //     // Would send a 0 if message j is in set S_0
	            //     //     x_vec[jdx] = 0;
	            //     // }
	            //     x_vec[jdx] = set_num[jdx];
	            // }
	            
	            // Decide what to transmit in this time slot
	            // Select the bit corresponding to the true message
	            // (Sends a 0 if the true message is in set S_0,
	            //  Send a 1 if the true message is in set S_1)
	            tx_bit = x_vec[W];
	            
	            // Corrupt the bit by noise (representing reception across the noisy
	            // channel)
	            // int rx_bit = (tx_bit + noise_vec[tdx]) % 2;

				// rand01 = rand() / double(RAND_MAX);				
	            // rx_bit = (rand01 < del) ? !tx_bit : tx_bit;

	            rx_bit = dist_bern(generator) ? !tx_bit : tx_bit;

	            // printf("\n");
	            // printf("%d",rx_bit);
	            
	            // update belief vector (posterior probabilities) according to eq. in App.
	            // A2 of [2]
	            // The receiver knows the encoding rule, so it can compute what the
	            // Tx outputs would be for each possible message
	            
	            update_probs(prob_vec_current, rx_bit, x_vec, del, M, prob_vec_next);
	            
	            // [Added 4/15/14 for PhD defense]
	            // if (msg_idx == 1) {
	            //     // Save the current belief vector for later                
	            //     prob_vec_evolution(tdx,:) = prob_vec_current;
	            //     // Also save the next belief vector in case this is the last
	            //     // step
	            //     prob_vec_evolution(tdx+1,:) = prob_vec_next;
	            //     // Save the true message as well
	            //     W_true = W;
	            // }
	            
	            // Determine if the receiver can stop (if message is sufficiently
	            // likely)
	            //[prob_max, W_hat] = max(prob_vec_next);
            
            	bool done = false;
			   	//printf("\n");
			   	for (unsigned long long i=0; i<M; ++i){
			   		//printf("%4.4f  ", prob_vec_next[i]);
					// if (prob_vec_next[i] > prob_max){
					// 	prob_max = prob_vec_next[i];
					// 	W_hat = i;
					// }
					if (prob_vec_next[i] > (1 - epsi)){
						// prob_max = prob_vec_next[i];
						W_hat = i;
						done = true;
						break;
					}
			    }
			    if(done)
			    	break;
			    // printf("\n%4.4f  ", prob_max);
	            // If probability of error is less than epsilon
	            // if (prob_max > (1 - epsi)){
	                
	            //     // Store the stopping time of this message
	            //     // tau = (tdx + 1);
	            //     // Don't send any more bits
	            //     break;
	            // }
	            
	            // // If probability of error is less than epsilon
	            // if (prob_vec_next[W] > (1 - epsi)){
	                
	            //     // Store the stopping time of this message
	            //     tau = (tdx + 1);
	            //     // Don't send any more bits
	            //     break;
	            // }
	            
	                        
	            // Store belief vector for next time
	            //prob_vec_current = prob_vec_next;    
			  //  	for (unsigned long long i=0; i<M; ++i){
					// prob_vec_current[i] = prob_vec_next[i];
			  //   }        
			    // memcpy( prob_vec_current, prob_vec_next, M );
			    // Swap arrays using pointers
			    temp = prob_vec_current;
			    prob_vec_current = prob_vec_next;
			    prob_vec_next = prob_vec_current;
	            
	        } // time tdx
	        
	        // if (tau == 0){
	        //     // We didn't find a stopping time before reaching the max number of
	        //     // bits
	        //     //error('Need to increase the max number of transmitted bits per message.');
	        //     printf("\n\n\n\n\nNeed to increase the max number of transmitted bits per message.\n\n\n\n\n");
	        // }
	        
	        // Update the number of message errors
	        if (W != W_hat) {
	            num_msg_errors = num_msg_errors + 1;
	        }
	        
	        // Update the cumulative stopping time (which will be averaged
	        // later)
	        sum_tau += tau;

	        
	    } // message msg_idx

        // delete [] prob_vec_next;
        // delete [] prob_vec_current;
	}    
	// elseif coding_scheme == GHBZ_determ_scheme || coding_scheme == GHBZ_random_scheme
	    
	//// Compute statistics of simulation

	// Compute probability of undetected error
	long double PrErr = (long double) (num_msg_errors) / num_msgs;


    if (print_output >= 1)
    {
		printf("\n--- M = %lu (k = %u) ---", M, (unsigned int) log2(M));
		printf("\nBSC crossover probability = %01.3f, target epsilon = %01.2e", del, epsi);
		printf("\n#{Msg. Errors} = %lu, out of %lu codewords simulated", num_msg_errors, num_msgs);
		printf("\nActual word-error probability is %01.3e", PrErr);

		// Compute the average rate
		long double Etau_sim = (long double) (sum_tau) / num_msgs;
		long double rate_sim = log2(M) / Etau_sim;
		printf("\nRate: %01.3f, latency: %01.2f", rate_sim, Etau_sim);
		// printf("\n");
	}
	// printf("\nsum_tau: %lu, num_msgs: %lu", sum_tau, num_msgs);
	*sum_tau_out = sum_tau;
	*_num_msg_errors = num_msg_errors;
}



// //function [Etau_sim, num_msg_errors, num_msgs] = active_coding_M(M, del, epsi, num_msgs, coding_scheme)
// void active_coding_M_thread(void *active_coding_M_args)
// // void active_coding_M_thread(unsigned long M, double del, double epsi, unsigned long num_msgs, int coding_scheme)
// {

// 	// function [Etau_sim, num_msg_errors, num_msgs] = active_coding_M(M, del, epsi, num_msgs, coding_scheme)
// 	//
// 	//   Sequential coding scheme for this particular number of messages M, BSC 
// 	//   crossover probability del, error probability epsi, and number of
// 	//   messages to simulate num_msgs.
// 	//
// 	//   coding_scheme is the particular coding scheme used
// 	//
// 	//   [2] Naghshvar & Javidi, "Optimal reliability over a class of
// 	//   binary-input channels with feedback", ITW 2012
// 	//   
// 	//   Includes an implementation of Algorithm 2 of [4] in partition_beliefs.m
// 	//   [4] Naghshvar & Javidi, "Extrinic Jensen-Shannon (EJS) Divergence", arXiv
// 	//   2013

// 	//// Simulation parameters

// 	// BSC channel capacity
// 	double capacity = 1 - binary_entropy(del); 

// 	// Debugging switch
// 	bool DBG_ON = false; 
// 	// DBG_ON = 1; 

// 	// Max. number of noise samples per message
// 	int num_bits_per_msg = ceil(log2(M)/capacity * 20);
// 	// printf("\ncapacity: %4.4f", capacity);
// 	// printf("\nnum_msgs: %u", num_msgs);
// 	// printf("\nnum_bits_per_msg: %u", num_bits_per_msg);

// 	// Initialize the number of message errors
// 	unsigned long num_msg_errors = 0;

// 	// Scheme from [2] & [4] for binary channels, based on partitioning messages
// 	// into two sets based on posterior probabilities (beliefs)
// 	int belief_sets_scheme = 1;
// 	// Generalized Horstein-Burnashev-Zigangirov scheme (deterministic)
// 	int GHBZ_determ_scheme = 2;
// 	// Generalized Horstein-Burnashev-Zigangirov scheme (randomized)
// 	int GHBZ_random_scheme = 3;

// 	//// Iterate through each of the messages to be sent
// 	// (so we don't have to allocate such a huge amount of memory for all the
// 	// arrays)
	       
// 	printf("\n");
// 	printf("\ncapacity: %4.4f",capacity);
// 	printf("\nnum_bits_per_msg: %d",num_bits_per_msg);
// 	//// Run the simulation
// 	// Previously mdx was the message index from 1 to num_msgs, but now it is a
// 	// constant
// 	unsigned long mdx = 1;
// 	unsigned long sum_tau = 0;

// 	if (coding_scheme == belief_sets_scheme) {	    
// 	    ////
// 	    // Deterministic message-partitioning scheme based on beliefs
// 	    //
// 	    double rand01;
// 	    double prob_vec_next[M], prob_vec_current[M];
// 	    unsigned long W, W_hat;
// 	        // Decoded message \hat W
// 	    unsigned long tau;
// 	    int x_vec[M];
// 	    int tx_bit, rx_bit;

// 	    for (unsigned long msg_idx = 1; msg_idx <= num_msgs; ++msg_idx) {	        
// 	        // Randomly assign the true message W
// 	        // unsigned long W;
// 	        if (DBG_ON) {
// 	            // Dummy message (always 1)
// 	            W = 1;
// 	        } else {
// 	            // Random assignment
// 	            rand01 = rand() / (double(RAND_MAX) + 0.01);
// 	            W = floor(rand01*M);
// 	        }
// 	        // Stopping times
// 	        tau = 0;	        
	        
// 	        // Start with a new message: initialize posterior probabilities at the
// 	        // receiver
// 	        for(int i=0;i<M;++i){
// 	        	prob_vec_current[i] = 1.0/M;
// 	        }
	        
// 	        //printf("... simulating message %u\n", msg_idx);
// 	        if ( msg_idx % 10000 == 0)
// 	            printf("\n... simulating message %lu", msg_idx);
	        
// 	        // Start sending symbols across the channel (tdx is the transmission
// 	        // index)
// 	        for(int tdx = 0; tdx < num_bits_per_msg; ++tdx){
	            
// 	            // Partition messages into two sets based on the current beliefs
// 	            // int x_vec[M];
// 	            //[set_S0, set_S1] = partition_beliefs(prob_vec_current, M, set_num);
// 				partition_beliefs(prob_vec_current, M, x_vec);
// 	            // Decide what to transmit in this time slot
// 	            // Select the bit corresponding to the true message
// 	            // (Sends a 0 if the true message is in set S_0,
// 	            //  Send a 1 if the true message is in set S_1)
// 	            tx_bit = x_vec[W];
	            
// 	            // Corrupt the bit by noise (representing reception across the noisy
// 	            // channel)
// 	            // int rx_bit = (tx_bit + noise_vec[tdx]) % 2;

// 				rand01 = rand() / double(RAND_MAX);				
// 	            rx_bit = (rand01 < del) ? !tx_bit : tx_bit;

// 	            // update belief vector (posterior probabilities) according to eq. in App.
// 	            // A2 of [2]
// 	            // The receiver knows the encoding rule, so it can compute what the
// 	            // Tx outputs would be for each possible message
	            
// 	            update_probs(prob_vec_current, rx_bit, x_vec, del, M, prob_vec_next);	            
            
//             	double prob_max = -1;
// 			   	//printf("\n");
// 			   	for (int i=0; i<M; ++i){
// 			   		//printf("%4.4f  ", prob_vec_next[i]);
// 					if (prob_vec_next[i] > prob_max){
// 						prob_max = prob_vec_next[i];
// 						W_hat = i;
// 					}
// 			    }
// 	            // If probability of error is less than epsilon
// 	            if (prob_max > (1 - epsi)){	                
// 	                // Store the stopping time of this message
// 	                tau = (tdx + 1);
// 	                // Don't send any more bits
// 	                break;
// 	            }	            
	                        
// 	            // Store belief vector for next time
// 	            //prob_vec_current = prob_vec_next;    
// 			   	for (int i=0; i<M; ++i){
// 					prob_vec_current[i] = prob_vec_next[i];
// 			    }        
	            
// 	        } // time tdx
	        
// 	        if (tau == 0){
// 	            // We didn't find a stopping time before reaching the max number of
// 	            // bits
// 	            //error('Need to increase the max number of transmitted bits per message.');
// 	            printf("\n\n\n\n\nNeed to increase the max number of transmitted bits per message.\n\n\n\n\n");
// 	        }
	        
// 	        // Update the number of message errors
// 	        if (W != W_hat) {
// 	            num_msg_errors = num_msg_errors + 1;
// 	        }
	        
// 	        // Update the cumulative stopping time (which will be averaged
// 	        // later)
// 	        sum_tau += tau;	        
// 	    } // message msg_idx
// 	}    

// 	// Compute probability of undetected error
// 	double PrErr = double(num_msg_errors) / num_msgs;

// 	printf("\n--- M = %lu (k = %u) ---", M, (unsigned int) log2(M));
// 	printf("\nBSC crossover probability = %01.3f, target epsilon = %01.2e", del, epsi);
// 	printf("\n#{Msg. Errors} = %lu, out of %lu codewords simulated", num_msg_errors, num_msgs);
// 	printf("\nActual word-error probability is %01.3e", PrErr);

// 	// Compute the average rate
// 	double Etau_sim = double(sum_tau) / num_msgs;
// 	double rate_sim = log2(M) / Etau_sim;
// 	printf("\nsum_tau: %lu, num_msgs: %lu", sum_tau, num_msgs);
// 	printf("\nRate: %01.3f, latency: %01.2f", rate_sim, Etau_sim);
// 	// printf("\n");
// }


int main(int argc, char* argv[]) 
{
	unsigned int k;
	// int num_threads;
	unsigned long num_msgs;
	if (argc <2)
		k = 1;
	else
		k = (unsigned int) atoi(argv[1]);
	// if (argc <3)
	// 	num_threads = 1;
	// else
	// 	num_threads = atoi(argv[2]);
	if (argc <3)
		seed = 1;
	else
		seed = atoi(argv[2]);
	if (argc < 4)
		num_msgs = 200000;
	else
		num_msgs = atoi(argv[3]);


	// num_threads = 1;

	printf("\nk: %u", k);
	// printf("\nnum_threads: %d", num_threads);
	printf("\n");

    // Debugging switch
	bool DBG_ON = false;

	// overall probability of error desired, epsilon
	// double epsi = pow(10,-1);           // for PhD defense animation
	long double epsi = pow(10,-3);       // for PhD thesis

	// BSC crossover probability delta
	//del = 0.11;   // in [PPV 2011]
	//del = 0.11009;  // for SNR = 1.77 dB, hard-decision decoding
	long double del = 0.05005;  // for SNR = 4.32 dB, hard-decision decoding, // for PhD thesis
	//double del = 0.1;  // for PhD defense animation


	// BSC channel capacity
	long double capacity = 1 - binary_entropy(del); 
	long double lr = del / (1 - del);   // likelihood ratio
	long double llr = log2(lr);           // log likelihood ratio

	// Maximum KL divergence (relative entropy), C_1
	long double C1 = (2*del - 1) * llr;
	// Maximum difference between Ui(t), Ui(t+1) in (13) of [2]
	long double C2 = -llr;
	// F(C,C1,C2) = K', a constant independent of M & epsilon
	long double Kprime = 3 * C2*C2 / (capacity*C1);
	  
	// k = # message bits
	// k = [2:2:20];
	//k = [2:12];
	// k = [2:8];  // active_coding_snr4.32_belief_sets_20140204
	//k = [9:11];
	//unsigned int k = 3;

	// GOURAV EDIT
	// M = # messages
	unsigned long long  M = round(pow(2,k));
	// Update k after choosing M that is an integer
	//k = log2(M);




	// Average length for each M
	//double Etau_sim = 0;
	// Number of messages (codewords) to simulate and number of errors obtained 
	// for each M
	// unsigned long num_msgs;
	// if (DBG_ON) {
	//     //num_msgs = 10000 * ones(size(M));
	//     num_msgs = 1;
	// }else{
	// //     num_msgs = 10^4 * ones(size(M));
	//     num_msgs = 2 * pow(10,5);
	// //     num_msgs = 1 * 10^6 * ones(size(M));
	// }

	unsigned long num_msg_errors = 0;

	// Possible sequential coding schemes

	// Scheme from [2] & [4] for binary channels, based on partitioning messages
	// into two sets based on posterior probabilities (beliefs)
	int belief_sets_scheme = 1;     
	// Generalized Horstein-Burnashev-Zigangirov scheme (deterministic)
	int GHBZ_determ_scheme = 2;
	// Generalized Horstein-Burnashev-Zigangirov scheme (randomized)
	int GHBZ_random_scheme = 3;

	// Pick the coding scheme to use
	int coding_scheme = belief_sets_scheme;
	// coding_scheme = GHBZ_determ_scheme;
	// coding_scheme = GHBZ_random_scheme;

	
	printf("\nCapacity: %4.4Lf", capacity);
	printf("\nC1: %4.4Lf", C1);
	printf("\nC2: %4.4Lf", C2);
	printf("\nKprime: %4.4Lf", Kprime);
	printf("\nk: %d", k);

	printf("\n");
	printf("\nM: %lu", M);
	printf("\ndel: %4.4Lf", del);
	printf("\nepsi: %4.4Lf", epsi);
	printf("\nnum_msgs: %lu", num_msgs);
	printf("\ncoding_scheme: %d", coding_scheme);
	//[Etau_sim(Mdx), num_msg_errors(Mdx)] = active_coding_M(M(Mdx), del, epsi, num_msgs(Mdx), coding_scheme);



  // code_to_time();

  // clock_t end = clock();

  	clock_t begin = clock();

	// active_coding_M(M, del, epsi, num_msgs, coding_scheme);

	unsigned long total_sum_tau = 0;
	int print_output = 0;
	active_coding_M(M, del, epsi, num_msgs, coding_scheme, &total_sum_tau, &num_msg_errors, print_output);
    
 //    //int num_threads = 4;

	// std::thread t[num_threads];
	// int print_output = 0;

	// unsigned long sum_taus[num_threads];

	// unsigned long num_msgs_per_thread = floor(num_msgs/num_threads);

	// //Launch a group of threads
	// for (int i = 0; i < num_threads; ++i) {
	// 	t[i] = std::thread(active_coding_M, M, del, epsi, num_msgs_per_thread, coding_scheme, &sum_taus[i], print_output);
	// }

	// // std::cout << "Launched from the main\n";

	// //Join the threads with the main thread
	// for (int i = 0; i < num_threads; ++i) {
	// 	t[i].join();
	// }
 
	// unsigned long total_sum_tau = 0;
	// for(int i = 0; i<num_threads; ++i){
	// 	total_sum_tau += sum_taus[i];
	// }

	// unsigned long total_num_msgs = num_msgs_per_thread*num_threads;
	unsigned long total_num_msgs = num_msgs;
	// Compute the average rate
	long double Etau_sim = (long double) (total_sum_tau) / total_num_msgs;
	long double rate_sim = log2(M) / Etau_sim;
	printf("\nsum_tau: %lu, num_msgs: %lu", total_sum_tau, total_num_msgs);
	printf("\nRate: %01.3Lf, latency: %01.2Lf", rate_sim, Etau_sim);
	printf("\nNumber of Message Errors: %lu", num_msg_errors);


	clock_t end = clock();
  	long double elapsed_secs = (long double) (end - begin) / CLOCKS_PER_SEC;
  	printf("\nElapsed Time: %0.0Lf\n", elapsed_secs);


	ofstream myfile ("./data/ASHT_k-" + to_string((long long unsigned int) k) + "_seed-" + to_string((long long unsigned int) seed) + ".txt");
	if (myfile.is_open())
	{
	myfile << "sum_tau\t" << total_sum_tau << endl;
	myfile << "num_msgs\t" << total_num_msgs << endl;
	myfile << "k\t" << k << endl;
	myfile << "Elapsed_Time\t" << elapsed_secs << endl;
	myfile << "Num_Msg_Errs\t" << num_msg_errors << endl;
	myfile.close();
	}
	else printf("Unable to open file");

    return 0;
}