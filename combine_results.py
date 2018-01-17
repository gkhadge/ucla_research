from __future__ import division
import csv
import glob
import numpy as np
import sys
k = int(sys.argv[1])

##filename = "./data/ASHT_k-5_seed-2.txt"
print "running"

sum_tau = 0
num_msgs = 0
taus = np.array([], dtype=np.uint32)
elapsed_times = np.array([], dtype=np.float)

for filename in glob.glob('./data/ASHT_k-' + str(k) + '_seed-*.txt'):
##    print filename
    with open(filename,'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        data = list(reader)
##        sum_tau = sum_tau + int(data[0][1])
        num_msgs = num_msgs + int(data[1][1])
        taus = np.append(taus, [int(data[0][1])])
        elapsed_times = np.append(elapsed_times, [float(data[3][1])])

##print sum_tau
##print taus
sum_tau = np.sum(taus);
tau_std = np.std(taus);
print "Sum Tau: ", sum_tau
print "Num Msgs: ", num_msgs
print "Tau st. dev: ", tau_std
print "Elapsed Time (avg): ", np.mean(elapsed_times)
Etau_sim = sum_tau / num_msgs;
rate_sim = k / Etau_sim;
print "Rate: ", rate_sim
print "Latency: ", Etau_sim
