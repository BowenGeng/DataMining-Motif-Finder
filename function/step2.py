import os
from sets import Set
from math import *
import time

global_IC_MAX = -1
predicted_sites_list_cands = []
predicted_PWM_list_cands = []

def step2_function():
    start_time = time.time()
    data_path = './data/'
    predicted_data_path = './predicted_data/'
    sequence_file = 'sequences.fa'
    motiflength_file = 'motiflength.txt'
    sequences_list = []
    ML = 0
    SL = 0

    global global_IC_MAX
    global predicted_sites_list_cands
    global predicted_PWM_list_cands
    global_IC_MAX = -1
    predicted_sites_list_cands = []
    predicted_PWM_list_cands = []
    SC = 0


    '''============ Input From Files ============='''
    with open(data_path + sequence_file, "r", 0) as f:
        line_count = 0
        for line in f.readlines():
            if line_count % 3 == 1 and len(line) > 0:
                sequences_list.append(line.rstrip())
            line_count += 1

    with open(data_path + motiflength_file, "r", 0) as f:
        for line in f.readlines():
            ML = int(line.rstrip())
            break

    SC = len(sequences_list) # 5, 10 ,20
    SL = len(sequences_list[0]) # 500


    '''============ Helper Functions ============='''
    # Calculating Information Content(IC)
    # PWM_list: [{'A': 2, 'C': 3, 'G': 1, 'T': 4}, {....}, ... , {....}]

    def list_to_PWM_array(PWM_list):
        PWM = [[] for i in range(len(PWM_list))]
        keys = ["A", "C", "G", "T"]
        for i in range(len(PWM_list)):
            elem_dict = PWM_list[i]
            total_num = 0
            for key in keys:
                total_num += elem_dict[key]
            for key in keys:
                p = elem_dict[key] * 1.0 / total_num
                PWM[i].append(p)
        return PWM

    def calculate_IC(PWM_list):
        IC = 0
        for elem_dict in PWM_list:
            total_num = 0
            for key in elem_dict:
                total_num += elem_dict[key]
            for key in elem_dict:
                p = elem_dict[key] * 1.0 / total_num
                if p == 0: # p * log(p * 4) == 0
                    IC += 0
                else:
                    IC += p * log(p * 4.0)
        return IC


    def copy_PWM_list(PWM_list):
        new_PWM = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for i in range(ML)]
        for n in range(ML):
            elem_dict = PWM_list[n]
            for key in elem_dict:
                new_PWM[n][key] = elem_dict[key]
        return new_PWM


    def copy_list(list):
        new_list = []
        for elem in list:
            new_list.append(elem)
        return new_list


    # DFS recusion for searching all possible optimal results
    def find_PWN_cands(offset, curr_sites, curr_PWM_list, first, second):
        global global_IC_MAX
        global predicted_sites_list_cands
        global predicted_PWM_list_cands

        # Visited before, skip and DFS to the next level
        if offset == first or offset == second:
            find_PWN_cands(offset + 1, curr_sites, curr_PWM_list, first, second)
            return

        # End level of DFS
        # print "SC when DFS end: "
        # print SC

        if offset == SC:
            IC = calculate_IC(curr_PWM_list)
            # print curr_sites
            if IC > global_IC_MAX:
                global_IC_MAX = IC
                predicted_sites_list_cands = [curr_sites]
                predicted_PWM_list_cands = [curr_PWM_list]
            elif IC == global_IC_MAX:
                predicted_sites_list_cands.append(curr_sites)
                predicted_PWM_list_cands.append(curr_PWM_list)
            return

        local_IC_MAX = -1
        local_PWM_list_cands = []
        local_sites_cands = []

        for m in range(0, SL - ML + 1):
            curr_seq = sequences_list[offset][m: m + ML]

            # Copy the PWM
            new_PWM = copy_PWM_list(curr_PWM_list)

            # Load the curr_seq, update new_PWM
            for n in range(ML):
                new_PWM[n][curr_seq[n]] += 1

            # Calculate IC
            curr_IC = calculate_IC(new_PWM)

            # Update local optimals
            if curr_IC > local_IC_MAX:
                local_IC_MAX = curr_IC
                new_sites = copy_list(curr_sites)
                new_sites.append((offset, m))
                local_sites_cands = [new_sites]
                local_PWM_list_cands = [new_PWM]
            elif curr_IC == local_IC_MAX:
                new_sites = copy_list(curr_sites)
                new_sites.append((offset, m))
                local_sites_cands.append(new_sites)
                local_PWM_list_cands.append(new_PWM)

        # DFS to the next level
        for i in range(len(local_sites_cands)):
            find_PWN_cands(offset + 1, local_sites_cands[i], local_PWM_list_cands[i], first, second)



    '''============ Main Function ============='''

    # initialization: match and mismatch between the 1st and 2nd sequence
    # print "SC before DFS: "
    # print SC

    for first in range(SC):
        for second in range(first + 1, SC):
            max_match_count = -1

            for m in range(0, SL - ML + 1):
                for n in range(0, SL - ML + 1):
                    match_count = 0
                    seq1 = sequences_list[first][m: m + ML]
                    seq2 = sequences_list[second][n: n + ML]
                    for k in range(ML):
                        if seq1[k] == seq2[k]:
                            match_count += 1

                    if match_count > max_match_count:
                        max_match_count = match_count

            for m in range(0, SL - ML + 1):
                for n in range(0, SL - ML + 1):
                    match_count = 0
                    seq1 = sequences_list[first][m: m + ML]
                    seq2 = sequences_list[second][n: n + ML]
                    for k in range(ML):
                        if seq1[k] == seq2[k]:
                            match_count += 1

                    if match_count == max_match_count:
                        new_PWM_list = [{'A': 0, 'C': 0, 'G': 0, 'T': 0} for i in range(ML)]
                        for k in range(0, ML):
                            new_PWM_list[k][seq1[k]] += 1
                            new_PWM_list[k][seq2[k]] += 1
                        local_sites = [(first, m), (second, n)]
                        find_PWN_cands(0, local_sites, new_PWM_list, first, second)


    '''============ Print Results ============='''
    elapsed_time = time.time() - start_time
    print elapsed_time

    def print_sites(sites_list):
        list.sort(sites_list, key = lambda x : x[0])
        predicted_sites_output = []
        for elem in sites_list:
            predicted_sites_output.append(elem[1])
        print "predicted_sites_output: "
        print predicted_sites_output

    def print_PWM(PWM_list):
        predicted_PWM_output = [[] for i in range(ML)]
        for n in range(ML):
            elem_dict = PWM_list[n]
            total_num = 0
            for key in elem_dict:
                total_num += elem_dict[key]
            for key in elem_dict:
                p = elem_dict[key] * 1.0 / total_num
                predicted_PWM_output[n].append(float("{:10.4f}".format(p)))
        print "predicted_PWM_output: "
        print predicted_PWM_output

    # Debug: print results
    # for i in range(len(predicted_sites_list_cands)):
    #     print "sites candidates: "
    #     print_sites(predicted_sites_list_cands[i])
    #
    #     print "PWM candidates: "
    #     print print_PWM(predicted_PWM_list_cands[i])
    #
    # print "len(predicted_sites_list_cands): "+ str(len(predicted_sites_list_cands))


    '''============ write predict sites to file ============='''
    output_site_location = predicted_sites_list_cands[0]
    list.sort(output_site_location, key = lambda x : x[0])
    temp_file_sites = open(predicted_data_path+'predictedsites.txt', 'w')
    # temp_file.close()
    # temp_file = open(predicted_data_path+'predictedsites.txt', 'a', 0)

    for i in range(len(output_site_location)):
        temp_file_sites.write(str(output_site_location[i][1]) + ' ')

    temp_file_sites.flush()
    os.fsync(temp_file_sites)
    temp_file_sites.close()


    # print "======== Predicted sites output: "
    # print output_site_location

    '''============ write predict motif to file  ============='''
    output_motif = list_to_PWM_array(predicted_PWM_list_cands[0])
    temp_file_motif = open(predicted_data_path+'predictedmotif.txt', 'w')
    # temp_file.close()
    # temp_file = open(predicted_data_path+'predictedmotif.txt', 'a', 0)

    temp_file_motif.write('>' + "Predicted Motif" + '\t' + str(ML) + '\n')

    for i in range(ML):
        for j in range(len(output_motif[0])):
            temp_file_motif.write(str(output_motif[i][j]) + '\t')
        temp_file_motif.write('\n')

    temp_file_motif.flush()
    os.fsync(temp_file_motif)
    temp_file_motif.close()
    return elapsed_time
