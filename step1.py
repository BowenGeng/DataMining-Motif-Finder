from random import *
import copy
import os
# from sympy import *
# from math import *


def step1_function(params):
    '''============ Input ============='''

    ICPC = params['ICPC']     # Default: 2,   ICPC = 1, 1.5, 2    -information content per column
    ML = params['ML']         # Default: 8,   ML = 6, 7, 8        -motif length
    SC = params['SC']         # Default: 10,  SC = 5, 10, 20      -sequence count
    SL = params['SL']         # Default: 500                      -sequence length
    DNA_dict = ['A', 'C', 'G', 'T']

    ICPC_list = {1: [0.5, 0.5, 0, 0], 2: [1, 0, 0, 0], 1.5: [0.02, 0.01, 0.94, 0.03]}

    output_sequence = [[] for x in range(SC)]
    output_site_location = []
    output_site = []
    output_motif = []
    output_motif_length = 0

    path = 'data/'
    motif_name = 'MOTIF1'

    '''============ step2 ============='''
    #  Generate SC random sequences (with uniform frequencies of A,C,G,T). Each random sequence has length SL.

    for i in range(SC):
        for j in range(SL):
            output_sequence[i].append(DNA_dict[randint(0, 3)])

    '''============ step3 ============='''
    #  Generate a random motif (position weight matrix) of length ML, with total information content being ICPC * ML.
    '''
    # Solution 1
    # use Sympy to solve the formula

    for i in range(ML):
        x1 = random()
        x2 = uniform(0, 1-x1)
        x3 = Symbol('x3')
        res = solve([x1*log(4*x1)+x2*log(4*x2)+x3*log(4*x3)+(1-x1-x2-x3)*log(4*(1-x1-x2-x3))-ICPC], [x3])
        output_motif[i].append(x1)
        output_motif[i].append(x2)
        output_motif[i].append(res['x3'])
        output_motif[i].append(1-x1-x2-res['x3'])
    '''

    # Solution 2
    temp_motif = ICPC_list[ICPC]

    for i in range(ML):
        shuffle(temp_motif)
        output_motif.append(copy.copy(temp_motif))

    output_motif_length = ML

    '''============ step4 ============='''
    # Generate SC strings of length ML each, by sampling from this random motif.

    for i in range(SC):
        temp_site = []
        for j in range(ML):
            temp_motif = output_motif[j]
            thres_A = temp_motif[0]
            thres_C = thres_A + temp_motif[1]
            thres_G = thres_C + temp_motif[2]

            temp_rand = random()

            if temp_rand <= thres_A:
                temp_site.append('A')
            elif temp_rand <= thres_C:
                temp_site.append('C')
            elif temp_rand <= thres_G:
                temp_site.append('G')
            else:
                temp_site.append('T')

        output_site.append(temp_site)

    '''============ step5 ============='''
    # "Plant" one sampled site at a random location in each random sequence generated in step 2.

    for i in range(SC):
        temp_location = randint(0, SL - ML)
        output_site_location.append(temp_location)
        for j in range(ML):
            output_sequence[i][temp_location + j] = output_site[i][j]

    '''============ step6 ============='''
    #  Write out the SC sequences into a FASTA format file called "sequences.fa"
    temp_file = open(path+'sequences.fa', 'w')

    for i in range(SC):
        temp_file.write('>' + str(i) + '\t' + str(SL) + '\n')
        for j in range(SL):
            temp_file.write(str(output_sequence[i][j]))
        #temp_file.write(str('\n<\n'))
        temp_file.write('\n\n')

    temp_file.flush()
    os.fsync(temp_file)
    temp_file.close()

    '''============ step7 ============='''
    #  In a separate text file (called "sites.txt") write down the location of the planted site in each sequence.

    temp_file = open(path + 'sites.txt', 'w')
    # temp_file.close()
    # temp_file = open(path + 'sites.txt', 'a')

    for i in range(len(output_site_location)):
        temp_file.write(str(output_site_location[i]) + ' ')

    # print "======== Step 1 generated sites:"
    # print output_site_location

    temp_file.flush()
    os.fsync(temp_file)
    temp_file.close()

    '''============ step8 ============='''
    # In a separate text file (called "motif.txt") write down the motif that was generated in step 3.

    temp_file = open(path + 'motif.txt', 'w')
    temp_file.write('>' + motif_name + '\t' + str(ML) + '\n')

    for i in range(ML):
        for j in range(len(output_motif[0])):
            temp_file.write(str(output_motif[i][j]) + '\t')
        temp_file.write('\n')

    temp_file.write('<')
    temp_file.flush()
    os.fsync(temp_file)
    temp_file.close()

    '''============ step9 ============='''
    # In a separate test file (called "motiflength.txt") write down the motif length.

    temp_file = open(path + 'motiflength.txt', 'w')
    temp_file.close()
    temp_file = open(path + 'motiflength.txt', 'a')

    temp_file.write(str(ML))
    temp_file.flush()
    temp_file.close()
