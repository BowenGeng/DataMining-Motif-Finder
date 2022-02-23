import numpy as np
import matplotlib.pyplot as plt


# Calculate KL-Divergence, position and sites matches for one run
def step3_calculate_metrics():

    motif = []
    predictedmotif = []
    sites_list = []
    predictedsites_list = []
    sequences = []

    #read input from motif.txt
    with open('./data/motif.txt','rt') as mymotif:
        for myline in mymotif:
            if myline[0] != '>' and myline[0] != '<':
                temp = [float(i) for i in myline.split()]
                motif.append(temp)

    #read input from predictedmotif.txt
    with open('./predicted_data/predictedmotif.txt','rt') as mypredictedmotif:
        for myline in mypredictedmotif:
            if myline[0] != '>' and myline[0] != '<':
                temp = [float(i) for i in myline.split()]
                predictedmotif.append(temp)

    # print "predictedmotif: "
    # print predictedmotif

    #read input from sites.txt
    sites_file = open('./data/sites.txt', 'r', 0)
    sites = sites_file.read()
    for i in range(len(sites.split())):
        sites_list.append(int(sites.split()[i]))
    #print ites_list

    #read input from predictedsites.txt
    predictedsites_file = open('./predicted_data/predictedsites.txt','r', 0)
    predictedsites = predictedsites_file.read()
    for i in range(len(predictedsites.split())):
        predictedsites_list.append(int(predictedsites.split()[i]))
    ##print predictedsites_list

    #read input from sequences.fa
    with open('./data/sequences.fa','rt', 0) as f:
        line_count = 0
        for line in f.readlines():
            if line_count % 3 == 1 and len(line) > 0:
                sequences.append(line.rstrip())
            line_count += 1
    #print sequences

    ML = len(motif)
    SC = len(sites_list)
    # print SC
    KLD = 0.0
    for i in range(ML):
    #for i in range(ML):
        for j in range(4):
            if motif[i][j] == 0:
                motif[i][j] = 0.00001
            if predictedmotif[i][j] == 0:
                predictedmotif[i][j] = 0.00001
            KLD += np.log(motif[i][j] * 1.0 /predictedmotif[i][j]) * motif[i][j]

    norm_KLD = KLD * 1.0 /ML
    # print 'KL Divergence = ',norm_KLD

    overlap_position = sum(1 for i in range(SC) if abs(predictedsites_list[i]-sites_list[i])< ML * 1.0 / 2)
    # print 'Number of overlapping positions = ',overlap_position

    overlap_sites = 0
    for i in range(SC):
        site = sequences[i][sites_list[i]:sites_list[i]+ML]
        predictedsite = sequences[i][predictedsites_list[i]:predictedsites_list[i]+ML]
        temp = sum(1 for j in range(ML) if site[j] == predictedsite[j])
        if temp >= ML * 1.0 / 2:
            overlap_sites += 1

    # print 'Number of overlapping sites = ',overlap_sites

    cal_result = {'KLD':norm_KLD, 'ovlp_sites':overlap_sites, 'ovlp_positions':overlap_position}
    # print "====== step3.a results"
    # print cal_result
    print

    return cal_result


# Plot line chart after 10 runs
def step3_plot_and_save(results_global):
    x = results_global["variant"]
    for y in results_global["avg"]:
        y_list = results_global["avg"][y]['y']
        x_list = results_global["avg"][y]['x']

        plt.plot(x_list, y_list, '-ok')
        plt.xlabel(x)
        plt.ylabel(y)
        for a,b in zip(x_list, y_list):
            value_str = "{:10.3f}".format(b)
            plt.text(a, b, value_str)
        plt.savefig("avg_" + x +'_'+ y + ".png")
        plt.close()

    for y in results_global["std"]:
        y_list = results_global["std"][y]['y']
        x_list = results_global["std"][y]['x']

        plt.plot(x_list, y_list, '-ok')
        plt.xlabel(x)
        plt.ylabel(y)
        for a,b in zip(x_list, y_list):
            value_str = "{:10.3f}".format(b)
            plt.text(a, b, value_str)
        plt.savefig("std_" + x +'_'+ y + ".png")
        plt.close()


# variant = 'ICPC'
# results_global = {'variant': variant,'data': {'KL': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},'ovlp_sites': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
# results_global = {
#     'variant': 'ICPC',
#     'avg': {
#         'KLD': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
#         'ovlp_sites': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
#         'ovlp_positions': {'x': [1, 1.5, 2], 'y': [10, 20, 30]}
#     },
#     'std': {
#         'KLD': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
#         'ovlp_sites': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
#         'ovlp_positions': {'x': [1, 1.5, 2], 'y': [10, 20, 30]}
#     }
# }
#   # 'ovlp_positions': {'x': [1, 1.5, 2], 'y': [10, 20, 30]}}}
#
# step3_calculate_metrics()
# step3_plot_and_save(results_global)
