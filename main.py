import numpy as np
import os
import time
from function.step1 import step1_function
from function.step2 import step2_function
from function.step3 import step3_calculate_metrics
from function.step3 import step3_plot_and_save

params_list = {'ICPC': [1, 1.5, 2], 'ML': [6, 7, 8], 'SC': [5, 10, 20], 'SL': 500}
random_sample_times = 10

def run(variant):
    params_input = {'ICPC': 2, 'ML': 8, 'SC': 10, 'SL': 500} # default params
    params_variants = params_list[variant]  # find the changeable param, get the list of value for traversal
    results_global = {                      #
        'variant': variant,
        'avg': {
            'KLD': [], #[(1, 10), (1.5, 20), (2, 20)]
            'ovlp_sites': [], #[(1, 10), (1.5, 20), (2, 20)]
            'ovlp_positions': [], #[(1, 10), (1.5, 20), (2, 20)]
            'running_time': []
        },
        'std': {
            'KLD': [], #[(1, 10), (1.5, 20), (2, 20)]
            'ovlp_sites': [], #[(1, 10), (1.5, 20), (2, 20)]
            'ovlp_positions': [], #[(1, 10), (1.5, 20), (2, 20)]
            'running_time': []
        }
    }

    for param_value in params_variants:
        params_input[variant] = param_value
        results_10_runs = {
                'KLD': [],
                'ovlp_sites': [], #[6,8,7,6,8,7,7,8,7,6]
                'ovlp_positions': [], #[6,8,7,6,8,7,7,8,7,6]
                'running_time': []
        }
        metrics = ['KLD', 'ovlp_sites', 'ovlp_positions', 'running_time']

        # random sampling for 10 times
        for i in range(random_sample_times):
            step1_function(params_input) # generate and write to motif.txt, sites.txt, sequences.fa
            running_time = step2_function() # 1. read from step 1 file) # 2. output to step 2 files: predictedmotif.txt and predictedsites.txt
            result_single_run = step3_calculate_metrics()# result_single_run: {
                                                         #    'KLD': 8.2,
                                                         #   'ovlp_sites': 2,
                                                         #   'ovlp_positions': 3
                                                         # }
            for metric_attr in metrics:
                if metric_attr == "running_time":
                    results_10_runs[metric_attr].append(running_time)
                    results_10_runs[metric_attr].append(running_time)
                else:
                    value = result_single_run[metric_attr]
                    results_10_runs[metric_attr].append(value)
                    results_10_runs[metric_attr].append(value)

        # calculate the avg and std based on 10 times' sampling
        # add to results_global for plotting
        print "========== RESULT AFTER 10 RUNS ==========="
        print variant + ": " + str(param_value)
        for metric_attr in metrics:
            curr_metric_mean = np.mean(results_10_runs[metric_attr])
            curr_metric_std = np.std(results_10_runs[metric_attr])
            results_global['avg'][metric_attr].append((param_value, curr_metric_mean))
            results_global['std'][metric_attr].append((param_value, curr_metric_std))

            #print avg, std results after 10 runs
            print "== AVG =="
            print "  " + metric_attr + ": " + str(curr_metric_mean)
            print "== STD =="
            print "  " + metric_attr + ": " + str(curr_metric_std)
            print


    # prepare the data in results_global for plotting the chart
    # form of "x:[], y:[]" for mathplotlib
    for metric_attr in metrics:
        curr_metric_avg_points = results_global['avg'][metric_attr]
        curr_metric_new = {'x': [], 'y': []}
        for point in curr_metric_avg_points:
            curr_metric_new['x'].append(point[0])
            curr_metric_new['y'].append(point[1])
        results_global['avg'][metric_attr] = curr_metric_new


    for metric_attr in metrics:
        curr_metric_std_points = results_global['std'][metric_attr]
        curr_metric_new = {'x': [], 'y': []}
        for point in curr_metric_std_points:
            curr_metric_new['x'].append(point[0])
            curr_metric_new['y'].append(point[1])
        results_global['std'][metric_attr] = curr_metric_new
        # results_global = {
        #     'variant': variant,
        #     'avg': {
        #         'KLD': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
        #         'ovlp_sites': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
        #         'ovlp_positions': {'x': [1, 1.5, 2], 'y': [10, 20, 30]}
        #         'running_time': {'x': [1, 1.5, 2], 'y': [10, 20, 30]}
        #     },
        #     'std': {
        #         'KLD': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
        #         'ovlp_sites': {'x': [1, 1.5, 2], 'y': [10, 20, 30]},
        #         'ovlp_positions': {'x': [1, 1.5, 2], 'y': [10, 20, 30]}
        #         'running_time': {'x': [1, 1.5, 2], 'y': [10, 20, 30]}
        #     }
        # }
    print "results_global: "
    print results_global
    step3_plot_and_save(results_global)

# Main Function
variants = ['SC']
for variant in variants:
    run(variant)
    break
