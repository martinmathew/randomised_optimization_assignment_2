import pandas as pd
import matplotlib.pyplot as plt

def plot_graph(filename, grapgname):
    fig, ax = plt.subplots()
    df = pd.read_csv(filename)
    df = df.T

    df.plot(ax = ax)
    ax.legend(["RHC", "SA", "GA", "MIMIC"]);
    plt.savefig(grapgname)

def plot_graph_single_opt(filename, grapgname, opt_algo, x_label,legend_title):
    fig, ax = plt.subplots()
    df = pd.read_csv(filename)
    df = df.T

    df.plot(ax = ax)
    ax.legend(opt_algo,title=legend_title);
    ax.set_xlabel(x_label)
    plt.savefig(grapgname)



plot_graph_single_opt('output/travelling_salesman_result_mimic_keep_count.csv','graphs/mimic_keep_count_itr.png' , [10,100, 200, 300, 400, 500, 600], "Iterations","Keep Count")
plot_graph_single_opt('output/four_peak_result_mimic_keep_count.csv','graphs/mimic_keep_count_itr_four_peaks.png' , [10,100, 200, 300, 400, 500, 600], "Iterations","Keep Count")
plot_graph_single_opt('output/n_queens_result_mimic_keep_count.csv','graphs/mimic_keep_count_itr_n_queen.png' , [10,100, 200, 300, 400, 500, 600], "Iterations","Keep Count")





plot_graph_single_opt('output/travelling_salesman_result_ga_cross_over_mutate.csv','graphs/ga_cross_mutant.png' , [20,40,60,80,100,120,140,160,180, 200], "Mutations","Crossover")
plot_graph_single_opt('output/four_peaks_result_ga_cross_over_mutate.csv','graphs/four_peaks_ga_cross_mutant.png' , [20,40,60,80,100,120,140,160,180, 200], "Mutations","Crossover")
plot_graph_single_opt('output/n_queens_result_ga_cross_over_mutate.csv','graphs/n_queens_ga_cross_mutant.png' , [20,40,60,80,100,120,140,160,180, 200], "Mutations","Crossover")
#


plot_graph_single_opt('output/travelling_salesman_result_sa_coolin_fitness.csv','graphs/sa_cooling_tsp_fitness.png' , [0.1, 0.3, 0.5, 0.7, 0.9, 1.0],"Iterations","Cooling Exponent")

plot_graph_single_opt('output/four_peaks_result_sa_coolin_fitness.csv','graphs/sa_cooling_foour_peaks_fitness.png' , [0.1, 0.3, 0.5, 0.7, 0.9, 1.0],"Iterations","Cooling Exponent" )
plot_graph_single_opt('output/n_queens_result_sa_coolin_fitness.csv','graphs/sa_cooling_n_queens_fitness.png' , [0.1, 0.3, 0.5, 0.7, 0.9, 1.0], "Iterations","Cooling Exponent" )


plot_graph_single_opt('output/four_peaks_result_rhc_fitness.csv','graphs/rhc_fourpeak_fitness.png' , ['RHC'],"Iterations", "Algorithim" )
plot_graph_single_opt('output/nquuen_result_rhc_fitness.csv','graphs/rhc_n_queen_fitness.png' , ['RHC'],"Iterations", "Algorithim" )




plot_graph_single_opt('output/travelling_salesman_result_runtime.csv', 'graphs/tsp_runtime.png', ['RHC', 'SA', 'GA', 'MIMIC'], 'Iterations', 'Algorithim')
plot_graph_single_opt('output/travelling_salesman_result_runtime.csv', 'graphs/tsp_runtime.png', ['RHC', 'SA', 'GA', 'MIMIC'], 'Iterations', 'Algorithim')

plot_graph_single_opt('output/fourpeak_result_accuracy.csv', 'graphs/four_square_accuract.png', ['RHC', 'SA', 'GA', 'MIMIC'], 'Iterations', 'Algorithim')
plot_graph_single_opt('output/fourpeak_result_runtime.csv', 'graphs/four_square_runtime.png', ['RHC', 'SA', 'GA', 'MIMIC'], 'Iterations', 'Algorithim')










