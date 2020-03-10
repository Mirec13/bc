import functions as f
from scipy.stats import chi2
import time
from collections import namedtuple


def search(file_name, initial_prob, number_of_iter, number_of_population, p_value, number_of_epi, repeat,
           min_value_for_df):
    # define structs and other variables
    non_dominated = []
    snp_data = namedtuple("snp_data", "sample_size snp_size data state")
    comb = pow(3, number_of_epi)
    flowers = []

    # init the population
    for i in range(number_of_population):
        flower = namedtuple("flower", "loci observed_value expected_value objective_function_score")
        flower.observed_value = [[0 for x in range(2)] for y in range(comb)]
        flower.expected_value = [0.0] * comb
        flower.loci = [0] * number_of_epi
        flower.objective_function_score = [0.0] * 2
        flowers.append(flower)

    while True:
        # Load the SNP data
        snp_data.data = f.load_file(file_name)

        # if we do not have any data end the loop
        if not snp_data.data:
            break

        for repeating in range(repeat):
            # initialize SNP data
            snp_data.snp_size = len(snp_data.data[0]) - 1
            snp_data.sample_size = len(snp_data.data)
            snp_data.state = []
            for row in snp_data.data:
                snp_data.state.append(row[snp_data.snp_size])
                del row[snp_data.snp_size]

            # initialize the first population
            vector = f.init_first_population(number_of_epi, snp_data.snp_size, number_of_population)
            print(vector)
            for i in range(number_of_population):
                for j in range(number_of_epi):
                    flowers[i].loci[j] = vector[i][j]
                f.compute_o_and_e_values(snp_data.data, snp_data.sample_size, flowers[i].loci, flowers[i].observed_value,
                                         flowers[i].expected_value, snp_data.state, number_of_epi)
                flowers[i].objective_function_score[0] = f.k2_score(flowers[i].observed_value, comb)
                flowers[i].objective_function_score[1] = f.gi_score(flowers[i].observed_value, comb,
                                                                    snp_data.sample_size)
            f.pareto_optimization(flowers, number_of_population, non_dominated)
            for i in flowers:
                print(i.objective_function_score)

            print("\n", non_dominated)
            # the number of generation we want to create
            for i in range(number_of_iter):
                glob = 0
                loc = 0
                switch_prob = f.prob_switch(initial_prob, number_of_iter, i)
                # population for each generation
                for j in range(number_of_population):
                    rand = f.random.uniform(0, 1)
                    if rand < switch_prob:
                        glob += 1
                    else:
                        loc += 1
                #print("Global ", glob, " Local ", loc)

        break


def start():
    starting_time = time.perf_counter()
    search("testFileNME.txt", 0.5, 50, 10, 0.1, 2, 1, 10)
    ending_time = time.perf_counter()
    print("\nCOMPUTATION TIME:", ending_time - starting_time, "seconds")


start()
