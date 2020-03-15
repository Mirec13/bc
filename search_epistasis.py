import functions as f
import numpy as np
import time

# TODO THINK THROUGH THE DUPLICATES IN NON_DOMINATED SOLUTION
def search(file_name, initial_prob, number_of_iter, number_of_population, p_value, number_of_epi, repeat,
           min_value_for_df, index_beta):
    # declare structs and other variables
    comb = pow(3, number_of_epi)
    flowers = []
    prev_flowers = []
    snp_data = f.declare_snp_data_struct()
    f.declare_flowers_struct(flowers, number_of_population, number_of_epi, comb)
    f.declare_flowers_struct(prev_flowers, number_of_population, number_of_epi, comb)

    while True:
        # Load the SNP data
        snp_data.data = f.load_file(file_name)

        # if we do not have any data files left end the loop
        if not snp_data.data:
            break

        for repeating in range(repeat):
            non_dominated = []
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

            # get the first non-dominated solution and the best loci
            non_dominated_tmp = f.pareto_optimization(flowers, number_of_population)
            prev_best = f.best_solution(non_dominated_tmp, non_dominated, min_value_for_df, comb, number_of_epi)

            # for i in flowers:
            #     print(i.objective_function_score)

            # print(non_dominated[0].g_dist, non_dominated[0].loci, prev_best)
            # print(f.global_search(1.5, flowers[1].loci, flowers[0].loci, snp_data.snp_size))
            # print(f.local_search(flowers[0].loci, flowers[1].loci, flowers[2].loci, snp_data.snp_size))

            # the number of generation we want to create
            for i in range(number_of_iter):
                glob = 0
                loc = 0
                # swap
                tmp = prev_flowers
                prev_flowers = flowers
                flowers = tmp
                switch_prob = f.prob_switch(initial_prob, number_of_iter, i)
                # population for each generation
                for j in range(number_of_population):
                    rand = f.random.uniform(0, 1)
                    if rand < switch_prob:
                        flowers[j].loci = f.global_search(index_beta, prev_best, prev_flowers[j].loci, snp_data.snp_size)
                        glob += 1
                    else:
                        prev_random_flower_one = int(round(np.random.uniform(0, number_of_population - 1)))
                        prev_random_flower_two = int(round(np.random.uniform(0, number_of_population - 1)))
                        flowers[j].loci = f.local_search(prev_flowers[j].loci, prev_flowers[prev_random_flower_one].loci,
                                                         prev_flowers[prev_random_flower_two].loci, snp_data.snp_size)
                        loc += 1
                    f.compute_o_and_e_values(snp_data.data, snp_data.sample_size, flowers[j].loci,
                                             flowers[j].observed_value, flowers[j].expected_value, snp_data.state,
                                             number_of_epi)
                    flowers[j].objective_function_score[0] = f.k2_score(flowers[j].observed_value, comb)
                    flowers[j].objective_function_score[1] = f.gi_score(flowers[j].observed_value, comb,
                                                                        snp_data.sample_size)

                # get the first non-dominated solution and the best loci
                non_dominated_tmp = f.pareto_optimization(flowers, number_of_population)
                prev_best = f.best_solution(non_dominated_tmp, non_dominated, min_value_for_df, comb, number_of_epi)
               # for k in non_dominated:
                 #   print(k.g_dist, k.loci)
                #print()
                #f.pareto_optimization(flowers, number_of_population, non_dominated_tmp)
                # print("Global ", glob, " Local ", loc)

            for i in non_dominated:
                print(i.g_dist, i.loci)

        break


def start():
    starting_time = time.perf_counter()
    search(file_name="testFileME.txt", initial_prob=0.5, number_of_iter=50, number_of_population=50, p_value=0.1,
           number_of_epi=2, repeat=1, min_value_for_df=10, index_beta=1.5)
    ending_time = time.perf_counter()
    print("\nCOMPUTATION TIME:", ending_time - starting_time, "seconds")


start()
