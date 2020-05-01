import functions as f
import numpy as np
import time
import random


def search(prefix_file_name, initial_prob, number_of_iter, number_of_population,  num_to_ban,
           best_n, p_value, number_of_epi, repeat, min_value_for_df, index_beta):
    # define the file counter
    file_counter = 0

    # declare structs and other variables
    snp_data = f.declare_snp_data_struct()
    comb = pow(3, number_of_epi)
    flowers = []
    prev_flowers = []
    ban_iter = 0
    tabu = []

    # these three variables are for storing the data evaluation results
    TP = 0
    FP = 0
    FN = 0
    while True:
        file_name = prefix_file_name + str(file_counter) + ".antesnp100.txt"
        print(file_name)
        # Load the SNP data
        snp_data.data = f.load_file(file_name)

        # if we do not have any data files left end the loop
        if not snp_data.data:
            break

        for repeating in range(repeat):
            # declare a non-dominated list
            non_dominated = []
            tabu = []
            ban_iter = 0

            # declare struct for flowers
            f.declare_flowers_struct(flowers, number_of_population, number_of_epi, comb)
            f.declare_flowers_struct(prev_flowers, number_of_population, number_of_epi, comb)

            # Load the SNP data
            snp_data.data = f.load_file(file_name)

            # initialize SNP data
            snp_data.snp_size = len(snp_data.data[0]) - 1
            snp_data.sample_size = len(snp_data.data)
            snp_data.state = []
            for row in snp_data.data:
                snp_data.state.append(row[snp_data.snp_size])
                del row[snp_data.snp_size]

            # code for shuffle the columns or just reverse it
            #for i in range(len(snp_data.data)):
                #random.seed(file_counter)
                #random.shuffle(snp_data.data[i])
                #snp_data.data[i].reverse()

            # initialize the first population
            vector = f.init_first_population(number_of_epi, snp_data.snp_size, number_of_population)

            # print the init population
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
            best = f.best_solution(non_dominated_tmp, non_dominated, min_value_for_df, comb, number_of_epi, tabu)
            ban_iter += 1
            GS = 0
            LS = 0
            # the number of generation we want to create
            for i in range(number_of_iter):
                # swap
                tmp = prev_flowers
                prev_flowers = flowers
                flowers = tmp
                # get the value for switching between global and local search
                switch_prob = f.prob_switch(initial_prob, number_of_iter, i)

                # population for each generation
                for j in range(number_of_population):
                    rand = f.random.uniform(0, 1)
                    if rand < switch_prob:
                        GS += 1
                        flowers[j].loci = f.global_search(index_beta, best, prev_flowers[j].loci, snp_data.snp_size)
                    else:
                        # get two random flowers from the previous population
                        LS += 1
                        prev_random_flower_one = int(round(np.random.uniform(0, number_of_population - 1)))
                        prev_random_flower_two = int(round(np.random.uniform(0, number_of_population - 1)))
                        flowers[j].loci = f.local_search(prev_flowers[j].loci, prev_flowers[prev_random_flower_one].loci,
                                                         prev_flowers[prev_random_flower_two].loci, snp_data.snp_size)

                    # check whether we have repeating snp in loci and if we do then change it
                    f.check_the_epistasis(flowers[j], number_of_epi, snp_data.snp_size)

                    # compute the fitting score via objective functions
                    f.compute_o_and_e_values(snp_data.data, snp_data.sample_size, flowers[j].loci,
                                             flowers[j].observed_value, flowers[j].expected_value, snp_data.state,
                                             number_of_epi)
                    flowers[j].objective_function_score[0] = f.k2_score(flowers[j].observed_value, comb)
                    flowers[j].objective_function_score[1] = f.gi_score(flowers[j].observed_value, comb,
                                                                        snp_data.sample_size)

                # get the first non-dominated solution and the best loci
                non_dominated_tmp = f.pareto_optimization(flowers, number_of_population)
                prev_best = best
                best = f.best_solution(non_dominated_tmp, non_dominated, min_value_for_df, comb, number_of_epi, tabu)

                # if the best solution was repeated then increment counter if not or the counter reached its limit
                # then reset the counter and add the epistasis in tabu table
                if (set(prev_best) == set(best)) and (ban_iter != (num_to_ban - 1)):
                    ban_iter += 1
                elif ban_iter == (num_to_ban - 1):
                    tabu.append(best)
                    ban_iter = 0

            # this entire section is for printing the results
            non_dominated_accepted = []

            for i in range(len(non_dominated)):
                print(non_dominated[i].g_dist, non_dominated[i].loci)

            found_snp = f.find_most_frequent_snp(non_dominated, number_of_iter, snp_data.snp_size)

            if found_snp != -1:
                non_dominated = []
                print("idem")
                f.get_comb_of_snp(non_dominated, found_snp, number_of_epi, snp_data.snp_size, snp_data.sample_size, comb, snp_data.data, snp_data.state, min_value_for_df)
            else:
                non_dominated = f.get_n_best(non_dominated, best_n)
                f.get_all_non_dominated_combinations(non_dominated, min_value_for_df, number_of_epi, comb,
                                                     snp_data.sample_size, snp_data.state, snp_data.data)

            p_value_final = p_value / f.comb_without_repetition(snp_data.snp_size, number_of_epi)

            for i in range(len(non_dominated)):
                print(non_dominated[i].g_dist, non_dominated[i].loci)

            for i in non_dominated:
                if p_value_final > i.g_dist:
                    f.add_to_non_dominated_accepted(i, non_dominated_accepted, number_of_epi)

            print("P_VALUE: ", p_value_final)

            print("-------------------------------------")
            for i in non_dominated_accepted:
                print(i.g_dist, i.loci)

            print("-------------------------------------")

            non_dominated_accepted = f.get_all_unique_epistasis(non_dominated_accepted, snp_data.snp_size)

            found = False
            for i in non_dominated_accepted:
                if set(i.loci) == set([0, 1]):
                    found = True
                else:
                    FP += 1
                print(i.g_dist, i.loci)

            if found:
                TP += 1
            else:
                FN += 1
            print("GS:", GS, "LS:", LS)
            print("\ncontrol:", "TP", TP, "FP", FP, "FN", FN, "\n")

        file_counter += 1

    print("\nfinal results:", "number of files tested:", file_counter, "each file was repeated", repeat, "times", "TP", TP, "FP", FP, "FN", FN, "\n")


def start():
    starting_time = time.perf_counter()
    search(prefix_file_name="ME76/76.1600.", initial_prob=0.5, number_of_iter=40, number_of_population=40,
           num_to_ban=5, best_n=30, p_value=0.05, number_of_epi=2, repeat=1, min_value_for_df=10, index_beta=1.5)
    ending_time = time.perf_counter()
    print("\nCOMPUTATION TIME:", ending_time - starting_time, "seconds")


start()
