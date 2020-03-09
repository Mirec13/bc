import functions as f
from scipy.stats import chi2
import time
from collections import namedtuple


def search(file_name, initial_prob, number_of_iter, number_of_population, p_value, number_of_epi, repeat, min_value_for_df):
    # define structs
    snp_data = namedtuple("snp_data", "sample_size snp_size data state")
    flower = [namedtuple("flower", "loci observed_value expected_value objective_function")] * number_of_population
    # number of possible combination for n correlations (epistasis)
    comb = pow(3, number_of_epi)
    for i in flower:
        i.observed_value = [[0 for x in range(2)] for y in range(comb)]
        i.expected_value = [0.0] * comb
        i.loci = [0] * number_of_epi
        i.objective_function = [0.0] * 2

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
                print("Global ", glob, " Local ", loc)

        break


def start():
    starting_time = time.perf_counter()
    search("testFileNME.txt", 0.5, 50, 50, 0.1, 2, 1, 10)
    ending_time = time.perf_counter()
    print("\nCOMPUTATION TIME:", ending_time - starting_time, "seconds")


start()

# init manually one flower for test purposes
# flower[0].loci[0] = 0
# flower[0].loci[1] = 1
# f.compute_o_and_e_values(snp_data.data, snp_data.sample_size, flower[0].loci,
#                          flower[0].observed_value, flower[0].expected_value, snp_data.state, number_of_epi)
# print(f.k2_score(flower[0].observed_value, comb))
# print(f.gi_score(flower[0].observed_value, comb, snp_data.sample_size))
# g = f.g_test(flower[0].observed_value, flower[0].expected_value, comb)
# print("g value", g)
# print("g-distribution: ", chi2.sf(g, 9))
# print("O and E values: ", flower[0].observed_value, flower[0].expected_value)
# print("P-value: ", 0.1/f.comb_without_repetition(100, number_of_epi))
# if chi2.sf(g, 9) < (0.1/f.comb_without_repetition(100, number_of_epi)):
#     print("yes")
# else:
#     print("no")
#  print(f.prob_switch(0.6, 50, 50))