import functions as f
from collections import namedtuple


def search(file_name, prob_switch, number_of_iter, number_of_population, number_of_repetition, p_value, num_of_epi):
    # define structs
    snp_data = namedtuple("snp_data", "sample_size snp_size data state")
    flower = [namedtuple("flower", "loci observed_value expected_value objective_function")] * number_of_population
    comb = pow(3, num_of_epi)
    for i in flower:
        i.observed_value = [[0 for x in range(2)] for y in range(comb)]
        i.expected_value = [0] * comb
        i.loci = [0] * num_of_epi
        i.objective_function = [0.0] * 2

    while True:

        # initialize SNP data
        snp_data.data = f.load_file(file_name)

        # if we do not have any data end the loop
        if not snp_data.data:
            break

        snp_data.snp_size = len(snp_data.data[0]) - 1
        snp_data.sample_size = len(snp_data.data)
        snp_data.state = []
        for row in snp_data.data:
            snp_data.state.append(row[snp_data.snp_size])
            del row[snp_data.snp_size]

        # init manually one flower for test purposes
        flower[0].loci[0] = 0
        flower[0].loci[1] = 1
        f.compute_o_and_e_values(snp_data.data, snp_data.sample_size, flower[0].loci,
                                 flower[0].observed_value, flower[0].expected_value, snp_data.state)
        print(flower[0].observed_value, flower[0].expected_value)

        break


def start():
    search("testFile.txt", 1, 2, 1, 4, 5, 2)


start()