import functions as f
from collections import namedtuple


def search(file_name, prob_switch, number_of_iter, number_of_population, number_of_repetition, p_value, num_of_epi):
    snp_data = namedtuple("snp_data", "sample_size snp_size data state")

    # initialize SNP data
    snp_data.data = f.load_file(file_name)
    snp_data.snp_size = len(snp_data.data[0]) - 1
    snp_data.sample_size = len(snp_data.data)
    snp_data.state = []
    for row in snp_data.data:
        snp_data.state.append(row[snp_data.snp_size])
        del row[snp_data.snp_size]
        print(row)

    print(snp_data.sample_size, snp_data.snp_size, snp_data.state)


def start():
    search("testFile.txt", 1, 2, 3, 4, 5, 6);


start();