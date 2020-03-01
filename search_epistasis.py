import functions as f
from collections import namedtuple

def search(file_name, prob_switch, number_of_iter, number_of_population, number_of_repetition, p_value, num_of_epi):
    snp_data = namedtuple("snp_data", "sample_size snp_size data");

    snp_data.data = f.load_file("testFile.txt");
    snp_data.snp_size = len(snp_data.data[0]);
    snp_data.sample_size = len(snp_data.data);

    print(snp_data.sample_size, snp_data.snp_size)


def start():
    search("testFile.txt", 1, 2, 3, 4, 5, 6);


start();