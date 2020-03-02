import csv


def load_file(file_name):
    try:
        # append the path to the file
        file_name = "datasets" + "/" + file_name

        # load the csv file in 2D list
        with open(file_name, newline='') as csv_file:
            data = list(csv.reader(csv_file))

        # delete the name of the columns from list
        del data[0]

        # convert data's components to int
        y = len(data)
        x = len(data[0])
        for i in range(y):
            for j in range(x):
                data[i][j] = int(data[i][j])

        return data
    except FileNotFoundError:
        data = False
        return data


def init_first_population():
    pass


def compute_o_and_e_values(data, sample_size, loci, observed_value, expected_value, state):
    comb_number = 0

    # compute the observed value for every combination's state
    for i in range(3):
        for j in range(3):
            for k in range(sample_size):
                if (data[k][loci[0]] == i) and (data[k][loci[1]] == j):
                    if state[k] == 0:
                        observed_value[comb_number][0] += 1
                    else:
                        observed_value[comb_number][1] += 1
            comb_number += 1

    # compute the expected value
    for i in range(comb_number):
        expected_value[i] = ((observed_value[i][0] + observed_value[i][1]) / 2)


def k2_score():
    pass


def gi_score():
    pass


def g_test():
    pass
