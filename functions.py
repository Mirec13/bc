import csv
import math
import numpy as np
import random
import itertools as c


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


def init_first_population(number_of_dimension, number_of_loci, population):
    random.seed(int(random.uniform(0, 100)))
    list_of_loci = [i for i in range(number_of_loci)]
    random.shuffle(list_of_loci)
    subset = []
    i = 0

    while i < population:
        epistasis = random.sample(list_of_loci, number_of_dimension)
        match = False
        for vector in subset:
            if set(vector) == set(epistasis):
                match = True
                break
        if match:
            continue
        else:
            subset.append(epistasis)
            i += 1

    return subset


def compute_o_and_e_values(data, sample_size, loci, observed_value, expected_value, state, number_of_epi):
    comb_number = 0
    match = 0

    # compute the observed value for every combination's state
    for i in range(sample_size):
        for j in c.product('012', repeat=number_of_epi):
            for k in range(number_of_epi):
                if data[i][loci[k]] == int(j[k]):
                    match += 1
            if match == number_of_epi:
                if state[i] == 0:
                    observed_value[comb_number][0] += 1
                else:
                    observed_value[comb_number][1] += 1
            match = 0
            comb_number += 1
        comb_number = 0

    comb_number = pow(3, number_of_epi)

    # compute the expected value
    for i in range(comb_number):
        expected_value[i] = ((observed_value[i][0] + observed_value[i][1]) / 2)


def k2_score(observed_value, comb):
    final_score = 0.0
    sub_score = [0.0] * 2

    for i in range(comb):
        num_o_value = observed_value[i][0] + observed_value[i][1]
        for j in range(num_o_value+1):
            final_score += math.log(j+1)
        for j in range(2):
            for k in range(observed_value[i][j]):
                sub_score[j] += math.log(k+1)
        final_score -= (sub_score[0] + sub_score[1])
        sub_score[0] = 0.0
        sub_score[1] = 0.0

    return final_score


def gi_score(observed_value, comb, sample_size):
    final_score = 0.0
    sub_score = [0.0] * 2

    for i in range(comb):
        num_o_value = observed_value[i][0] + observed_value[i][1]
        for j in range(2):
            if num_o_value > 0:
                sub_score[j] = (observed_value[i][j] / num_o_value)
                sub_score[j] = pow(sub_score[j], 2)
        final_score += (num_o_value / sample_size) * (1 - (sub_score[0] + sub_score[1]))
        sub_score[0] = 0.0
        sub_score[1] = 0.0

    return final_score


def pareto_optimization(flowers, population, non_dominated):
    for i in range(population):
        dominated = False
        for j in range(population):
            if (((flowers[j].objective_function_score[0] < flowers[i].objective_function_score[0]) and
                 (flowers[j].objective_function_score[1] < flowers[i].objective_function_score[1])) or
                ((flowers[j].objective_function_score[0] == flowers[i].objective_function_score[0]) and
                 (flowers[j].objective_function_score[1] < flowers[i].objective_function_score[1])) or
                ((flowers[j].objective_function_score[0] < flowers[i].objective_function_score[0]) and
                 (flowers[j].objective_function_score[1] == flowers[i].objective_function_score[1]))):
                dominated = True
                break
            if (((flowers[j].objective_function_score[0] == flowers[i].objective_function_score[0]) and
                (flowers[j].objective_function_score[1] == flowers[i].objective_function_score[1])) and
               (i < j)):
                dominated = True
                break
        if not dominated:
            non_dominated.append(flowers[i].objective_function_score)


def g_test(observed_value, expected_value, comb):
    final_score = 0.0
    sub_score = [0.0] * 2

    for i in range(comb):
        for j in range(2):
            if expected_value[i] > 0.0:
                prob = observed_value[i][j] / expected_value[i]
            else:
                prob = 0
            if prob != 0:
                sub_score[j] = (observed_value[i][j]) * (math.log(prob))
        final_score += (sub_score[0] + sub_score[1])
        sub_score[0] = 0.0
        sub_score[1] = 0.0

    return 2 * final_score


def subtract_df(observed_values, min_value, comb):
    df = comb
    for o in observed_values:
        if (o[0] < min_value) and (o[1] < min_value):
            df -= 1

    return df


def comb_without_repetition(number, number_rep):
    numerator = number - (number - number_rep)
    denominator = number_rep
    result_numerator = 1
    result_denominator = 1

    for i in range(numerator):
        result_numerator *= number
        number -= 1

    for i in range(denominator):
        result_denominator *= number_rep
        number_rep -= 1

    final_result = result_numerator / result_denominator

    return final_result


def prob_switch(initial_prob, number_of_iteration, actual_iteration):
    prob = initial_prob + (0.1 * ((number_of_iteration - actual_iteration) / number_of_iteration))
    return prob


def local_search():
    pass


def levy_flight(index_beta):
    exponent = 1 / index_beta
    numerator = math.gamma(1 + index_beta) * math.sin((math.pi * index_beta) / 2)
    denominator = math.gamma((1 + index_beta) / 2) * index_beta * math.pow(2, (index_beta - 1) / 2)
    base = math.pow(numerator / denominator, exponent)

    u = np.random.sample(0, base)
    v = np.random.sample(0, 1)

    random_step = u / math.pow(v, 1 / index_beta)

    return random_step


def global_search(index_beta, prev_best_flower, prev_flower):
    pass






