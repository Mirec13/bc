import csv


def load_file(file_name):
    # append the path to the file
    file_name = "datasets" + "/" + file_name

    # load the csv file in 2D list
    with open(file_name, newline='') as csv_file:
        data = list(csv.reader(csv_file))

    # delete the name of the columns from list
    del data[0]

    return data
