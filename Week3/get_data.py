import csv
import os
import shutil

import numpy as np
import pandas as pd

if __name__ == '__main__':
    for directory in ['470ohm', '470ohm-b', '470ohm-b']:
        for file in os.listdir(directory):
            if os.path.isdir(directory + '/' + file):
                continue
            df = pd.read_csv(directory + '/' + file).values
            time_lists = []
            pos_lists = []
            scales = []
            index = 3
            while index < len(df[0]):
                time_lists.append(df[:, index])
                pos_lists.append(df[:, index + 1])
                scales.append(df[7][index - 2])
                index += 6
            if os.path.isdir(directory + '/' + file[:-4]):
                shutil.rmtree(directory + '/' + file[:-4])
            os.mkdir(directory + '/' + file[:-4])
            for i in range(len(time_lists)):
                with open(f'{directory}/{file[:-4]}/meas{i + 1}.csv', 'w', newline='') as newfile:
                    wr = csv.writer(newfile)
                    wr.writerow(['time', 'position', 'scale'])
                    to_write = [list([time_lists[i][j], pos_lists[i][j], float(scales[i])]) for j in range(len(time_lists[i]))]
                    wr.writerows(to_write)

