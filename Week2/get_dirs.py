import os

if __name__ == '__main__':
    for directory in ['PartA', 'PartB', 'PartB-LowPass']:
        lst = []
        for file in os.listdir(directory):
            if os.path.isdir(directory + '/' + file):
                lst.append(file)
        print(lst)



