import os

if __name__ == '__main__':
    for directory in ['470ohm', '470ohm-b', '470ohm-c']:
        lst = []
        for file in os.listdir(directory):
            if os.path.isdir(directory + '/' + file):
                lst.append(file)
        print(lst)



