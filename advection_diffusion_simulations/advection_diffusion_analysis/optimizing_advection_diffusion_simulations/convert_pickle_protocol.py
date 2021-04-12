import pickle

filename2 = 'new_pkl_file.pkl'

with open (filename2, 'rb') as g:
    data_list = []
    done = False
    while not done:
        try:
            data2 = pickle.load(g)
            print (data2)
        except EOFError as err:
            done = True
            break
        data_list.append(data2)

for i, data in enumerate(data_list):
    print (i)
    for key in data[0]:
        print (key)
    print ()
