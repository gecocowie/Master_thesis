#!/usr/bin/env python
# coding: utf-8

# In[80]:



"""
Read all Mseries data into one big dataframe
Find the path to files in subfolders using os.listdir
"""
def read_M_series_subfolders(path):
    # create a list of file and sub directories 
    # names in the given directory 
    subfolders = os.listdir(path)
    file_list = []
    # Iterate over all the entries
    for names in subfolders:
        # Create full path
        fullPath = os.path.join(path, names)
        # If entry is a directory then get the list of files in this directory 
        temp_path = fullPath#for glob
        if os.path.isdir(fullPath):
            file_list = file_list + read_subfolders(fullPath)   
        else:
            txt_files = glob.glob(path + "/M*.txt")#only read txt files
            if (txt_files != []):
                file_list.append(txt_files)
    return file_list

path = "/Users/georgecowie/Documents/Master/Masteroppgave/data/Drifters/Supraglacial channel"#common path
filenames = np.array(read_subfolders(path),dtype=object)
filenames = np.unique(filenames)#remove duplicate files

k = 0
for i in range(len(filenames)):
    for j in range(len(filenames[i])):
        print(filenames[i][j])
print(k)

