import os
import concurrent.futures
import hashlib
import re
import sys
import pathlib
import pandas as pd
from Bio.Seq import Seq
from collections import deque

# packages installed:
# pip install Biopython
# pip install pandas
# pip intall openpyxl


# command line arguments:
# 1: excel file containing sample id's and indices

# recursively find all files (excluding hidden files) from current directory
def scan_files(path):
    
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            yield from scan_files(entry.path)
        elif entry.name.startswith('.'):
            continue
        else:
            yield entry
            
# group files by sample. Each sample will have it's own list with relative path to files
def group_smp_files(paths):
    
    smp_file_ls = []
    tmp_smp_file_ls = []
    # length of run identifier + i7 + i5
    common_path_length = 24
    
    for entry in paths:
        
        # if the list currently contains a path
        if tmp_smp_file_ls:
            
            # calculate the common path size between entry and previous path in list
            common = os.path.commonpath([entry, tmp_smp_file_ls[-1]])
            # print("common: ", len(common))
            
            # return the current list, create a new empty list, add entry
            if len(common) < common_path_length:
                smp_file_ls.append(tmp_smp_file_ls)
                tmp_smp_file_ls = []
                tmp_smp_file_ls.append(entry)
            else:
                tmp_smp_file_ls.append(entry)
                
        # the list does not currently contain a path, add entry        
        else:
            tmp_smp_file_ls.append(entry)
    
    # return list of paths for final sample        
    if smp_file_ls:
        smp_file_ls.append(tmp_smp_file_ls)
    
    return smp_file_ls
        
# returns i7/i5 indices from file name
# retuns as list with i7 in pos 0 and i5 in pos 1
def get_indices(file_paths):

    # for each of the files extract the i7 and i5 indices, save as list of lists
    indices = []
    
    for f in file_paths:
        
        file = os.path.basename(f[0])
        sample = []
        
        try:
            i7 = re.search('_[ATCGU]{8}', file).group(0)
        except AttributeError:
            print("i7 index not found in file name, exiting")
            quit()
        
        sample.append(i7[1:])      
        
        try:
            i5 = re.search('-[ATCGU]{8}', file).group(0)
        except AttributeError:
            print("i5 index not found in file name, exiting")
            quit()
        
        sample.append(i5[1:])
        indices.append(sample)
    
    
    if indices[0] != indices[1]:
        print("sample indices do not match, exiting")
        quit()
        
    return indices[0]       
    
    
# rname fastq files according to corresponding sample_id determined from i7/i5
def rename_fastq_files(file_paths, sample_ids_file):
    
    # entries will be a tuple -> (file path, file name)
    rename_files = []
    
    # list of updated files to return
    updated_file_paths = deque() 
    
    for file in file_paths:
        if "fastq" in file:
            path = file
            file_name = os.path.basename(file)
            rename_files.append((path, file_name))
        else:
            updated_file_paths.append(os.path.abspath(file))
            

    # returns a list
    # indices[0] -> i7, indices[1] -> i5
    indices = get_indices(rename_files)
    
    # resolve indices and their reverse complement
    i7 = Seq(indices[0]) 
    i5 = Seq(indices[1])
    i5_rc = i5.reverse_complement()
    i7_rc = i7.reverse_complement()
    i5 = "-" + i5
    i5_rc = "-" + i5_rc
    
    # read in excel file with sample info as df
    samples_df = pd.read_excel(sample_ids_file)
    samples_df = samples_df[['sample', 'i7', 'i5']]

    # find corresponding sample id given i5/i7
    result_df = samples_df.loc[((samples_df['i7'] == i7) | (samples_df['i7'] == i7_rc)) & ((samples_df['i5'] == i5) | (samples_df['i5'] == i5_rc))]
    
    # check that a single corresponding sample id was found
    if result_df.empty or result_df.shape[0] > 1:
        print("issue resolving sample id for given indices, i7: ", i7, " i5: ", i5)
        quit()
    
    # rename files, update file and file path
    for i, file in enumerate(rename_files):

        # determine new file name and path
        a = re.sub(r'\w*[ATCGU]{8}.[ATCGU]{8}', result_df.iloc[0]['sample'], file[1])
        start = re.search(result_df.iloc[0]['sample'] + r'_\d', a).group(0)
        end = re.search(r'.fastq\w*', file[1]).group(0)
        new_name = start + end
        new_path = re.sub(file[1], new_name, file[0])

        # rename file        
        os.rename(file[0], new_path)
        # update list 
        rename_files[i] = (new_path, new_name)
    
    updated_file_paths.extendleft((os.path.abspath(rename_files[0][0]), os.path.abspath(rename_files[1][0]), result_df.iloc[0]['sample']))
    return list(updated_file_paths)

# create text file containing a checksum for each file present        
def md5_checksum(files):

    # read file in 64kb chunks  
    BUF_SIZE = 64 * 2^10
    dir = os.path.commonpath(files)
    name = dir + "/md5_checksum.txt"
    hashes = []
    
    # create hashes for each file in increments 
    for file in files:
        
        # create hashes using md5 hash function
        md5 = hashlib.md5()
        
        with open(file, "rb") as f:
            while True:
                data = f.read(BUF_SIZE)
                if not data:
                    break
                md5.update(data)
        
        # add complete hash for file
        hashes.append(md5.hexdigest()) 
    
    # create a write to file to store checksums
    with open(name, "w") as checksum_file:
        for index, file in enumerate(files):
            checksum_file.write(hashes[index] + "\t" + os.path.basename(file) + "\n")
    
    return name
    
def main():
    
    # recursively find all files from current directory using a generator function
    # must use relative path until after "group_smp_files" fn call
    file_ls = scan_files("./")
    
    # create a list of all files
    path_ls = []
    for entry in file_ls:
        path_ls.append(entry.path)

    # group files together into seperate lists for each sample
    file_ls_grouped = group_smp_files(path_ls)
    print("\n\n\nfile_ls_grouped before checksum and renaming: \n", file_ls_grouped)

    # generate check sum text file, rename fastq files and create a list of seperate deques for each sample
    # first entry of deque contains sample id, other entries are the relative file paths 
    for index, entry in enumerate(file_ls_grouped):
        file_ls_grouped[index] = rename_fastq_files(entry, sys.argv[1])
        checksum_file_path = md5_checksum(file_ls_grouped[index][1:])
        file_ls_grouped[index].append(os.path.abspath(checksum_file_path))
        
    print("\n\n\nfile_ls_grouped after checksum and renaming: \n", file_ls_grouped)
    
if __name__ == "__main__":
    main()

        
    