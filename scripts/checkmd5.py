#%%
import hashlib
from collections import defaultdict
import os
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description ="Check MD5")
    parser.add_argument('-i', '--input_dir', action = 'store', metavar = 'FILE',
                        help = 'Directory with files to ged MD5')
    parser.add_argument('-m', '--manifest', action = 'store', metavar = 'FILE',
                        help = 'Text file with file name and MD5 of original files')
    return parser.parse_args()
#%%
def main():

    opts = parse_arguments()

    # manifest = '/mnt/data/MANIFEST'

    hash_dict = defaultdict()
    with open(opts.manifest, 'r') as in_f:
        for line in in_f:
            file_name, hash_str = line.rstrip(' \n').split(' ')
            hash_dict[file_name] = hash_str


    #list all files in directory
    onlyfiles = [os.path.join(opts.input_dir, f) for f in os.listdir(opts.input_dir) if os.path.isfile(os.path.join(opts.input_dir, f))]

    for file_path in onlyfiles:
        check_md5(file_path, hash_dict)
#%%

def check_md5(file_path, hash_dict):
    '''Generate file md5 and compare it to one in hash dict
    '''
    file_name = os.path.basename(file_path)

    hash_md5 = hashlib.md5()
    with open(file_path, 'rb') as f_to_check:
        for chunk in iter(lambda: f_to_check.read(4096), b""):
            hash_md5.update(chunk)
    md5_returned = hash_md5.hexdigest()
        
    original_md5 = hash_dict.get(file_name, "Hash not in manifest")

    if original_md5 == md5_returned:
        print(file_name, "md5 matches original")
    else:
        print(file_name, original_md5, md5_returned)

if __name__ == "__main__":
    main()