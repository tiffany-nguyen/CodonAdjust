#-------------------------------------------------------------------------------
# Name:        mk_iupac_nt_input
# Purpose:     make nucleotide input of all 3375 IUPAC code
#
# Author:      Tiffany
#
# Created:     10/11/2019
# Copyright:   (c) Tiffany 2019
#-------------------------------------------------------------------------------
import pandas as pd
import numpy as np
import csv
import sys
from mkdir_p import mkdir_p # library for command mkdir -p (create folder if not existed)

iupac_code = {"A","C","G","T","W","S","M","K","R","Y","B","D","H","V","N"}
length = sys.argv[1]

OUTPUT = "./IUPAC_input/"
mkdir_p(OUTPUT)

#read csv
def readcsv(infile):
    print("reading csv...")
    df = pd.read_csv(infile, index_col=0)
    return df

def gen_iupac_nt_input():
    for code1 in iupac_code:
        for code2 in iupac_code:
            for code3 in iupac_code:
                output_df = pd.DataFrame()
                iucode = code1 + code2 + code3
                outfile = OUTPUT + "nt_" + iucode + ".csv"
                for i in range(1, (int(length) + 1)):
                    output_df = pd.concat([output_df, load_iupac_nt(iucode)], axis=1)
                output_df.to_csv(outfile, index = False, header = False)
                
def load_iupac_nt(iucode):
    filename = "iupac_code/nucleotide_freq_" + iucode + ".csv"
    df = readcsv(filename)
    # sort dataframe by index (A, C, G , T)
    df.sort_index(inplace=True)

    output_df = pd.DataFrame()
    output_arr = np.array([])
    for i in range(1, 4):
        col_name = "pos" + str(i)
        tmp_arr = df[col_name].values
        tmp_arr = np.append(output_arr, tmp_arr)
        tmp_s = pd.Series(tmp_arr)
        output_df = output_df.append(pd.DataFrame(tmp_s))

    return output_df

def main():
    print("Start reading csv...")
    gen_iupac_nt_input()
    print("done")

if __name__ == '__main__':
    main()

