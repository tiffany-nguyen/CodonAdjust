print("importing libraries...")

import pandas as pd
import numpy as np
import csv
import sys
import logging
from logging import getLogger, StreamHandler, FileHandler, Formatter
from mkdir_p import mkdir_p # library for command mkdir -p (create folder if not existed)

INDIR = sys.argv[1]
LEN = int(sys.argv[2])

OUTDIR  = INDIR + "/BEST/"
OUT_AA = OUTDIR + "aa_opt.all.best.csv"
OUT_NT = OUTDIR + "nt_opt.all.best.csv"
OUT_MSE = OUTDIR + "MSE_opt.all.best.csv"
OUT_CODE = OUTDIR + "code_opt.all.best.csv"

mkdir_p(OUTDIR)

IUPAC = ["A", "C", "G", "T", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"]

# ----- initialize best results ------------
best_mse = pd.Series(1.0, index = range(LEN))
best_code = pd.Series("NNN", index = range(LEN))
best_aa = pd.DataFrame(np.zeros(shape = (20, LEN)))*1.0
best_nt = pd.DataFrame(np.zeros(shape = (12, LEN)))*1.0
# ----- initialize best results ------------

# ----- initialize logger ------------
logger = getLogger("IUPAC")
logger.setLevel(logging.DEBUG)
handler_format = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

stream_handler = StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(handler_format)

file_handler = FileHandler(OUTDIR + 'iupac.log', 'a')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(handler_format)

logger.addHandler(stream_handler)
logger.addHandler(file_handler)
# ----- initialize logger ------------

#read csv
def readcsv(infile):
    df = pd.read_csv(infile, index_col = 0, header = None)
    return df

def find_best_matches():
    for code1 in IUPAC:
        for code2 in IUPAC:
            for code3 in IUPAC:
                param = code1 + code2 + code3
                curr_mse = INDIR + "/" + param + "/MSE_opt.all.csv"
                df_mse = readcsv(curr_mse)

                for i in range(LEN):
                    if best_mse[i] > df_mse.iloc[0, i]:
                        logger.info("Found better IUPAC code for position %s. Current code: %s will be updated to new code %s.", str(i+1), best_code[i], param)
                        curr_aa = INDIR + "/" + param + "/aa_opt.all.csv"
                        curr_nt = INDIR + "/" + param + "/nt_opt.all.csv"
                        df_aa = readcsv(curr_aa)
                        df_nt = readcsv(curr_nt)
                        best_aa.index = df_aa.index
                        best_nt.index = df_nt.index

                        best_mse[i] = df_mse.iloc[0, i]
                        best_code[i] = param
                        best_aa.iloc[:, i] = df_aa.iloc[:, i]
                        best_nt.iloc[:, i] = df_nt.iloc[:, i]

    best_aa.to_csv(OUT_AA, header = False)
    best_nt.to_csv(OUT_NT, header = False)

    df_mse = pd.DataFrame(best_mse).T
    df_mse.index = ["best_mse"]
    df_mse.to_csv(OUT_MSE, header = False)

    df_code = pd.DataFrame(best_code).T
    df_code.index = ["best_code"]
    df_code.to_csv(OUT_CODE, header = False)


def main():
    try:
        find_best_matches()
        print("done")
    except:
        logger.error("Exception occurred.")

if __name__ == '__main__':
    main()
