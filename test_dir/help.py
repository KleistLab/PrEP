#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys


def process(file):
    df = pd.read_csv(file, index_col=0)
    df = df.apply(lambda x: np.round(x, 8), axis=1)
    df.to_csv(file)


if __name__ == '__main__':
    process(sys.argv[1])
