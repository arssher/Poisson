#! /usr/bin/env python2.7

from os import listdir, path
import json
import argparse
from subprocess import check_call

import numpy as np

def collect_data():
    data = []

    for f in listdir(res_dir):
        f_path = path.join(res_dir, f)
        if path.isfile(f_path):
            print "Importing data from %s" % f_path
            with open(f_path, "r") as df:
                data.extend(json.load(df))

    print data;

def correct_solution(x, y):
    return np.sin(x * y) + 1.0

def draw_correct():
    n = 100
    fname = "correct.data"

    with open(fname, "w") as gpltf:
        dots = np.linspace(0, 2, num=n)
        for x in dots:
            for y in dots:
                gpltf.write(str(x) + " " + str(y) + " " + str(correct_solution(x, y)) + "\n")
            gpltf.write("\n")

    check_call(["gnuplot", "-e",
                "dataf='%s'; outf='correct_sol.pdf'; tit='Correct solution'; zl='u(x, y)'" % fname,
                "draw.gpt"])



if __name__ == '__main__':
    draw_correct()
