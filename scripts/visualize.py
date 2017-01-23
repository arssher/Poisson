#! /usr/bin/env python2.7

from os import listdir, path, makedirs
import csv
import argparse
from subprocess import check_call

import numpy as np

# each file in res_dir must be comma-separated csv file with x, y, u columns
def collect_data(res_dir):
    data = []

    for f in listdir(res_dir):
        f_path = path.join(res_dir, f)
        if path.isfile(f_path):
            print "Importing data from %s" % f_path
            with open(f_path, "rb") as csvf:
                freader = csv.reader(csvf, delimiter=',')
                for row in freader:
                    data.append(row)
    return data

def comparator(row1, row2):
    return row1[0] < row2[0]

def run_gnuplot(dataf, outf, tit, zl):
    if not path.exists("drawings"):
        makedirs("drawings")

    check_call(["gnuplot", "-e",
                "dataf='%s'; outf='drawings/%s'; tit='%s'; zl='%s'" %
                (dataf, outf, tit, zl),
                "scripts/draw.gpt"])

def draw_mine():
    fname = "mine.data"
    data = collect_data("res")
    # we sort the data by x and insert newline when x changes so that gnuplot
    # can distinguish lines
    data.sort(cmp=comparator)
    prev_x = 42.0
    with open(fname, "w") as gpltf:
        for row in data:
            if prev_x != row[0]:
                prev_x = row[0]
                gpltf.write("\n")
            gpltf.write(str(row[0] + " " + str(row[1]) + " " + str(row[2]) + "\n"))

    run_gnuplot(fname, "mine_sol.pdf", "Incorrect solution", 'u(x, y)')


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

    run_gnuplot(fname, "correct_sol.pdf", "Correct solution", "u(x, y)")


if __name__ == '__main__':
    draw_mine()
    draw_correct()
