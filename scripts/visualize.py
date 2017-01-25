#! /usr/bin/env python2.7

from os import listdir, path, makedirs
import csv
import argparse
from subprocess import check_call

import numpy as np

# each file in res_dir must be comma-separated csv file with x, y, u columns
# merge them and return [[x, y, u], [x, y, u], ...]
def collect_data(res_dir):
    data = []

    for f in listdir(res_dir):
        f_path = path.join(res_dir, f)
        if path.isfile(f_path):
            print "Importing data from %s" % f_path
            with open(f_path, "rb") as csvf:
                freader = csv.reader(csvf, delimiter=',')
                for row in freader:
                    data.append([float(row[0]), float(row[1]), float(row[2])])
    return data

# run scripts/draw.gpt with given args and put the result to drawings/'outf'
def run_gnuplot(dataf, outf, tit, zl):
    if not path.exists("drawings"):
        makedirs("drawings")

    check_call(["gnuplot", "-e",
                "dataf='%s'; outf='drawings/%s'; tit='%s'; zl='%s'" %
                (dataf, outf, tit, zl),
                "scripts/draw.gpt"])

# Pulls *.csv files with my solution from res/ dir and draws it.
def draw_mine():
    fname = "/tmp/mine.data"
    data = collect_data("res")
    # we sort the data by x and y and insert newline when x changes so that gnuplot
    # can distinguish lines
    data_sorted = sorted(data, key=lambda row: (row[0], row[1]))
    prev_x = data_sorted[0][0]
    with open(fname, "w") as gpltf:
        for row in data_sorted:
            if prev_x != row[0]:
                prev_x = row[0]
                gpltf.write("\n")
            gpltf.write(str(row[0]) + " " + str(row[1]) + " " + str(row[2]) + "\n")

    run_gnuplot(fname, "mine_sol.pdf", "Parallel solution", 'u(x, y)')


def correct_solution(x, y):
    return np.sin(x * y) + 1.0

def draw_correct():
    n = 100
    fname = "/tmp/correct.data"

    with open(fname, "w") as gpltf:
        dots = np.linspace(0, 2, num=n)
        for x in dots:
            for y in dots:
                gpltf.write(str(x) + " " + str(y) + " " + str(correct_solution(x, y)) + "\n")
            # insert newline when x changes so that gnuplot can distinguish lines
            gpltf.write("\n")

    run_gnuplot(fname, "correct_sol.pdf", "Exact solution", "u(x, y)")

# draw absolute error
def draw_error():
    fname = "/tmp/error.data"

    data = collect_data("res")
    # we sort the data by x and y and insert newline when x changes so that gnuplot
    # can distinguish lines
    data_sorted = sorted(data, key=lambda row: (row[0], row[1]))
    prev_x = data_sorted[0][0]
    with open(fname, "w") as gpltf:
        for row in data_sorted:
            if prev_x != row[0]:
                prev_x = row[0]
                gpltf.write("\n")
            x, y = row[0], row[1]
            err = abs(correct_solution(x, y) - row[2])
            gpltf.write(str(x) + " " + str(y) + " " + str(err) + "\n")

    run_gnuplot(fname, "abs_error.pdf", "Absolute error", '')

# calculate max error without any drawings. We need this because gnuplot couldn't
# swallow 2000x2000 dots on my machine.
def calc_max_error():
    data = collect_data("res")
    max_error = 0.0
    for row in data:
        err = abs(correct_solution(row[0], row[1]) - row[2])
        if err > max_error:
            max_error = err

    print "Max error is %s" % max_error

if __name__ == '__main__':
    draw_mine()
    draw_correct()
    draw_error()
    calc_max_error()
