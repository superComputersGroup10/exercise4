#!/usr/bin/python3

import argparse
import struct


def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


def main():
    parser = argparse.ArgumentParser(description="Converts the input matrix from ASCII to binary")
    parser.add_argument("-i", "--input", help="Input ascii matrix file")
    parser.add_argument("-o", "--output", help="Output binary matrix file")
    args = parser.parse_args()

    with open(args.output, "wb") as bin_file:
        with open(args.input, "r") as ascii_file:
            first_line = ascii_file.readline()
            dimensions = [int(n) for n in first_line.replace('\n', '').split('\t') if n.isdecimal()]
            for n in dimensions:
                bin_file.write(struct.pack('i', n))
            for line in ascii_file:
                numbers = [float(n) for n in line.replace("\n", "").split("\t") if isfloat(n)]
                for n in numbers:
                    bin_file.write(struct.pack('d', n))




if __name__ == "__main__":
    main()
