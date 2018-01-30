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
    parser.add_argument("-i", "--input", help="Input binary matrix file")
    parser.add_argument("-o", "--output", help="Output ASCII matrix file")
    args = parser.parse_args()

    bin_file = open(args.input, "rb")
    num_rows = struct.unpack('i', bin_file.read(4))[0]
    num_columns = struct.unpack('i', bin_file.read(4))[0]

    ascii_file = open(args.output, "w")
    ascii_file.write("%i\t%i\n" % (num_rows, num_columns))
    for i in range(0, num_rows):
        for j in range(0, num_columns):
            num = struct.unpack("d", bin_file.read(8))[0]
            ascii_file.write("%f\t" % num)
            
        ascii_file.write("\n")

    bin_file.close()
    ascii_file.close()

if __name__ == "__main__":
    main()
