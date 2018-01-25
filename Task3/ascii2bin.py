#!/usr/bin/python3

import argparse
import struct


def main():
    parser = argparse.ArgumentParser(description="Converts the input matrix from ASCII to binary")
    parser.add_argument("-i", "--input", help="Input ascii matrix file")
    parser.add_argument("-o", "--output", help="Output binary matrix file")
    args = parser.parse_args()

    with open(args.output, "wb") as output_file:
        output_file.write(struct.pack('i', 1000))




if __name__ == "__main__":
    main()
