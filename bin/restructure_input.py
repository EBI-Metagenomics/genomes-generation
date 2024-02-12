#!/usr/bin/env python3
import argparse
import os
import shutil

def main():
    parser = argparse.ArgumentParser(
        description="Restructure input folders"
    )
    parser.add_argument("-i", "--input", dest="input", help="Step input")
    parser.add_argument("-o", "--output", dest="output", help="Final output directory")

    args = parser.parse_args()
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    for root, dirs, files in os.walk(args.input):
        for filename in files:
            # Check if the file is a .fa file
            if filename.endswith('.fa'):
                # Construct the source and destination paths
                src_path = os.path.join(root, filename)
                dest_path = os.path.join(output_dir, filename)
                # Copy the file to the output directory
                shutil.copy(src_path, dest_path)
    print(len(os.listdir(output_dir)))

if __name__ == "__main__":
    main()
