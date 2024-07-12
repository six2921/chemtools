#! /home/siu/anaconda/bin/python

import argparse

def count_dollars(filename):
    try:
        with open(filename, 'r', encoding='utf-8') as file:
            content = file.read()
            count = content.count('$$$$')
            return count
    except FileNotFoundError:
        print(f"File {filename} not found.")
        return None

def main():
    parser = argparse.ArgumentParser(description="Count occurrences of $$$$ in a text file.")
    parser.add_argument('filename', type=str, help='The name of the text file to read.')
    args = parser.parse_args()

    count = count_dollars(args.filename)
    if count is not None:
        print(f"The number of occurrences of $$$$ is: {count}")

if __name__ == '__main__':
    main()

