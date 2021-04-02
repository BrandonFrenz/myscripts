#!/usr/bin/env python
import argparse
import subprocess
import amino_acids

def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--workflow', help='The single chain workflow')
    parser.add_argument('-f', '--fasta', help='The fasta file')
    parser.add_argument('-b', '--bucket_name', help='The name of the bucket')
    args = parser.parse_args()
    return args

def parse_fasta(fasta_file):
    sequence = ''
    with open(fasta_file, 'r') as inf:
        for line in inf:
            if line.startswith('>'):
                continue
            else:
                sequence+=line.strip()
    return sequence

def validate_sequence(sequence):
    for aa in sequence:
        if aa not in amino_acids.one_letter_names.keys():
            return False
    return True

def submit_argo(workflow, bucket, sequence):
    command = ['argo', 'submit', workflow, '-p', f'bucket-name={bucket}', '-p', f'sequence={sequence}']
    subprocess.call(command)

def main():
    args = parseargs()
    sequence = parse_fasta(args.fasta)
    if not validate_sequence(sequence):
        print('The sequence is invalid. Check for non standard characters')
        exit()
    submit_argo(args.workflow, args.bucket_name, sequence)

if __name__ == '__main__':
    main()
