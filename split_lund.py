#!/usr/bin/env python3
"""
split_lund.py  --  Split a LUND file into smaller files.

Usage:
    python split_lund.py input.lund [-n 5000] [-o output_dir] [-prefix tag]

Arguments:
    input.lund      Input LUND file
    -n N            Events per output file (default: 5000)
    -o DIR          Output directory (default: same as input file)
    -prefix STR     Output filename prefix (default: input basename without .lund)

Output files are named:  <prefix>_part001.lund, <prefix>_part002.lund, ...
"""

import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('input',           help='Input LUND file')
parser.add_argument('-n',   default=5000, type=int, help='Events per output file')
parser.add_argument('-o',   default='',  help='Output directory')
parser.add_argument('-prefix', default='', help='Output filename prefix')
args = parser.parse_args()

infile  = args.input
n_per   = args.n
outdir  = args.o if args.o else os.path.dirname(os.path.abspath(infile))
prefix  = args.prefix if args.prefix else \
          os.path.splitext(os.path.basename(infile))[0]

os.makedirs(outdir, exist_ok=True)

def open_part(part_num):
    fname = os.path.join(outdir, f'{prefix}_part{part_num:03d}.lund')
    return open(fname, 'w'), fname

part_num  = 1
ev_count  = 0   # events in current file
total_ev  = 0
fout, fname = open_part(part_num)

with open(infile) as fin:
    while True:
        header = fin.readline()
        if not header:
            break                   # EOF
        tok = header.split()
        if not tok:
            continue
        try:
            npart = int(tok[0])
            if npart not in (6, 7):
                raise ValueError
        except (ValueError, IndexError):
            # Not a valid event header — skip line
            continue

        # Read the npart particle lines
        particles = []
        for _ in range(npart):
            line = fin.readline()
            if not line:
                break
            particles.append(line)

        # Write event to current output file
        fout.write(header)
        fout.writelines(particles)
        ev_count += 1
        total_ev += 1

        # Roll over to next file if chunk is full
        if ev_count >= n_per:
            fout.close()
            print(f'  Wrote {ev_count} events -> {fname}')
            part_num += 1
            ev_count  = 0
            fout, fname = open_part(part_num)

# Close last file (may have fewer than n_per events)
fout.close()
if ev_count > 0:
    print(f'  Wrote {ev_count} events -> {fname}')
else:
    os.remove(fname)   # empty file — remove it
    part_num -= 1

print(f'\nDone: {total_ev} events split into {part_num} files of <= {n_per} events each.')
