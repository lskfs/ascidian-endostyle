#!/usr/bin/env python3
"""usage: {} [options] -h

options:
    --mask          <STR>       label mask
    --bam           <STR>       bam file
    --out           <STR>       out.bam

    --shiftx        <INT>       mask2gem shift X
    --shifty        <INT>       mask2gem shift Y
        coords in mask + shift == coords in gem, you can get
        those two from the result of registration

    --offsetx       <INT>       gem2bam offset X
    --offsety       <INT>       gem2bam offset Y
        coords in gem + offset == coords in bam, you can get 
        those two from the header of original GEM file

    --xmin, --xmax
    --ymin, --ymax
        contour box positions

NOTE: you should always kown the meaning of shift and offset!!!
    MUST be caution to use, otherwise you will get wrong data
"""

import os
import sys
import getopt
import numpy as np
import pandas as pd
import pysam

class Options:
    sargs = 'h'
    largs = ['help', 'mask=', 'bam=', 'out=', 'offsetx=', 'offsety=', 'shiftx=', 'shifty=', 
            'xmin=', 'xmax=', 'ymin=', 'ymax=']
    def __init__(self, cmdline):
        self.cmdline = cmdline
        self.label_mask = None
        self.bam_file = None
        self.out_file = None

        self.shiftx = 0
        self.shifty = 0
        self.offsetx = 0
        self.offsety = 0

        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        try:
            self.options, self.not_options = getopt.gnu_getopt(self.cmdline, 
                    self.sargs, self.largs)
        except getopt.GetoptError as err:
            sys.stderr.write(f'***Error: {err}\n')
            sys.exit()
        if not self.options:
            sys.stdout.write(self.usage)
            sys.exit()
        for opt, arg in self.options:
            if opt in ['-h', '--help']:
                sys.stdout.write(self.usage)
                sys.exit()
            elif opt in ['--mask']:
                self.label_mask = arg
            elif opt in ['--bam']:
                self.bam_file = arg
            elif opt in ['--out']:
                self.out_file = arg
            elif opt in ['--shiftx']:
                self.shiftx = int(arg)
            elif opt in ['--shifty']:
                self.shifty = int(arg)
            elif opt in ['--offsetx']:
                self.offsetx = int(arg)
            elif opt in ['--offsety']:
                self.offsety = int(arg)
            elif opt in ['--xmin']:
                self.xmin = int(arg)
            elif opt in ['--xmax']:
                self.xmax = int(arg)
            elif opt in ['--ymin']:
                self.ymin = int(arg)
            elif opt in ['--ymax']:
                self.ymax = int(arg)
    @property
    def usage(self):
        return __doc__.format(__file__)

def mask2dict(mask, offsetx=0, offsety=0, shiftx=0, shifty=0):
    xy2label = dict()
    for (y, x), label in np.ndenumerate(mask):
        if label == 0:
            continue
        x = x + offsetx + shiftx
        y = y + offsety + shifty
        xy2label[(x, y)] = str(int(label))
    return xy2label

def main():
    opts = Options(sys.argv[1:])

    if opts.label_mask:
        label_mask = np.loadtxt(opts.label_mask)
    
        xy2label = mask2dict(label_mask, offsetx=opts.offsetx, offsety=opts.offsety, 
                shiftx=opts.shiftx, shifty=opts.shifty)
    else:
        xy2label = None

    if opts.xmin and opts.xmax and opts.ymin and opts.ymax:
        xmin, xmax = opts.xmin, opts.xmax
        ymin, ymax = opts.ymin, opts.ymax

    sam = pysam.AlignmentFile(opts.bam_file, 'rb')
    out_bam = pysam.AlignmentFile(opts.out_file, 'wb', template=sam)
    for read in sam.fetch():
        if not read.has_tag('GE'):
            continue
        if read.is_duplicate:
            continue
        if read.is_qcfail:
            continue
        x = read.get_tag('Cx')
        y = read.get_tag('Cy')
        
        if xy2label is not None:
            if not xy2label.get((x, y)):
                continue
            label = xy2label[(x, y)]
            read.set_tag('CB', label)
        else:
            if not (xmin <= x <= xmax and ymin <= y <= ymax):
                continue
        out_bam.write(read)

    out_bam.close()
    sam.close()

if __name__ == '__main__':
    main()

