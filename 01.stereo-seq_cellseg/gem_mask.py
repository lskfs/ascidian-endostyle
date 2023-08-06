"""usage: {} -g <input.gem> -m <mask.txt> -p <prefix> [-a <affine.txt>] [-i <image>] [--binary] [--crop]

options:
    -g, --gem       <STR>       gem file name
    -m, --mask      <STR>       mask matrix TXT file used to segmente gem file
    -p, --prefix    <STR>       output file name prefix
    -a, --affine    <STR>       affine tranform matrix
    -i, --image     <STR>       image
    --binary        <BOOL>      force set mask matrix as binary mask
"""
#    --crop          <BOOL>      crop mask and image to exclude blank region

import os
import sys
import getopt
import numpy as np
sys.path.append('/dellfsqd2/ST_OCEAN/USER/hankai/software/SpatialTranscript/site-packages')

from stomics.gem import Gem
from stomics.mask import Mask

class Options:
    sargs = 'hg:m:i:a:p:'
    largs = ['help', 'gem=', 'mask=', 'prefix=', 'image=', 'affine=', 'crop', 'binary']
    def __init__(self, cmdline):
        self.cmdline = cmdline
        self.gem = None
        self.mask = None
        self.prefix = 'out'
        self.image = None
        self.crop = False
        self.affine = None
        self.binary = False
        self.invert = False
        try:
            self.options, self.not_options = getopt.gnu_getopt(self.cmdline,
                    self.sargs, self.largs)
        except getopt.GetoptError as err:
            sys.stderr.write(f'*** Error: {err}\n')
            sys.exit()
        if (not self.options) or self.not_options:
            sys.stderr.write(self.usage)
            sys.exit()
        for opt, arg in self.options:
            if opt in ['-h', '--help']:
                sys.stderr.write(self.usage)
                sys.exit()
            elif opt in ['-g', '--gem']:
                self.gem = arg
            elif opt in ['-m', '--mask']:
                self.mask = arg
            elif opt in ['-p', '--prefix']:
                self.prefix = arg
            elif opt in ['-i', '--image']:
                self.image = arg
            elif opt in ['-a', '--affine']:
                self.affine = arg
            elif opt in ['--binary']:
                self.binary = True
            elif opt in ['--invert']:
                self.invert = True
            elif opt in ['--crop']:
                self.crop = True
    @property
    def usage(self):
        return __doc__.format(__file__)

def main():
    opts = Options(sys.argv[1:])
        
    gem = Gem.readin(opts.gem)
    heatmap = gem.to_img()
    #print(heatmap.shape)

    mask = Mask(opts.mask)

    if opts.image:
        image = Mask(opts.image, mode=0)
        assert image.matrix.shape == mask.matrix.shape, 'mask and image should have the same shape'

    if opts.affine:
        affine = np.loadtxt(opts.affine)
        mask = mask.affine_transform(
                affine, 
                reference=heatmap,
                outfile=f'{opts.prefix}.transformed_mask.txt'
                )
        if opts.image:
            image = image.affine_transform(
                    affine,
                    reference=heatmap,
                    outfile=f'{opts.prefix}.transformed.tif'
                    )
            border = mask.overlayoutlines(
                    image=image.matrix,
                    prefix=f'{opts.prefix}.transformed'
                    )
    
    count = len(np.unique(mask.matrix))
    if count > 2:
        sys.stdout.write(f'detect mask as labeled with {count-1} objects\n')
        label_object=True
    else:
        sys.stdout.write('detect mask as binary\n')
        label_object=False
    if opts.binary:
        sys.stdout.write('force using mask as binary with --binary\n')
        label_object=False

    gem, (offset_x, offset_y) = gem.mask(
            mask.matrix, 
            outfile=f'{opts.prefix}.gem', 
            label_object=label_object,
            )
    sys.stderr.write(f'offsetx={offset_x} offsety={offset_y}\n')
    
    """
    if opts.crop:
        heatmap = gem.to_img()
        print(heatmap.shape)
        if opts.image:
            #print(image.matrix.shape)
            image = image.crop(
                    mask=mask.matrix, 
                    outfile=f'{opts.prefix}.cropped.tif'
                    )
            print(image.matrix.shape)
            #assert heatmap.shape == image.matrix.shape
        #print(mask.matrix.shape)
        mask = mask.crop(
                outfile=f'{opts.prefix}.cropped_mask.txt'
                )
        print(mask.matrix.shape)
        #assert heatmap.shape == mask.matrix.shape
    """
    return

if __name__ == '__main__':
    main()


