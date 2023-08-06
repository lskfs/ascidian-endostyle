#!/usr/bin/env python3
"""usage: {} --image ssDNA.tif --outdir ./ --prefix out

options:
    -i, --image     <STR>       tiff image filename
    -o, --outdir    <STR>       directory saving the output, default ./
    -p, --prefix    <STR>       output prefix, default 'out'
    --cpi           <STR>       custom cellprofile pipeline file
"""
import os
os.environ['OPENCV_IO_MAX_IMAGE_PIXELS'] = pow(2, 40).__str__()
import sys
import getopt
from pathlib import Path
import pandas as pd
import numpy as np
from itertools import repeat
#from multiprocessing import set_start_method
#set_start_method('spawn')
#import multiprocessing
#from multiprocessing import get_context
import cv2
import skimage
import scipy.sparse

import bioformats.formatreader
import cellprofiler_core.pipeline
import cellprofiler_core.preferences
import cellprofiler_core.utilities.zmq
import cellprofiler_core.utilities.java
from cellprofiler_core.setting.subscriber import LabelSubscriber
from cellprofiler_core.setting.range import IntegerRange

def cp_main(image, name='DNA', cpi=None, skip=False, psize=None, saved_object='IdentifySecondaryObjects'):

    # load pipeline from cpi file
    sys.stdout.write(f'Load pipeline from {cpi}\n')
    sys.stdout.flush()
    pipeline = cellprofiler_core.pipeline.Pipeline()
    pipeline.load(cpi)

    # get modules list
    modules = pipeline.modules()

    # setup image_set
    image_set = cellprofiler_core.image.ImageSet(0, {'name':name}, name)
    x = cv2.imread(str(image), 0)
    x[x > 230] = 230
    image_x = cellprofiler_core.image.Image(x, path_name=image.parent, file_name=image.name)
    image_set.add(name, image_x)

    # init empty object_set
    object_set = cellprofiler_core.object.ObjectSet()

    measurements = cellprofiler_core.measurement.Measurements()

    workspace = cellprofiler_core.workspace.Workspace(
            pipeline, 
            modules,
            image_set,
            object_set,
            measurements,
            [image_set]
            )

    for module in modules:
        sys.stdout.write(f'\t... starting {module.module_name}\n')
        sys.stdout.flush()
        module.run(workspace)

    objects = workspace.object_set.get_objects(saved_object)
    try:
        celloutlines = workspace.image_set.get_image('CellOutlines')
    except:
        sys.stderr.write('cell outlines not get\n')
        celloutlines = None

    return objects, celloutlines

class Options:
    sargs = 'hi:n:o:p:m:'
    largs = ['help', 'image=', 'name=', 'outdir=', 'prefix=', 
            'cpi=', 'object=', 'mask=']
    cpi = str(Path(__file__).parent.resolve() / 'default.cppipe')
    
    def __init__(self, cmdline):
        self.cmdline = cmdline
        self.image = None
        self.image_name = 'DNA'
        self.outdir = Path('./').resolve()
        self.prefix = 'out'
        self.saved_object = 'IdentifySecondaryObjects'
        self.mask = None
        try:
            self.options, self.not_options = getopt.gnu_getopt(self.cmdline, 
                    self.sargs, self.largs)
        except getopt.GetoptError as err:
            sys.stderr.write('***Error: {}\n'.format(err))
            sys.exit()
        if not self.options:
            sys.stderr.write(self.usage)
            sys.exit()
        for opt, arg in self.options:
            if opt in ['-h', '--help']:
                sys.stderr.write(self.usage)
                sys.exit()
            elif opt in ['-i', '--image']:
                self.image = Path(arg).resolve()
            elif opt in ['-m', '--mask']:
                self.mask = arg
            elif opt in ['-n', '--name']:
                self.image_name = arg
            elif opt in ['-o', '--outdir']:
                self.outdir = Path(arg).resolve()
            elif opt in ['-p', '--prefix']:
                self.prefix = arg
            elif opt in ['--cpi']:
                self.cpi = arg
            elif opt in ['--object']:
                self.saved_object = arg
        
    @property
    def usage(self):
        return __doc__.format(__file__)

def main():
    opts = Options(sys.argv[1:])
    if not opts.outdir.exists():
        opts.outdir.mkdir(parents=True, exist_ok=True)

    if opts.image and not opts.mask:
        try:
            cellprofiler_core.preferences.set_headless()
            cellprofiler_core.preferences.set_temporary_directory(opts.outdir)
            cellprofiler_core.preferences.set_default_output_directory(opts.outdir)

            cellprofiler_core.utilities.java.start_java()

            print('Starting cellprofiler identify ...', flush=True)
            objects, celloutlines = cp_main(
                    opts.image, name=opts.image_name, 
                    cpi=opts.cpi, skip=opts.skip, psize=opts.pixel_range, 
                    saved_object=opts.saved_object
                    )
            print('Cell objects and outlines generated', flush=True)
        except Exception as err:
            sys.stderr.write('***Error: {}\n'.format(err))
        finally:
            cellprofiler_core.utilities.zmq.join_to_the_boundary()
            bioformats.formatreader.clear_image_reader_cache()
            cellprofiler_core.utilities.java.stop_java()

        print('Saving labled cells ...', flush=True)
        mask_file = opts.outdir / '{}_mask.txt'.format(opts.prefix)
        np.savetxt(mask_file, objects.segmented, fmt='%d')
        mask = objects.segmented
    elif opts.image and opts.mask:
        mask = np.loadtxt(opts.mask)
        celloutlines = None

    celloutlines_file = str(opts.outdir / f'{opts.prefix}_CellOutlines.txt')
    if celloutlines is not None:
        print('Saving cell outlines ...', flush=True)
        b, g, r = cv2.split(celloutlines.pixel_data)
        np.savetxt(celloutlines_file, b, fmt='%d')
    else:
        print('Detecting cell outlines ...', flush=True)
        import skimage.segmentation
        image = cv2.imread(str(opts.image))
        outlines = skimage.segmentation.mark_boundaries(
                image,
                mask,
                color=(1, 0, 0),
                mode='inner',
                )
        b, g, r = cv2.split(outlines)
        np.savetxt(celloutlines_file, b, fmt='%d')
    print('Saving masked image ...', flush=True)
    b = np.isin(b, [1])
    image = cv2.imread(str(opts.image))
    image[b] = (255, 0, 0)
    celloutlines_png = str(opts.outdir / '{}_CellOutlines.png'.format(opts.prefix))
    cv2.imwrite(celloutlines_png, image)
    
    return 

if __name__ == '__main__':
    main()

