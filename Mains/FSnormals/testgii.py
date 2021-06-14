#! /usr/bin/env python

import sys
import numpy as np
import nibabel

print("nibabel.__version__ =", nibabel.__version__)

idx = np.arange(10)

da = nibabel.gifti.GiftiDataArray(idx,
                                  intent = 'NIFTI_INTENT_NODE_INDEX',
                                  datatype = 'NIFTI_TYPE_INT32')

gim = nibabel.gifti.gifti.GiftiImage(darrays = [da])
nibabel.save(gim, "idx.gii")

gim = nibabel.load("idx.gii")
print(gim.darrays[0].data)

da = nibabel.gifti.GiftiDataArray(idx, encoding = 'GIFTI_ENCODING_ASCII',
                                  intent = 'NIFTI_INTENT_NODE_INDEX',
                                  datatype = 'NIFTI_TYPE_INT32')

gim = nibabel.gifti.gifti.GiftiImage(darrays = [da])
nibabel.save(gim, "idx.gii")

gim = nibabel.load("idx.gii")
print(gim.darrays[0].data)
