import glob
import numpy as np
import os,sys
from PIL import Image


def permute( rin, Ni=7, Nj=4 ):
   r = []
   for xi in range(Ni):
      for yi in range(Nj):
         r.append( rin[(Ni-xi-1)+Nj*yi] )
   return r

files = glob.glob("grouped_sample_*/*.png")
files.sort()
#files = permute(files)

file_str = " ".join(files)
im0name = files[0]
image = Image.open(im0name)
width,height = image.size
in_str = file_str
out_str = "montage_ss.png"
os.system(f"montage -tile 4x7 -geometry {width}x{height}+10+10 {in_str} {out_str}")
os.system(f"convert -trim -flip -flop {out_str} {out_str}")
os.system(f"convert -rotate 90 {out_str} {out_str}")


