import glob
import numpy as np
import os,sys
from PIL import Image

groupname = os.getcwd().split('/')[-1]
simname = groupname.split("grouped_")[-1]
sims = glob.glob(f"../{simname}-*")
N = len(sims)
in_str = ""

for sim in sims:
   basename = sim.split("../")[-1]
   print(basename)
   pngname = sim + "/" + f"{basename}--001.png"
   image = Image.open(pngname)
   width, height = image.size
   in_str += pngname+" "
in_str = in_str[:-1]
out_str = f"{groupname}.png"
os.system(f"montage -tile {N}x3 -geometry {width}x{height}+4+0 {in_str} {out_str}")
os.system(f"convert -trim {out_str} {out_str}")



