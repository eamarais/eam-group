##############################################################################################
##############################################################################################
# Generate animated GIFs for research presentations using Python
# Avoid the reliance on some online tools. They may not be always safe. 
# And there are always limitations in maximum number of files, file size or pixel number or output file aspect ratio.
# When processing high resolution images, online tools can also be time consuming!
##############################################################################################
# Step 1: find the images needed for the GIFs

# move to the directory
import os
os.chdir("Y:\Study\Research_Data\TROPOMI\project_1_AF_SEAPORT")

# list all files within this directory
all_files = os.listdir()

# select your images
# there are many ways to select the file, the key is to play around the "string"
# first get the string for each file name, then use "str.startswith", "str.endswith" or if statements like: if "xxx" in str

# loop over all files, find the images that you need
image_files = []

for file in all_files:
    # get the strings for each file
    filename = os.fsdecode(file) 
    #  add the image to your list when the file name starts with "Cape Town_NO2"
    if filename.startswith('Cape Town_NO2') == True: 
        image_files.append(filename)

# sort the images to decide the order 
image_files.sort() 

# check if the files and the order are as expected
for file in image_files:
    print(file,sep='\n')
##############################################################################################
# Step 2: read each image and save out as GIFs

# Solution 1: "imageio"
import imageio

# read each image
images = [imageio.imread(file) for file in image_files]

# save out as a GIFs with duration of each frame specified in seconds
imageio.mimsave(os.path.join('CapeTown_TROPI_NO2_test_gif1.gif'),images,duration = 1) 
###################################
# Solution 2: "PIL" 
from PIL import Image  

img,*imgs = [Image.open(f) for f in image_files]
img.save(fp='CapeTown_TROPI_NO2_test_gif2.gif', format='GIF', append_images=imgs,save_all=True,duration=1,loop=0)
###################################
# Compare the speed (Soluiton 2 is faster than 1 on my PC)

import timeit
start = timeit.default_timer()
# your codes here
stop = timeit.default_timer()
print('Time: ', stop - start)  
##############################################################################################
# Remarks:
# both are very basic codes to generate gifs
# pros: easy, fast and safe (online tools need time to upload the files, and I personally think that your files are no longer yours as soon as you upload them to some unknow online tools)
# cons: little control of the output (image quality,output file size, output file pixel number...)

# Known issue: The color bars become inconsistent among the frames, but you can overlay a static color bar image in the powerpoint to cover this. 

# End
##############################################################################################
