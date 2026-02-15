import os
import sys
import glob
import re

def main():
    generateMovies()


def extract_number(filename):
    # Find all sequences of digits in the string
    numbers = re.findall(r'\d+', filename)
    # Take the last sequence found and convert to int
    if numbers:
        return int(numbers[-1])
    return 0 # Fallback if no numbers found



def generateMovies():
    print('Getting slice directories')
    #getting slices directories
    sliceDirs = glob.glob(os.path.join('postProcessing','images','*slice*'))
    
    for dir in sliceDirs:
        print('\nCreating movies:')
        #gathering the prefixes
        #print(dir)
        slices = glob.glob(os.path.join(dir,'*.png'))
        prefixes = []
        for slice in slices:
            prefix = "_".join(slice.split('/')[-1].split('_')[:-2])
            if prefix in prefixes:
                continue
            else:
                prefixes.append(prefix)

        for prefix in prefixes:
            print('\t%s' % (prefix))
            images = glob.glob(os.path.join(dir,'%s*.png' % (prefix)))
            images = sorted(images, key=extract_number)
            
            
            output_txt_file = "%s/%s_images.txt" % (dir,prefix)
            with open(output_txt_file, "w") as f:
                for image in images:
                    print('\t\t%s' % (image))
                    f.write(f"file '{image}'\n")
                

            os.system("ffmpeg -y -framerate 10 -f concat -safe 0 -i '%s' -vf scale=1920:1080 -c:v libx264 -pix_fmt yuv420p %s.mp4 >> log.createMovies" % (output_txt_file,
                                                                                                                        os.path.join(dir,prefix)))
            
        


main()