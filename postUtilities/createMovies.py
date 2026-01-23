import os
import sys
import glob


def main():
    generateMovies()


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
            os.system("ffmpeg -y -framerate 10 -pattern_type glob -i '%s*.png' -vf scale=1920:1080 -c:v libx264 -pix_fmt yuv420p %s.mp4 >> log.createMovies" % (os.path.join(dir),
                                                                                                                        os.path.join(dir,prefix)))
            break
        break


main()