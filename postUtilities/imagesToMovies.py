"""
Copyright (c) 2026 Zhi Sheng Chen

Licensed under the PolyForm Noncommercial License 1.0.0 (the "License").
You may not use this file except in compliance with the License. Any
commercial use requires a separate commercial license from the copyright
holder.

You may obtain a copy of the License at:
    https://polyformproject.org/licenses/noncommercial/1.0.0

Or see the LICENSE file in the root of this repository.

This software is provided "as is", without warranty of any kind, express or
implied. See the License for the specific language governing permissions
and limitations.
"""

#makes one mp4 per image directory (postProcessing/images/<panel>) with ffmpeg
#fps from --fps, else [MOVIES] FPS in pvPostInternalSetup, else 10
#usage (in the case directory): python3 imagesToMovies.py [--fps N] [--imageDir dir] [--movieDir dir]

import os
import sys
import glob
import shutil
import argparse
import subprocess
import configparser

templateLoc = os.path.dirname(os.path.realpath(__file__))


def parseArgs(argv=None):
    parser = argparse.ArgumentParser(prog='imagesToMovies', description='mp4 videos from the post-processing images')
    parser.add_argument('--fps', type=float, default=None, help='frames per second (default: pvPostInternalSetup [MOVIES] FPS)')
    parser.add_argument('--imageDir', default=os.path.join('postProcessing', 'images'), help='one sub directory per video')
    parser.add_argument('--movieDir', default=os.path.join('postProcessing', 'movies'))
    parser.add_argument('--force', action='store_true', help='make the videos even if [MOVIES] ENABLE is False')
    return parser.parse_args(argv)


def getMovieSetup(caseDir):
    #getting the [MOVIES] settings from the case pvPostInternalSetup, else the default file
    setup = {'ENABLE': 'True', 'FPS': '10'}
    for path in (os.path.join(templateLoc, 'default', 'pvPostInternalSetup'), os.path.join(caseDir, 'pvPostInternalSetup')):
        if not os.path.isfile(path):
            continue
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(path)
        if config.has_section('MOVIES'):
            setup.update(dict(config['MOVIES'].items()))
    return setup


def findFrames(imageDir):
    #finding the frames of each panel, in file name order
    panels = {}
    for panelDir in sorted(glob.glob(os.path.join(imageDir, '*'))):
        if not os.path.isdir(panelDir):
            continue
        frames = sorted(glob.glob(os.path.join(panelDir, '*.png')) + glob.glob(os.path.join(panelDir, '*.jpeg'))
                        + glob.glob(os.path.join(panelDir, '*.jpg')))
        if frames:
            panels[os.path.basename(panelDir)] = frames
    return panels


def ffmpegCommand(listFile, fps, moviePath):
    #ffmpeg command, one image per frame at fps, even size for yuv420p
    return ['ffmpeg', '-y', '-loglevel', 'error', '-f', 'concat', '-safe', '0', '-i', listFile,
            '-vf', 'setpts=N/(%g*TB),scale=trunc(iw/2)*2:trunc(ih/2)*2' % (fps), '-r', '%g' % (fps),
            '-c:v', 'libx264', '-pix_fmt', 'yuv420p', moviePath]


def makeMovie(frames, fps, moviePath):
    #writing the frame list for ffmpeg
    listFile = moviePath + '.frames.txt'
    with open(listFile, 'w') as f:
        for frame in frames:
            f.write("file '%s'\n" % (os.path.abspath(frame).replace('\\', '/')))
    result = subprocess.run(ffmpegCommand(listFile, fps, moviePath), capture_output=True, text=True)
    os.remove(listFile)
    if result.returncode != 0:
        print('\t\tffmpeg failed for %s:\n%s' % (moviePath, result.stderr))
        return False
    return True


def main(argv=None):
    args = parseArgs(argv)
    caseDir = os.getcwd()
    caseName = os.path.basename(caseDir)
    print('\n\tMaking videos from %s...' % (args.imageDir))
    setup = getMovieSetup(caseDir)
    if setup['ENABLE'].lower() != 'true' and not args.force:
        print('\t\t[MOVIES] ENABLE is False, no videos')
        return 0
    fps = args.fps if args.fps is not None else float(setup['FPS'])
    if fps <= 0:
        sys.exit('ERROR! FPS must be positive')
    if shutil.which('ffmpeg') is None:
        sys.exit('ERROR! ffmpeg not found (load the ffmpeg module or install it)')
    panels = findFrames(args.imageDir)
    if not panels:
        print('\t\tno images in %s' % (args.imageDir))
        return 0
    os.makedirs(args.movieDir, exist_ok=True)
    failed = 0
    for panel, frames in panels.items():
        moviePath = os.path.join(args.movieDir, '%s_%s.mp4' % (caseName, panel))
        if makeMovie(frames, fps, moviePath):
            print('\t\t%s (%d frames, %g fps)' % (moviePath, len(frames), fps))
        else:
            failed += 1
    print('\t\t%d videos in %s%s' % (len(panels) - failed, args.movieDir, ', %d failed' % (failed) if failed else ''))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
