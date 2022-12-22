# orbital-maneuvers
  Calculate time and delta-v requirements for changing orbits within LEO

## Installation
  - Install python 3.11 or greater (3.10 down untested, but may work)
  - Install packages `$ py -m pip install -r requirements.txt`

## Run
  - Run `get.py` with credentials (space-track.org login) to save current debris data
  - `catch.py` is the primary calculation script.
    Ensure your current directory is `src` when running.
  - Constants at the bottom of `catch.py` can be modified
  - Other files are used to generate graphs
