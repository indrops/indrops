# This script submits each subprocess of adrian's dropseq.py script at 
# once using dependencies

import os
import sh

sh.bsub("python dropseq.py filter")


sh.bsub("python dropseq.py filter")


sh.bsub("python dropseq.py filter")

sh.bsub("python dropseq.py filter")