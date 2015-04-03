# This script submits each subprocess of adrian's dropseq.py script at 
# once using dependencies

import os
import sh
import time

#define job identifiers
now = time.strftime("%y%m%d-%H%M")
expt = "2KCL_b2"
tag = "%s_%s_" % (now, expt)
job = tag + "%s"

#path to scripts
scripts = "/home/man36/dbseq"

#path to log file
logs = os.path.join(
	"/groups/neuroduo/Aurel/dawnchorus/dawnchorus_data/dropseq/150402_b2/logs",
	job % "barcodeprocessing_LOG.txt")

#submit filtering job
bsub1 = sh.bsub.bake(H=True, J=job % "filter", 
	o=logs,
	W="10:00", q="priority", R="select[mem>=16000] && rusage[mem=16000]", N=True)

bsub1("python dropseq.py filter")

time.sleep(5)

bsub2 = sh.bsub.bake(J=job % "histo", 
	o=logs,
	w="done(%s_%s_filter)" % (date, expt),
	W="10:00", q="priority", R="select[mem>=16000] && rusage[mem=16000]", N=True)
bsub2("python dropseq.py histo")

time.sleep(5)

bsub3 = sh.bsub.bake(J=job % "choose_barcodes", 
	o=logs,
	w="done(%s_%s_histo)" % (date, expt),
	W="10:00", q="priority", R="select[mem>=16000] && rusage[mem=16000]", N=True)
bsub3("python dropseq.py choose_barcodes")

time.sleep(5)

bsub4 = sh.bsub.bake(J=job % "split", 
	o=logs,
	w="done(%s_%s_choose_barcodes)" % (date, expt),
	W="10:00", q="priority", R="select[mem>=16000] && rusage[mem=16000]", N=True)
bsub4("python dropseq.py split")

sh.bresume(J=job % "filter")