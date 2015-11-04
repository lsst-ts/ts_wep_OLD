import luigi
import numpy as np
from astropy.io import fits
from CWFS.wcs import wcs
from SRCProcessingTask import SRCProcessingTask

class CWFSTask(luigi.Task):

	def output(self):
		return

	def requires(self):
		return SRCProcessingTask()

	def run(self):
		#paraxial test
		#    intra = "../TestImages/F1.23_1mm_v61/z7_0.25_intra.txt"
		#    extra ="../TestImages/F1.23_1mm_v61/z7_0.25_extra.txt"
		#    m = wcs(intra, 0, 0, extra, 0, 0, 'lsst.param', 'fft.algo', 'paraxial')

		#on axis test: FFT, z7
		#    intra = "../TestImages/LSST_C_SN26/z7_0.25_intra.txt"
		#    extra ="../TestImages/LSST_C_SN26/z7_0.25_extra.txt"
		#    m = wcs(intra, 0, 0, extra, 0, 0, 'lsst.param', 'fft.algo', 'onAxis')

		#on axis test: series expansion, z7
		#    intra = "../TestImages/LSST_C_SN26/z7_0.25_intra.txt"
		#    extra ="../TestImages/LSST_C_SN26/z7_0.25_extra.txt"
		#    m = wcs(intra, 0, 0, extra, 0, 0, 'lsst.param', 'exp.algo', 'onAxis')

		#off axis test: FFT, z11
		#    intra = "../TestImages/LSST_NE_SN25/z11_0.25_intra.txt"
		#    extra ="../TestImages/LSST_NE_SN25/z11_0.25_extra.txt"
		#    m = wcs(intra, 1.185, 1.185, extra, 1.185, 1.185, 'lsst.param', 'fft.algo', 'offAxis')

		#off axis test: series expansion, z11
		intra = "../test_images/LSST_NE_SN25/z11_0.25_intra.txt"
		extra = "../test_images/LSST_NE_SN25/z11_0.25_extra.txt"
		instruFile = "../conf/lsst.param"
		algoFile = "../conf/exp.algo"
		model = "offAxis"
		I1fldx = 1.185
		I1fldy = 1.185
		I2fldx = 1.185
		I2fldy = 1.185
		m = wcs(intra, I1fldx, I1fldy, extra, I2fldx, I2fldy, instruFile, algoFile, model)

		np.savetxt('python_zc.txt',m.converge[3:,-1])
