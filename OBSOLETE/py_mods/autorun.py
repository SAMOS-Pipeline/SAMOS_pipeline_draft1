from InitializeSAMOS import initialize_SAMOS
from OverscanAndTrim import overscan_and_trim
from NormDivFlats import norm_div_flats
from OutlineSlits import outline_slits
from PyrafIdentify import pyraf_identify
import PyrafTransform



pyraf_trans = PyrafTransform.pyraf_trnsf
def auto_run():

	fuel_init = initialize_SAMOS('2017-11-30','LMask2')
	fuel_ost = overscan_and_trim(fuel_init)
	fuel_ndf = norm_div_flats(fuel_ost)
	fuel_os = outline_slits(fuel_ndf)
	fuel_id = pyraf_identify(fuel_os)
	fuel_tr = pyraf_trans(fuel_id)

	return fuel_init,fuel_ost,fuel_ndf,fuel_os,fuel_id,fuel_tr


