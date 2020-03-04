def set_aidpars_calibration(aidpars):
    #noao>twodspec>longslit>aidpars
    aidpars.reflist = 'linelist' #Reference coordinate list
    aidpars.refspec = '' #Reference spectrum
    aidpars.crpix = 695 #Coordinate reference pixel
    aidpars.crquad = 0 #Quadratic pixel distortion at reference pixel
    aidpars.cddir = 'sign' #Dispersion direction -- Best left alone
    aidpars.crsearch = 'INDEF' #Coordinate value search radius
    #INDEF translates to -0.1 which corresponds to a search radius of 10%
    #of the estimated dispersion range
    aidpars.cdsearch = 'INDEF' #Coordinate interval search radius
    aidpars.ntarget = 20 #Number of target features
    aidpars.npattern = 7 #Number of lines in patterns
    aidpars.nneighbors = 10 #Number of nearest neighbors in patterns
    aidpars.nbins = 6 #Maximum number of search bins
    aidpars.ndmax = 500 #Maximum number of dispersions to evaluate
    aidpars.aidord = 3 #Dispersion fitting order 3 = Quadratic
    aidpars.maxnl = 0.035 #Maximum non-linearity
    aidpars.nfound = 6 #Minimum number of lines in final solution
    aidpars.sigma = 0.05 #Sigma of line centering (pixels)
    aidpars.minratio = 0.1 #Minimum spacing ratio to use
    aidpars.rms = 0.25 #RMS goal (fwidths)
    aidpars.fmatch = 0.1 #Matching goal (fraction unmatched)
    aidpars.debug = '' #Print debugging information -- for developer
    return

def set_autoidentify_calibration(autoidentify):
    #noaotwodspec>longslit>autoidentify
    #autoidentify.images = None Images containing features to be identified
    autoidentify.crval = 8000
    autoidentify.cdelt = 1.175
    autoidentify.coordlist = 'linelist_full'#'linelists$fear.dat' # Coordinate list
    autoidentify.units = 'angstroms' #Coordinate units
    autoidentify.interactive = 'NO' # Examine identifications interactively?
    autoidentify.aidpars = 'aidpars'
    autoidentify.section = 'middle line' #Section to apply to 2D images
    autoidentify.nsum = 15 # Number of lines/columns/bands to sum in 2D/3D images
    autoidentify.ftype = 'emission' #Feature type
    autoidentify.fwidth = 15 # Feature width in pixels
    autoidentify.cradius = 5.0 #Centering radius in pixels
    autoidentify.threshold = 10 # Feature threshold for centering
    autoidentify.minsep = 2 #Minimum pixel separation
    autoidentify.match = 20 # Coordinate list matching limit
    autoidentify.function = 'spline3' # Coordinate function
    autoidentify.order = 3 # Order of coordinate function
    autoidentify.sample = '*' #Coordinate sample regions
    autoidentify.niterate = 0 #Rejection iterations
    autoidentify.low_reject = 3.0 #Lower rejection sigma
    autoidentify.high_reject = 3.0 #Upper rejection sigma
    autoidentify.grow = 0.0 #Rejection growing radius
    autoidentify.dbwrite = 'YES' # Write results to database?
    autoidentify.overwrite = 'YES' #Overwrite existing database entries?
    autoidentify.database = 'database' #Database in which to record feature data
    autoidentify.verbose = 'NO' # Verbose output?
    return

def set_identify_calibration(identify):
    #noao>twodspec>longslit>identify
    identify.section = 'middle line'
    identify.databas = 'database'
    identify.coordli = 'linelist_full'
    identify.nsum = 10
    identify.match =20
    identify.maxfeat = 30
    identify.zwidth = 100
    identify.ftype = 'emission'
    identify.fwidth = 15
    identify.cradius = 5
    identify.thresho = 10
    identify.minsep = 2
    identify.function = 'spline3'
    identify.order = 3
    identify.niterat = 0
    identify.low_rej = 3
    identify.high_re = 3
    identify.grow = 0
    return


def set_reidentify_calibration(reidentify):
    #noao>twodspec>longslit>reidentify
    reidentify.interac = 'no'
    reidentify.section = 'middle line'
    reidentify.newaps = 'Yes'
    reidentify.overrid = 'Yes'
    reidentify.refit = 'Yes'
    reidentify.trace = 'Yes'
    reidentify.step = 50
    reidentify.nsum = 15
    reidentify.shift = 0
    reidentify.search = 0
    reidentify.nlost = 10
    reidentify.cradius = 5
    reidentify.thresho = 0
    reidentify.addfeat = 'No'
    reidentify.coordli = 'linelist_full'
    reidentify.match = 30
    reidentify.maxfeat = 30
    reidentify.minsep = 2
    reidentify.databas = 'database'
    reidentify.verbose = 'no'
    return

def set_fitcoords_calibration(fitcoords):
    #noao>twodspec>longslit>fitcoords
    fitcoords.interac = 'No'
    fitcoords.combine = 'Yes'
    fitcoords.databas = 'database'
    fitcoords.functio = 'chebyshev'
    fitcoords.xorder = 3
    fitcoords.yorder = 3
    return


def set_transform(transform):
    #noao>twodspec>longslit>transform
    transform.database = 'database'
    transform.interpt = 'spline3'
    transform.xlog = 'No'
    transform.flux = 'Yes'
    return
