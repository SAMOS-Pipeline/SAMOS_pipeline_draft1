#!/usr/bin/env python

"""
defineobs.py

An alternative Python version of defineobs
with a simple command line interface

Written by Edward Villanueva
Last updated 2014-02-27
"""

import os,sys

class obsdict:

      def __init__(self):
	  self.setDictOptions()

      def setDictOptions(self):
	  """Define dictionaries for fields with options"""
	  self.instrument_options = { "IMACS" : "IMACS", \
	  			      "LDSS2" : "LDSS2", \
				      "LDSS3" : "LDSS3" }

#	  self.camera_options = { "f/2 Mosaic1 (SITE)" : "SITE", \
#	                          "f/2 Mosaic1 (SITE2)" : "SITE2", \
#	                          "f/2 Mosaic2 (E2V)" : "E2V", \
#	                          "f/4 Mosaic1" : "SITE", \
#	                          "f/4 Mosaic3" : "Mos3" }

	  self.camera_options = { "SHORT f/2" : "SHORT", \
	                          "LONG f/4" : "LONG" }

	  self.mode_options = { "Direct" : "DIRECT", \
	  		        "Spectroscopic" : "SPEC" }

	  self.ns_options = { "Yes" : True, \
	  		      "No" : False }

	  self.grating_options = { '150 l'  : 'IMACS_150A', \
	                           '300 l'  : 'IMACS_300A', \
	                           '600 l'  : 'IMACS_600A', \
	                           '1200 l' : 'IMACS_1200A', \
	                           'MOE'    : 'MOE' }

	  self.imacs_grism_options = { '150 l'  : 'IMACS_grism_150', \
	                               '200 l'  : 'IMACS_grism_200', \
	                               '300 l'  : 'IMACS_grism_300', \
	                               '300R l' : 'IMACS_grism_300R', \
	                               '600 l'  : 'IMACS_grism_600' }

	  self.ldss2_grism_options = { 'low'       : 'LDSS3_low', \
	                               'med blue'  : 'LDSS3_med_blue', \
	                               'med red'   : 'LDSS3_med_red', \
	                               'high blue' : 'LDSS3_high_blue', \
	                               'high red'  : 'LDSS3_high_red' }

	  self.ldss3_grism_options = { 'med red'   : 'LDSS3_med_red', \
	                               'high blue' : 'LDSS3_high_blue', \
	                               'high red'  : 'LDSS3_high_red', \
	                               'VPH red'   : 'LDSS3_VPHred', \
	                               'VPH blue'  : 'LDSS3_VPHblue', \
	                               'VPH all'   : 'LDSS3_VPHall' }

      def generateMainDict(self):
	  """ keys to have for each field:
	      1. label = display name for field
	      2. options = choices to select from (sub-dict with value ?)
	      3. out_label = name of field output in obsdef file
	      4. value = store input/choice from user?
	  """

	  D = {}

	  D['date'] = { "label" : "Observation Date" }
	  D['instrument'] = { "label" : "Instrument", \
	                      "options" : self.instrument_options }
	  D['mask'] = { "label" : "Mask" }
	  D['dewoff'] = { "label" : "Dewar Offset File" }
	  D['camera'] = { "label" : "Camera", \
	                  "options" : self.camera_options }
	  D['mode'] = { "label" : "Mode", \
	                "options" : self.mode_options }
	  D['ns'] = { "label" : "Nod & Shuffle Mode", \
	              "options" : self.ns_options }
	  D['grating'] = { "label" : "Grating", \
	                   "options" : self.grating_options }
	  D['gr_order'] = { "label" : "Grating Order", \
	                    "value" : 0 }
	  D['gr_angle'] = { "label" : "Grating Angle", \
	                    "value" : 0 }
	  D['grism'] = { "label" : "Grism", \
	                 "options" : { "LDSS2" : self.ldss2_grism_options, \
	                      	       "LDSS3" : self.ldss3_grism_options, \
	                      	       "IMACS" : self.imacs_grism_options }
		       }
	  D['d_alignrot'] = { "label" : "Disperser Misalignment", \
	  		      "value" : 0.0 }
	  D['dewar'] = {}
	  D['obsdef'] = { "label" : "Observation Definition File" }

	  return D

class obsinfo:

      def __init__(self,obsdef=None,test=None,verb=None,help=None):
	  if help: self.printHelp()
	  self.yr = None
	  self.mo = None
	  self.test = test
	  self.verb = verb
	  self.obsdef = obsdef
	  self.printTips()
	  Dobj = obsdict()
          self.MD = Dobj.generateMainDict()
	  if self.obsdef: self.loadObsdef()
          main_fields = ['date','instrument','mask','dewoff']
          self.out_fields = ['date','instrument','mask','dewoff','camera', \
	  	             'mode','dewar','grating','gr_order', \
			     'gr_angle','d_alignrot']
	  if verb:
	     print self.MD
	  if self.test:
	     print main_fields
	  [self.getInput(field) for field in main_fields]
	  if self.MD["instrument"]["value"] == "IMACS":
	     [self.getInput(field) for field in ["camera","mode"]]
	     if self.MD["mode"]["value"] == "SPEC":
	        self.MD["gr_order"]["value"] = 1
	        if self.MD["camera"]["value"] == "SHORT":
		   next_fields = ["grism","ns","gr_order","d_alignrot","obsdef"]
	        else:
		   next_fields = ["grating","gr_order","gr_angle","d_alignrot","obsdef"]
	     else:
	        if self.MD["camera"]["value"] == "SHORT":
	           self.MD["grating"]["value"] = "IMACS_direct_grism"
		else:
	           self.MD["grating"]["value"] = "IMACS_direct"
		next_fields = ["obsdef"]
	  else:
	     self.getInput("mode")
	     self.MD["camera"]["value"] = ""
	     if self.MD["mode"]["value"] == "SPEC":
	        self.MD["gr_order"]["value"] = 1
	        next_fields = ["grism","gr_order","d_alignrot","obsdef"]
	     else:
	        self.MD["grating"]["value"] = "LDSS_direct_grism"
		next_fields = ["obsdef"]
	  if self.test:
	     print next_fields
	  [self.getInput(field) for field in next_fields]
	  self.writeObsdef()

      def printTips(self):
	  tips = "* Press <Ctrl+C> at any time to quit *"
	  print len(tips)*"*"
	  print tips
	  print len(tips)*"*"

      def getDate(self):
	  curyr,curmo = self.yr,self.mo
	  curval = self.checkCurrentVal("self.yr")
	  q = "Enter year of observation (e.g. 2012)%s: " % (curval)
	  self.yr = None
	  while not self.yr:
	        tmpval = raw_input(q)
	        if tmpval != "":
	           try:
	              self.yr = int(tmpval)
	     	      if self.yr < 2001:
	     	         print "Invalid year !"
	     	         self.yr = None
	     	         continue
	           except:
	     	      print "Please input an integer !"
	              continue
	        else:
		   if len(curval): self.yr = curyr
	  curval = self.checkCurrentVal("self.mo")
	  q = "Enter month of observation (e.g. '5' for May)%s: " % (curval)
	  self.mo = None
	  while not self.mo:
	        tmpval = raw_input(q)
	        if tmpval != "":
	           try:
	              self.mo = int(tmpval)
	           except:
	     	      print "Please input an integer !"
	              continue
	           if self.mo not in range(13)[1:]:
	              self.mo = None
	     	      print "Please type a number from 1 to 12 !"
	              continue
	        else:
		   if len(curval): self.mo = curmo
#	  self.MD["date"]["value"] = "%d-%02d" % (self.yr,self.mo)
	  self.MD["date"]["value"] = "%d/%d" % (self.mo,self.yr)
	  if self.test:
	     print "Setting %s to %s" % ("date",self.MD["date"]["value"])
	  print 20*"="

      def getInput(self,field):
	  if field == "date":
	     self.getDate()
	     return
	  if "options" not in self.MD[field].keys():
	     curval = self.checkCurrentVal(field)
	     q = "Enter %s%s: " % (self.MD[field]["label"],curval)
	     tmpval = ""
	     while tmpval == "":
	           tmpval = raw_input(q)
		   if "value" in self.MD[field].keys():
#		      tmpval = ""
		      break
	     if tmpval != "":
	        self.MD[field]["value"] = tmpval
	     if self.test:
		if tmpval == "":
	           print "Keeping %s as %s" % (self.MD[field]["label"],self.MD[field]["value"])
		else:
	           print "Setting %s to %s" % (self.MD[field]["label"],self.MD[field]["value"])
	  else:
	     ans = -1
	     opts = sorted(self.MD[field]["options"].keys())
	     if field == "grism":
	        opts = sorted(self.MD[field]["options"][self.MD["instrument"]["value"]].keys())
	     curval = self.checkCurrentVal(field)
#	     if self.test: print str(curval[2:-1])
	     if curval == "":
	        curind = ""
	     else:
	        curind = " (%d)" % (opts.index(str(curval)[2:-1])+1)
	     while ans not in range(len(opts)):
		   ans = -1
                   print "Choose %s:" % (self.MD[field]["label"])
	           for n,opt in enumerate(opts):
	               print "(%d) %s" % (n+1,opt)
		   q = "Enter number%s: " % (curind)
		   ans = raw_input(q)
		   try:
	              ans = int(ans) - 1
		   except:
		      if curval != "" and not len(ans):
			 ans = -1
		         break
		   if ans not in range(len(opts)):
	              print "Invalid choice! Try again."
	     if field == "grism":
		if ans > -1:
	           self.MD["grating"]["value"] = self.MD[field]["options"][self.MD["instrument"]["value"]][opts[ans]]
	           if self.test:
	              print "Setting %s to %s" % (self.MD["grating"]["label"],self.MD["grating"]["value"])
		else:
		   if self.test:
	              print "Keeping %s as %s" % (self.MD["grating"]["label"],self.MD["grating"]["value"])
	     else:
	        if ans > -1:
		   self.MD[field]["value"] = self.MD[field]["options"][opts[ans]]
	           if self.test:
	              print "Setting %s to %s" % (self.MD[field]["label"],self.MD[field]["value"])
	        else:
		   if self.test:
	              print "Keeping %s as %s" % (self.MD[field]["label"],self.MD[field]["value"])
	  print 20*"="

      def getDewar(self):
	  """ f/2: SITE if date < Aug 2005
	           SITE2 if Aug 2005 < date <= Mar 2008
	           E2V if > Mar 2008
	      f/4: SITE if date < Oct 2011
	           Mos3 if date >= Oct 2011
	  """
          if self.MD["instrument"]["value"] == "IMACS":
	     if self.MD["camera"]["value"] == "SHORT":
	        if self.dateToFloat() < 2005.08:
		   self.MD["dewar"]["value"] = "SITE"
	        else:
		   if self.dateToFloat() < 2008.04:
		      self.MD["dewar"]["value"] = "SITE2"
		   else:
		      self.MD["dewar"]["value"] = "E2V"
	     else:
	        if self.dateToFloat() < 2011.10:
		   self.MD["dewar"]["value"] = "SITE"
	        else:
		   self.MD["dewar"]["value"] = "Mos3"
          elif self.MD["instrument"]["value"] == "LDSS3":
	     self.MD["dewar"]["value"] = "LDSS3"
	     if self.dateToFloat() > 2014.03:
	        self.MD["dewar"]["value"] = "LDSS3C"
          else:
	     self.MD["dewar"]["value"] = self.MD["instrument"]["value"]
	  if self.test:
	     print "Setting dewar to",self.MD["dewar"]["value"]

      def writeObsdef(self):
#	  print self.MD["ns"]["value"],self.mo,self.yr
	  if self.MD["camera"]["value"] == "SHORT":
	     self.checkNSchange()
	  self.getDewar()
	  out = "%s.obsdef" % (self.MD["obsdef"]["value"])
	  print "Writing the following to %s:" % (out)
	  print
	  f = open(out,"w")
          for field in self.out_fields:
	      val = self.MD[field]["value"]
	      if test:
	         sys.stdout.write("%s "%(field));sys.stdout.flush()
	      print >> f, "%-15s" % (field.upper()),val
	      print "%-15s" % (field.upper()),val
	  f.close()
	  if test: print
	  print
	  print "Wrote",out

      def checkNSchange(self):
	  """After May 2007, Nod & Shuffle mode was done by
	     rotating the disperser instead of the dewar"""
	  if self.MD["ns"]["value"] and self.dateToFloat() > 2007.05:
	     if "_NS" not in self.MD["grating"]["value"]:
	        self.MD["grating"]["value"] += "_NS"

      def dateToFloat(self):
          fdate = float(self.yr) + (0.01*float(self.mo))
	  return fdate

      def checkCurrentVal(self,field):
	  getOpts = lambda f: self.MD[f]["options"]
	  getVal = lambda f: self.MD[f]["value"]
	  curval = ""
	  if field in self.MD.keys():
	     if "value" in self.MD[field].keys():
	        if "options" not in self.MD[field].keys():
	           curval = " (%s)" % (str(getVal(field)))
	        else:
		   if self.test:
		      print getOpts(field).keys()
		   curkey = [k for k in getOpts(field).keys() \
		   	       if getOpts(field)[k] == getVal(field)][0]
	           curval = " (%s)" % (curkey)
	     elif field == "grism" and "value" in self.MD["grating"].keys():
		inst = getVal("instrument")
		if self.test:
		   print getOpts(field)[inst].keys()
#		   print getVal("grating")
		curkey = [k for k in getOpts(field)[inst].keys() \
			    if getOpts(field)[inst][k] == getVal("grating").rstrip("_NS")][0]
#	        curval = " (%s)" % (str(getVal("grating")))
	        curval = " (%s)" % (curkey)
	     else:
	        curval = ""
	  elif field in ["self.yr","self.mo"]:
	     if eval(field): curval = " (%s)" % (str(eval(field)))
	  return curval

      def loadObsdef(self):
	  if ".obsdef" not in self.obsdef:
	     self.obsdef += ".obsdef"
	  if not os.path.exists(self.obsdef):
	     print "%s does not exist!" % (self.obsdef)
	     sys.exit()
          f = open(self.obsdef,"r")
	  ls = [l.strip() for l in f.readlines()]
	  f.close()
	  for l in ls:
	      try:
	         fld,v = l.split()
	      except:
	         fld,v = l,""
	      self.MD[fld.lower()]["value"] = v
	  if "value" in self.MD["date"].keys():
#	     self.yr,self.mo = map(int,self.MD["date"]["value"].split("-"))
	     self.mo,self.yr = map(int,self.MD["date"]["value"].split("/"))
	  if "_NS" in self.MD["grating"]["value"]:
	     self.MD["ns"]["value"] = True
	     self.MD["grating"]["value"] = self.MD["grating"]["value"].rstrip("_NS")
	  else:
	     self.MD["ns"]["value"] = False
	  self.MD["obsdef"]["value"] = self.obsdef.split(".")[0]
	  print "Loaded",self.obsdef

      def printHelp(self):
	  print 20*"#"
          print "Optional arguments:"
	  print 20*"#"
	  print "\t-o [obsdef name], --obsdef [obsdef name]\n\t\tload existing obsdef file"
	  print
	  print "\t-h, --help\n\t\tprint this help file"

	  print 10*"#"
	  print "Examples:"
	  print 10*"#"
	  print "\t%s" % (sys.argv[0])
	  print "\t%s -o mymask" % (sys.argv[0])
	  print "\t%s -o mymask.obsdef" % (sys.argv[0])
	  sys.exit()

if __name__ == "__main__":
   obsdef = None
   help = 0
   test = 0
   for n,arg in enumerate(sys.argv):
       if arg in ("-t","--test"): test = 1
       elif arg in ("-o","--obsdef"): obsdef = sys.argv[n+1]
       elif arg in ("-h","--help"): help = 1
   O = obsinfo(obsdef=obsdef,test=test,help=help)
