#!/usr/bin/env python2.7
from subprocess import Popen, call
from string import join

def FormatSeqForScapBackbone(rawSeq):
  ## This function takes a sequence in Youval's format,
  ## i.e.PLSRKH|E(R)HVGDL|GNVTAD|-KDGVA|-DVSIE, records
  ## the positions of loops (*) and cuts -. Then loops
  ## are removed and Glycine is put in for cut,-, symbols.
  ## An example is given below:
  ## input:  rawSeq = "PLSRKH|E(R)HVGDL|GNVTAD|-KDGVA|-DVSIE"
  ## output: formatSeq = "PLSRKHEHVGDLGNVTADGKDGVAGDVSIE"
  ## 		   loopcutData = [[7,'R'],[19,0],[25,0]]
  ##  	   startPos = 1

  loopcutData = []
  formatSeq = ''
  loopflag = 0
  i = 0
  firstTurn = 0
  startPos = 1
  looplength = 0
  seqtemp = ''
  k = 0
  valLoop = []
  
  j = 0
  while firstTurn == 0:
    if rawSeq[j] == "|":
      firstTurn = 1
    else:
      startPos = 6 - j
    j += 1

  
  ## Loop through chars in the input string.
  for x in rawSeq:
    ## If NOT currently in a loop.
    if loopflag == 0:
      if x == '|':
        pass
      elif x == '-':
		## Denote cut with 0. Insert glycine
		## placeholder.
        loc = i
        loopcutData.append([loc,0])
        #formatSeq += 'G'
        #i+=1  ## THIS IS THE TROUBLE LINE.
        continue
      elif x == '(':
		## Denote the beginning of new loop.
		## Outermost if/else will now hit the else.
        temploop = ''
        loopflag = 1
        loc = i
        continue
      else:
        if (k + startPos)%6 == 3 or (k + startPos)%6 == 5:
          formatSeq += 'V'
          seqtemp += 'V'
          valLoop.append([k,'G'])
          i += 1
          k += 1
        else: 
          formatSeq += 'G'
          seqtemp += 'G'
          i += 1
          k += 1

    elif loopflag == 1:
	  ## This is entered when residues are in a loop.
      if x == ')':
		## Reached the end of a loop. Will now 
		## hit the outermost 'if'.

        loopflag = 0
        loopcutData.append([loc,temploop])
      else:
        temploop += 'G'
        #seqtemp += x
        i += 1

  #loopcutData.extend(valLoop)

  return formatSeq, seqtemp, startPos, loopcutData

def FormatSeqForScap(rawSeq):
  ## This function takes a sequence in Youval's format,
  ## i.e.PLSRKH|E(R)HVGDL|GNVTAD|-KDGVA|-DVSIE, records
  ## the positions of loops (*), cuts -, and mutations
  ## Then loops are removed and Glycine is put in for 
  ## cut,-, symbols. A flag is given for cut=0,loop=1,
  ## or mutation=2.
  ## An example is given below:
  ## input:  rawSeq = "PLSRKH|E(R)HVGDL|GNVTAD|-KDGVA|-DVSIE"
  ## output: formatSeq = "PLSRKHEHVGDLGNVTADGKDGVAGDVSIE"
  ## 		   loopcutData = [[1,7,'R'],[0,19,0],[0,25,0]]
  ##  	   startPos = 1

  loopcutData = []
  formatSeq = ''
  loopflag = 0
  i = 0
  firstTurn = 0
  startPos = 1
  looplength = 0
  seqtemp = ''
  
  ## Loop through chars in the input string.
  for x in rawSeq:
    ## If NOT currently in a loop.
    if loopflag == 0:
      if x == '|':
		## Assigns startPos at first turn symbol.
		## Then ignore all subsequent '|'s.
        if firstTurn == 0:
          startPos = 7 - len(formatSeq)
          firstTurn = 1
        else:
          continue
      elif x == '-':
		## Denote cut with 0. Insert glycine
		## placeholder.
        loc = i
        loopcutData.append([0,loc,0])
        formatSeq += 'G'
        #i+=1  ## THIS IS THE TROUBLE LINE.
        continue
      elif x == '(':
		## Denote the beginning of new loop.
		## Outermost if/else will now hit the else.
        temploop = ''
        loopflag = 1
        loc = i
        continue
      else:
        formatSeq += x
        seqtemp += x
        i += 1
    elif loopflag == 1:
	  ## This is entered when residues are in a loop.
      if x == ')':
		## Reached the end of a loop. Will now 
		## hit the outermost 'if'.

        loopflag = 0
        loopcutData.append([1,loc,temploop])
      else:
        temploop += x
        #seqtemp += x
        i += 1

  return formatSeq, seqtemp, startPos, loopcutData
     
def MakeLoopyList(thread, loopdata):

  temp = []
  length = len(thread)
  for i in range(length-2):
    temp.append([2, i, thread[i]])
  temp.append([3,i+1,thread[-2]])
  temp.append([4,i+2,thread[-1]])

  temp.extend(loopdata)
  print temp
  return temp
     
def MakeVariableLengthTemplate(length,startingPos):
  ## MasterTemplate.pdb starts at L2. The 'length' variable
  ## is the sequence length not counting residues in loops.
  ## The 'startingPos' variable is the desired starting 
  ## position {L1,L2,L3,L4,L5,L6}. 
  ## Example:
  ## 	input:  length = 46					startingPos = 1
  ##	output:	gname = 'tplt_46_L3.pdb'	begin = 5

  fname = '/home/alex/projects/Threader/MasterTemplate.pdb'
  gname = '/home/alex/projects/Threader/' 
  gfile = 'tplt_' + str(length) + '_L' + str(startingPos) + '.pdb'
  gpath = gname + gfile
  f = open(fname, 'r')
  g = open(gpath, 'w')

  ## atomLines trims off header and footer;
  ## it only has the coordinate data.
  
  allLines = f.readlines()
  atomLines = allLines[:-2]
  
  resLines = []
  previndex = 0
  tempRes = []
  
  ## This chunk of code sorts the lines into residues.
  ## This is redundant and a little wasteful to do 
  ## every time, so this could be a point to improve.
  for k in range(len(atomLines)):
    temp2 = atomLines[k]
    resindex = int(temp2[23:27])
  
    if resindex == previndex:
      tempRes.append(atomLines[k])
    else:
      previndex = resindex
      resLines.append(tempRes)
      tempRes = []
      tempRes.append(atomLines[k])
  
  ## Adjusting the starting position.
  begin = (startingPos-1+3) % 6
  
  ## Lines written to file 'g'.
  k = 0
  m = 1
  for res in resLines[begin:begin+length]:
    for line in res:
      resNum = str(k)
      atomNum = str(m)
      newAtomLine = "%s%6s%s%3s%s\n" % (line[0:5],atomNum,line[11:23], resNum, line[27:])
      g.write(newAtomLine)
      m +=1
    k += 1
    
  f.close()
  g.close()

  return gfile

def ThreadToScapFile(thread, protName):
  ## This takes the formatted thread from FormatSeqForScap,
  ## the starting index from MakeVariableLengthTemplate,
  ## and the protein name. This function simply creates
  ## a text file used in the scap program. Returns
  ## formatfname, which has a funky extension.
  
  formatfname = protName + '_scap.list'
  f = open(formatfname,"w")
  
  i = 0
  for x in thread:
      f.write("A," + str(i) + ',' + x + '\n')
      i += 1
  f.close()
  
  return formatfname
    
def RunScap(fname,tmpltname):
  ## This calls the Jackal program Scap with the python subprocess
  ## module. Scap mutates the structure template, tmpltname, using
  ## the residue list in fname. Parameters for Scap are currently 
  ## hardwired in, however it wouldn't be hard to change this in 
  ## the future.
  prm = '2'
  ini = '10'
  min = '2'
  rtm = '1'
  

  scaplistfile = (fname.split('/'))[-1]
  task = "scap -prm " + prm + " -min " + min + " -rtm " + rtm + " -ini " + ini + " " + tmpltname + " " + scaplistfile
  newtask = task.split()
  call(newtask)
  
  temp = tmpltname.split('.')
  newfname = temp[0] + '_scap.pdb'

  newtask = task.split()
  call(newtask)

  err = CheckForSuccess(newfname,fname)
    
  call(['rm',fname, tmpltname])
 
  
  return newfname, err
 
def RunLoopy(fname, loopData, seq, seq2,mastseq, structNum, name, fdest, failList):
  ## This function runs the Jackal program Loopy using the Python
  ## subprocess module. Loopy takes a structure and inserts loops.
  ## The output is a new file with the added extension '.loopy'.
  prm = '3'
  ini = '2000' 
  loopcount = 0
  cleanup = ["rm"]
  attempts = 1
  endfile = name + '.pdb'
  m = 1
  #print loopData, '\n' # TEMPORARY
  
  #call(['rm','tplt*']) # TEMPORARY
  
  orfname = fname

  for flag, loc, loop in loopData:
    loopcount = AddLoopCut(flag, loc,loop,loopcount, seq2,mastseq, prm, ini, fname)

    temp = (fname.split('.pdb'))[0] + '_loopy.pdb'
    
    err = CheckForSuccess(temp, name) # TEMPORARY
    #err = 0 # DEBUGGING

    if err == 0:
      cleanup.append(fname)
      fname = orfname[:-4] + "_" + str(m)  +'_loopy.pdb'
      call(["mv",temp,fname])
      print fname
    else:
      call(['rm','tplt*'])  # TEMPORARY
      break
    m += 1


  cleanup.append(fname)
  #err = 1 # DEBUGGING
  if err == 0:
    cutCall = "cut -c0-55 " + fname + " > " + endfile 
    call(cutCall, shell=True)

    mvtask = "mv " + endfile + " " + fdest
    call(mvtask.split())
    
    call(cleanup)
  else:
    pass

  return endfile, err, failList


def AddLoopCut(flag,loc,loop,loopcount, seq2, mastseq, prm, ini, fname):

  loc2 = loc - loopcount
  if flag == 0:
    ## This is entered when a cut needs to be entered.
    print "Now inserting a cut at residue ", loc
    loctemp = str(loc - 2) + "-" +  str(loc + 2)
    restemp = seq2[loc2 - 2: loc2] + seq2[loc2 : loc2 + 2]

  elif flag == 1:
    print "Now inserting a loop of ", loop, " at residue ", loc
    loctemp = str(loc - 2) + "-" + str(loc + 1)
    restemp = seq2[loc2 - 2: loc2 ] + loop + seq2[loc2 : loc2 + 2]
    loopcount += len(loop)

  elif flag == 2 or flag == 3 or flag == 4:
    print "Now mutating residue ", loc, " to ", loop
    if loc == 0:
      loctemp = str(loc+1) + "-" + str(loc + 3)
      restemp = mastseq[loc2+1] + loop + mastseq[loc2 + 2: loc2 + 4]

    elif loc == 1:
      loctemp = str(loc) + "-" + str(loc + 3)
      restemp = mastseq[loc2 - 1] + seq2[loc2 - 1] + loop + mastseq[loc2 + 1]

    elif flag == 4:
      loctemp = str(loc) + "-" + str(loc + 3)
      restemp = seq2[loc2 - 2: loc2 ] + loop + mastseq[loc2+1]

    else:
      loctemp = str(loc + 1) + "-" + str(loc + 3)
      restemp = seq2[loc2 - 1] + loop + mastseq[loc2 + 1]

  task = "loopy -prm " + prm + " -ini " + ini + " -obj " + loctemp + " -cid A -res " + restemp + " " + fname
  
  print task, '\n'
  newtask = task.split()
  call(newtask)

  return loopcount

def CheckForSuccess(temp, name):

  err = call(['tail', temp])

  if err == 1:
    print "Operation unsuccessfull. Skipping sequence " + name
  else:
    pass
 
  return err


def ReformatSeq(seq, start):
  ## This functin justs takes a thread being processed and puts the 
  ## corner symbols back in for readability.
  i = start
  new = ''
  loopflag = 0
  for x in seq:
    if loopflag == 1:
      new += x
      if x == ')':
        loopflag = 0
      else:
        pass
    else:
      if i == 6:
        new += x + '|'
        i = 1
      elif x == '(':
        loopflag = 1
        new += x 
      else:
        new += x
        i += 1

  return new

def ThreadByScapThenLoopy(thread,name, structNum, fdest,failList):

  mastseq = 'IDKSAFVHPTAIVEEGASIGANAHIGPFCIVGPHVEIGEGTVLKSHVVVNGHTKIGRDNEIYQFASIGGGVEIGDRNRIRESVTIGGGGVVGSDNLLMINAHIAHDCTVGNRCILANNATLAGHVSVDDFAIIGGMTAVHQFCII'
  formatSeq, seqtemp, startPos, loopData = FormatSeqForScap(thread)
  template = MakeVariableLengthTemplate(len(formatSeq),startPos)
  scapList = ThreadToScapFile(formatSeq, name)
  scapfile, err = RunScap(scapList, template) 

  if err == 1:
    finalfile = '0'
  elif err == 0 and loopData != []:
    finalfile, err2, failList = RunLoopy(scapfile,loopData ,formatSeq, seqtemp,mastseq, structNum, name, fdest, failList)
  elif err == 0 and loopData == []:
    finalfile = scapfile

  return finalfile, failList


def ThreadByLoopy(thread,name,structNum,fdest,failList):

  mastseq = 'IDKSAFVHPTAIVEEGASIGANAHIGPFCIVGPHVEIGEGTVLKSHVVVNGHTKIGRDNEIYQFASIGGGVEIGDRNRIRESVTIGGGGVVGSDNLLMINAHIAHDCTVGNRCILANNATLAGHVSVDDFAIIGGMTAVHQFCII'
  formatSeq, seqtemp, startPos, loopData = FormatSeqForScap(thread)
  newLoopData =  MakeLoopyList(formatSeq,loopData)

  if startPos == 0:
    startPos = 5
  elif startPos == 1:
    startPos = 6
  else:
    startPos = startPos - 2
  templateName = MakeVariableLengthTemplate(len(thread) + 2,startPos)
  finalfile,err2,failList =  RunLoopy(templateName, newLoopData, formatSeq, seqtemp, mastseq, structNum,name, fdest,failList)


  return finalfile,failList


def ThreadToStructure(thread, name, structNum, fdest, failList):
  ## This function is a wrapper for the entire process. It simply takes in 
  ## a sequence and the protein name. The name is just used to rename the
  ## final structure file.

  #formatSeq, seqtemp, startPos, loopData = FormatSeqForScapBackbone(thread)

  #finalfile, failList = ThreadByLoopy(thread,name,structNum,fdest,failList)
  finalfile, failList = ThreadByScapThenLoopy(thread,name,structNum,fdest,failList)

  return finalfile, failList

def GetInputThreads():
  ## This function retreives threads from the input file. It ignores any
  ## lines starting with a '#' and blank lines.
  inputFile = open("ThreadFile.txt","r")

  threadList = []
  for line in inputFile.readlines():
    if line == '' or line == '\n' or line == 'EOF' or line[0] == '#':
      #print 'pass'
      pass
    else:
      name = (line.split())[0]
      thread = (line.split())[1]
      threadList.append([name,thread])

  inputFile.close()

  return threadList


if __name__ == "__main__":
  ## The threads should be in a file named 'ThreadFile.txt', each line 
  ## formatted as follows:
  ## name GGGGGG|TTTTTT|(RR)DJKLKD|GGGGGG

  threadList = GetInputThreads()
  structNum = 1
  fdest = "/home/alex/projects/threads/LHBYeast"
  #fdest = "/export/home/shared/yeast"
  failList = []

  for name,thread in threadList:
    print "Target: ", thread
    fname, failList = ThreadToStructure(thread,name, structNum, fdest, failList)
    #thread = StripSequenceToBackbone(thread)
    print "Done with ", name, ". Located in file: ", fname
    structNum += 1
    print "*******************************************************************************\n"

