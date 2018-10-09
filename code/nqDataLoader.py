# -*- coding: utf-8 -*-

# set modules  dir
import numpy as np
import sys, os, re, datetime


class NqDataLoader:
    FLT_NO_MOUSE = 1 << 0
    FLT_NO_LETTERS = 1 << 1
    FLT_NO_BACK = 1 << 2
    FLT_NO_SHORT_META = 1 << 3    # space, enter, arrows, etc.
    FLT_NO_LONG_META = 1 << 4 # shift, control, alt, ect.
    FLT_NO_PUNCT = 1 << 5
    
    def __init__(self):
        self.dataKeys = None
        self.dataHT = None
        self.dataTimeStart = None
        self.dataTimeEnd = None
        pass
    

    def sanityCheck( self ):
        """
        Filter out keystrokes variables in the member variables. 
        Eliminate anything < 0.
        returns the number of elements removed
        """
        assert( self.dataKeys is not None and len(self.dataKeys) > 0 )
        assert( self.dataHT is not None and len(self.dataHT) > 0 )
        assert( self.dataTimeStart is not None and len(self.dataTimeStart) > 0 )
        assert( self.dataTimeEnd is not None and len(self.dataTimeEnd) > 0 )
        
        badLbl = self.dataTimeStart <= 0
        badLbl = np.bitwise_or( badLbl,  self.dataTimeEnd <= 0)
        badLbl = np.bitwise_or( badLbl,  self.dataHT < 0)
        badLbl = np.bitwise_or( badLbl,  self.dataHT >= 5)
        #----- remove non consecutive start times
        nonConsTmpLbl = np.zeros( len(self.dataTimeStart) ) == 0 # start with all True labels
        nonConsLbl = np.zeros( len(self.dataTimeStart) ) > 0 # start with all False labels
        startTmpArr = self.dataTimeStart.copy()
        while ( np.sum( nonConsTmpLbl ) > 0 ):
			# find non consecutive labels
			nonConsTmpLbl = np.append([False], np.diff(startTmpArr)<0)                
   			# keep track of the indeces to remove
			nonConsLbl = np.bitwise_or( nonConsLbl,  nonConsTmpLbl)
               # changes value in the temporary array
			indecesToChange = np.arange(len(nonConsTmpLbl))[nonConsTmpLbl]
			startTmpArr[indecesToChange] = startTmpArr[indecesToChange-1]

        badLbl = np.bitwise_or( badLbl,  nonConsLbl)
        #-----
		
        # invert bad labels
        goodLbl = np.bitwise_not(badLbl)
        
        self.dataKeys = self.dataKeys[goodLbl]
        self.dataHT = self.dataHT[goodLbl]
        self.dataTimeStart = self.dataTimeStart[goodLbl]
        self.dataTimeEnd = self.dataTimeEnd[goodLbl]
             
        
        return sum(badLbl)
    
    def loadDataFile(self, fileIn, autoFilt=True, impType=None, debug=False):  
        """
        Load raw data file
        """      
        errorStr = ''
        try:
            data = []
            
#            if data.dtype == np.int64: # Sleep inertia format
            if impType =='si':
                data = np.genfromtxt(fileIn, dtype=long, delimiter=',', skip_header=0)
                data = data - data.min()
                data = data.astype(np.float64) / 1000
                self.dataTimeStart = data[:,0]  
                self.dataTimeEnd = data[:,1]
                self.dataHT = self.dataTimeEnd - self.dataTimeStart
                #TO REMOVE
                self.dataKeys = np.zeros(len(self.dataHT))#Just to make sanity work
                remNum = self.sanityCheck()
                #print remNum
            else: # PD format
                data = np.genfromtxt(fileIn, dtype=None, delimiter=',', skip_header=0)
                # load
                self.dataKeys = data['f0']
                self.dataHT = data['f1']  
                self.dataTimeStart = data['f3']  #No CHANGED 2<->3
                self.dataTimeEnd = data['f2']
                remNum = self.sanityCheck()
                #print '{:}, {:} %'.format( remNum, 1.0*remNum/len(self.dataHT) )
                
                if (debug):
                    print 'removed ', str(remNum), ' elements'

                if( autoFilt ):
                    self.filtData(self.FLT_NO_MOUSE  | self.FLT_NO_LONG_META )
            
            # load flight time
            self.dataFT = np.array([ self.dataTimeStart[i]-self.dataTimeStart[i-1]  for i in range(1,self.dataTimeStart.size) ])
            self.dataFT = np.append(self.dataFT, 0)
            
            
            
            return True
        except IOError:
            errorStr = 'file {:s} not found'.format(fileIn)
            return errorStr
    def loadDataArr(self, lstArr):
        self.dataKeys = np.zeros((len(lstArr),1), dtype='S30')
        self.dataHT = np.zeros((len(lstArr),1))
        self.dataTimeStart = np.zeros((len(lstArr),1))  
        self.dataTimeEnd =np.zeros((len(lstArr),1))
        i = 0
        for row in lstArr:
            tok = row.split(',')
            self.dataKeys[i] = str(tok[0])
            self.dataHT[i] = str(tok[1])
            self.dataTimeStart[i] = str(tok[2])
            self.dataTimeEnd[i] = str(tok[3]) 
            i += 1
            
        #self.loadDataFile(lstArr.toString())
    

    def filtData(self, flags):
        """
        Filter data
        return (fltKeys, fltHT, fltTimeStart, fltTimeEnd)
        """
        #-- filters
        pMouse=re.compile('("mouse.+")')
        pChar=re.compile('(".{1}")')
        pBack=re.compile('("BackSpace")')
        pLongMeta=re.compile('("Shift.+")|("Alt.+")|("Control.+")')
        pShortMeta=re.compile('("space")|("Num_Lock")|("Return")|("P_Enter")|("Caps_Lock")|("Left")|("Right")|("Up")|("Down")')
        pPunct=re.compile('("more")|("less")|("exclamdown")|("comma")|("\[65027\]")|("\[65105\]")|("ntilde")|("minus")|("equal")|("bracketleft")|("bracketright")|("semicolon")|("backslash")|("apostrophe")|("comma")|("period")|("slash")|("grave")')
        #--

        #-- create mask labels        
        lbl = np.ones(len( self.dataKeys ))==1
        if( flags & self.FLT_NO_MOUSE ):
            lblTmp = [ pMouse.match( k ) is None for k in self.dataKeys]
            lbl = lbl & lblTmp
        if( flags & self.FLT_NO_LETTERS ):
            lblTmp = [ pChar.match( k ) is None for k in self.dataKeys]
            lbl = lbl & lblTmp
        if( flags & self.FLT_NO_BACK ):
            lblTmp = [ pBack.match( k ) is None for k in self.dataKeys]
            lbl = lbl & lblTmp
        if( flags & self.FLT_NO_SHORT_META ):
            lblTmp = [ pShortMeta.match( k ) is None for k in self.dataKeys]
            lbl = lbl & lblTmp
        if( flags & self.FLT_NO_LONG_META ):
            lblTmp = [ pLongMeta.match( k ) is None for k in self.dataKeys]
            lbl = lbl & lblTmp
        if( flags & self.FLT_NO_PUNCT ):
            lblTmp = [ pPunct.match( k ) is None for k in self.dataKeys]
            lbl = lbl & lblTmp
        #--
        
        self.lbl = lbl        
        
        self.dataKeys = self.dataKeys[lbl]
        self.dataHT = self.dataHT[lbl]
        self.dataTimeStart = self.dataTimeStart[lbl]
        self.dataTimeEnd = self.dataTimeEnd[lbl]        
        
    def getStdVariablesFilt( fileIn, impType=None ):
        """
        Receives as parameter the location of the raw typing file
        Return filtered variables (i.e. no mouse clicks, no long meta buttons, no backspaces) 
        format returned (array of keys, array of hold times, array of press events timestamps, array of release events timestamps )
        """
        nqObj = self
        res = nqObj.loadDataFile( fileIn, False, impType)
        # remove delete button
        nqObj.filtData(nqObj.FLT_NO_MOUSE  | nqObj.FLT_NO_LONG_META | nqObj.FLT_NO_BACK )
        assert(res==True) # make sure the file exists
        dataKeys = nqObj.dataKeys
        dataHT = nqObj.dataHT
        dataTimeStart = nqObj.dataTimeStart
        dataTimeEnd = nqObj.dataTimeEnd
        
        return dataKeys, dataHT, dataTimeStart, dataTimeEnd


def getDataFiltHelper( fileIn, impType=None ):
    """
    Helper method to load filtered keypress data from given file
    :param fileIn: path to csv keypress file 
    :param impType: format of the csv file ('si': for sleep inertia data, None for PD data)
    :return: list of array with dataKeys, dataHT, dataTimeStart, dataTimeEnd
    """
    nqObj = NqDataLoader()
    res = nqObj.loadDataFile( fileIn, False, impType)
    # remove delete button
    nqObj.filtData(nqObj.FLT_NO_MOUSE  | nqObj.FLT_NO_LONG_META | nqObj.FLT_NO_BACK )
    assert(res==True) # make sure the file exists
    dataKeys = nqObj.dataKeys
    dataHT = nqObj.dataHT
    dataTimeStart = nqObj.dataTimeStart
    dataTimeEnd = nqObj.dataTimeEnd
    
    return dataKeys, dataHT, dataTimeStart, dataTimeEnd
    
    
def genFileStruct( dataDir, maxRepNum=4 ):
    '''
    Generate a dictionary with the NQ file list and test date (legacy method)
    :param dataDir: base directory containing the CSV files
    :param maxRepNum: integer with the maximum repetition number
    :return: two dictionaries: fMap, dateMap = NQ file/date list[pID][repID][expID]
    '''
    fMap = {} # data container
    dateMap = {}
    files = os.listdir( dataDir )    
    p = re.compile( '([0-9]+)\.{1}([0-9]+)_([0-9]+)_([0-9]+)\.csv' )
    for f in files:
        m = p.match( f )
        
        if( m ): # file found
            timeStamp = m.group(1)
            pID = int(m.group(2))
            repID = int(m.group(3))
            expID = int(m.group(4))
            # store new patient
            if( not fMap.has_key(pID) ):
                fMap[pID] = {}
                dateMap[pID] = {}
                for tmpRid in range(1, maxRepNum+1):
                    fMap[pID][tmpRid] = {}
                    dateMap[pID][tmpRid] = {}
                # fMap[pID] = {1: {}, 2: {}, 3: {}, 4:{}}
            # store data
            fMap[pID][repID][expID] = dataDir + f
            dateMap[pID][repID][expID] = datetime.datetime.fromtimestamp(int(timeStamp))
        else:
            print f, ' no'
            
    return fMap, dateMap
