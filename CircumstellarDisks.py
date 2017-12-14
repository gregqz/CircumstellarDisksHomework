import argparse
import time
import string
import os
import re
import StringIO
import cookielib
import urllib
from urllib2 import HTTPError
import mechanize
import logging
import math
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import random
import sys
from bs4 import BeautifulSoup
from urlparse import urlsplit
import requests
import itertools
from urlparse import urlparse
import posixpath
import ntpath
import ssl
from socket import error as SocketError
import errno

debugView = False
debugFileBool = False
debugLogFile = ""
outputLogFile = ""


def path_parse( path_string, normalize = True, module = posixpath ):
    result = []
    if normalize:
        tmp = module.normpath( path_string )
    else:
        tmp = path_string
    while tmp != "/":
        ( tmp, item ) = module.split( tmp )
        result.insert( 0, item )
    return result

def dump_array( array ):
    string = "[ "
    for index, item in enumerate( array ):
        if index > 0:
            string += ", "
        string += "\"{}\"".format( item )
    string += " ]"
    return string

def test_url( url, normalize = True, module = posixpath ):
    url_parsed = urlparse( url )
    path_parsed = path_parse( urllib.unquote( url_parsed.path ),
        normalize=normalize, module=module )
    return path_parsed
    sys.stdout.write( "{}\n  --[n={},m={}]-->\n    {}\n".format( 
        url, normalize, module.__name__, dump_array( path_parsed ) ) )



def testErrorCodes(br,theCodes):
    for x in theCodes:

        br.select_form(nr=0)

        theAction = br.action
        payload = {'code': x}

        response = requests.post(theAction, data=payload)
        print response.status_code

def parse_args():
    parser = argparse.ArgumentParser(description="Download Disk Data")
    parser.add_argument('-h', help="Shows this")   
    return parser.parse_args()

def debug_print(stringError=""):
    print stringError
    r = str(raw_input('Enter to Continue   '))
    if r == "q" : sys.exit()

def debugOutput(stringOutput=""):
    global debugView
    global debugFileBool
    global debugLogFile
    if debugView and debugFileBool: 
        debugLogFile.write(stringOutput) 
    elif debugView:
        print stringOutput
    else:
        i = 1
    return 0

globalCookieJar = "";
def createBrowserObject():
    global globalCookieJar
    #initalize the Cookie Jar and Browser
    br = mechanize.Browser()
    #print globalCookieJar
    if globalCookieJar == "":
        globalCookieJar = cookielib.LWPCookieJar()
    
    br.set_cookiejar(globalCookieJar)
    
    br.set_handle_equiv(True)
    listOfUserAgents = 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/42.0.2311.90 Safari/537.36]]Mozilla/5.0 (Windows NT 6.1; WOW64; rv:37.0) Gecko/20100101 Firefox/37.0]]Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_3) AppleWebKit/600.5.17 (KHTML, like Gecko) Version/8.0.5 Safari/600.5.17]]Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2272.118 Safari/537.36]]Mozilla/5.0 (Windows NT 6.3; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/42.0.2311.90 Safari/537.36]]Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/42.0.2311.135 Safari/537.36]]Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/42.0.2311.90 Safari/537.36]]Mozilla/5.0 (Windows NT 6.3; WOW64; rv:37.0) Gecko/20100101 Firefox/37.0]]Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1'.split("]]")
    br.addheaders=[('user-agent',random.choice(listOfUserAgents))]
    br.set_handle_gzip(True)
    br.set_handle_equiv(True)
    br.set_handle_redirect(True)
    br.set_handle_referer(True)
    br.set_handle_robots(False)
    #br.set_debug_http(True)
    #br.set_debug_redirects(True)
    #br.set_debug_responses(True)
    br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
    
    return br

class galaxy:
    class RA:
        def __init__(self,string):
            string = string.replace("  "," ")
            self.hr = float(string.split(" ")[0])
            self.min = float(string.split(" ")[1])
            self.sec = float(string.split(" ")[2])
            
        def __str__(self):
            return "" + str(self.hr) + "h " + str(self.min) + "m " + str(self.sec) + "s"
        
        def toradians(self):
            return math.pi/12 * self.hr + math.pi/720 * self.min + math.pi/43200 * self.sec
        
        def todegrees(self):
            return 180/math.pi * self.toradians()
        
    class Dec:
        def __init__(self,string):
            string = string.replace("  "," ")
            self.deg = float(string.split(" ")[0])
            self.arcmin = float(string.split(" ")[1])
            self.arcsec = float(string.split(" ")[2])
            
        def __str__(self):
            return "" + str(self.deg) + " deg " + str(self.arcmin) + " ' " + str(self.arcsec) + " ''"
                
        def toradians(self):
            return math.pi/180 * self.todegrees()
        
        def todegrees(self):
            return self.deg + self.arcmin/60 + self.arcsec/3600
    
    class angmomvect:
        def __init__(self):
            self.start = []
            self.end = []
            self.flag = 0
        
        def setstartend(self,startArr,endArr):
            self.start = startArr
            self.end = endArr
            self.flag = 1
        
        def getflag(self):
            return self.flag
        
        def getstartend(self):
            return [item for sublist in [self.start,self.end] for item in sublist]
        
        def __str__(self):
            return "Starting Point : " + dump_array(self.start) + "; Ending Point : " + dump_array(self.end)
    
    def __init__(self,galname,galRA,galDec,galinc,galpa,galDis):
        galpa = galpa.replace(u'\u2212', '-').replace(u'\u2013', '-')
        if galpa == 'N/A':
            galpa = 'nan'
        self.name = galname #Name
        self.RA = galaxy.RA(galRA)
        self.dec = galaxy.Dec(galDec)
        self.inc = float(galinc)
        self.dis = float(galDis)
        self.pa = float(galpa)
        self.gall = float('nan')
        self.galb = float('nan')
        self.galx = float('nan')
        self.galy = float('nan')
        self.galz = float('nan')
        self.angmoment = galaxy.angmomvect()

    def __str__(self):
        if math.isnan(self.galx) or math.isnan(self.galy) or math.isnan(self.galz):
            return "Name : " + self.name + ", RA : " + str(self.RA) +", Dec : " + str(self.dec) + ", Inclination : " + str(self.inc) + ", Distance : "  + str(self.dis)
        else:
            return "Name : " + self.name + ", RA : " + str(self.RA) +", Dec : " + str(self.dec) + ", Inclination : " + str(self.inc) + ", Distance : "  + str(self.dis) + ", (x,y,z) :" + str(self.getcartcoord()) + "\r\n"
    
    def calcgalcoord(self):
        a0 = 192.8595 * math.pi/180
        d0 = 27.1287 * math.pi/180
        l0 = 122.9320 * math.pi/180
        a = self.RA.toradians()
        d = self.dec.toradians()
        self.gall = (l0 - math.atan((math.cos(d)*math.sin(a-a0))/(math.sin(d)*math.cos(d0)-math.cos(d)*math.sin(d0)*math.cos(a-a0))))* 180/math.pi
        self.galb = math.asin(math.sin(d)*math.sin(d0) + math.cos(d)*math.cos(d0)*math.cos(a-a0))*180/math.pi
    
    def calccartcoord(self):
        theta = self.gall * math.pi/180
        phi = (90 - self.galb) * math.pi/180
        r = self.dis
        self.galx = r*math.sin(phi)*math.cos(theta)
        self.galy = r*math.sin(phi)*math.sin(theta)
        self.galz = r*math.cos(phi)
    
    def getgalcoord(self):
        if math.isnan(self.gall) or math.isnan(self.galb):
            self.calcgalcoord()
        return (self.gall, self.galb)

    def getcartcoord(self):
        if math.isnan(self.galx) or math.isnan(self.galy) or math.isnan(self.galz):
            self.getgalcoord()
            self.calccartcoord()
        return (self.galx, self.galy, self.galz)
    
    def getangularmomentvect(self, lenr = 1):
        if self.angmoment.getflag() == 0:
            x0,y0,z0 = self.getcartcoord()
            if True in map(math.isnan,[self.inc,self.pa]):
                self.angmoment.setstartend([float('nan'),float('nan'),float('nan')],[float('nan'),float('nan'),float('nan')])
                return (self.angmoment) 
            if self.inc == 90:
                self.pa = self.pa - 90
            print " inc " + str(self.inc) + " ; pa " + str(self.pa) + " ; gall " + str(self.gall) + " ; galb " + str(self.galb)
            self.galb *= math.pi/180
            self.gall *= math.pi/180
            self.inc *= math.pi/180
            self.pa *= math.pi/180
            print " inc " + str(self.inc) + " ; pa " + str(self.pa) + " ; gall " + str(self.gall) + " ; galb " + str(self.galb)
            angx=(math.cos(self.galb) * math.cos(self.gall)**(2) +math.sin(self.gall)**(2) )*(math.cos(self.inc) * math.cos(self.gall)-math.sin(self.inc) * math.sin(self.pa) * math.sin(self.gall))+(-math.cos(self.galb) * math.sin(self.gall) * math.cos(self.gall)-math.sin(self.gall) * math.cos(self.gall))*(math.cos(self.inc) * math.sin(self.gall)**(2) +math.sin(self.inc) * math.sin(self.pa) * math.cos(self.gall))+(math.sin(self.galb) * math.sin(self.gall))*(math.sin(self.inc) * math.cos(self.pa))
            angy=(math.cos(self.galb) * math.cos(self.gall) * math.sin(self.gall)-math.cos(self.gall) * math.sin(self.gall))*(math.cos(self.inc) * math.cos(self.gall)-math.sin(self.inc) * math.sin(self.pa) * math.sin(self.gall))+(math.cos(self.galb) * math.sin(self.gall)**(2) +math.cos(self.gall)**(2) )*(math.cos(self.inc) * math.sin(self.gall)**(2) +math.sin(self.inc) * math.sin(self.pa) * math.cos(self.gall))+(-math.sin(self.galb) * math.sin(self.gall))*(math.sin(self.inc) * math.sin(self.pa))
            angz = math.sin(self.galb) * math.cos(self.inc) * (math.cos(self.gall)**(2) +math.sin(self.gall)**(3) )+math.cos(self.galb) * math.sin(self.inc) * math.cos(self.pa)
            self.galb *= 180/math.pi
            self.gall *= 180/math.pi
            self.inc *= 180/math.pi
            self.pa *= 180/math.pi            start = [x0,y0,z0]
            end = [x0+angx,y0+angy,z0+angx]
            self.angmoment.setstartend(start,end)
        return (self.angmoment)

def findString(obj):
    for descendant in obj.descendants:
        if isinstance(descendant, basestring):
            return descendant


def gatherData():
    global debugView
    
    br = createBrowserObject()

 
    #open up the webpage
    rMain = br.open("https://webdisks.jpl.nasa.gov/")
    
    
    #parse the webpage using BeautifulSoup and the lxml parser
    mainHtml = rMain.read()
    mainSoup = BeautifulSoup(mainHtml,"lxml")
    print "Got Webpage and Initalized Everything"
    listOfLinks = mainSoup.find('table',{"id":"object_table"}).find_all('a')
    listofgals = []
    n = 0
    for link in listOfLinks:
        url = 'http://webdisks.jpl.nasa.gov/' + link.get('href')
        url = url.split("&")[0]
        print url
        rGal = br.open(url)
        galHTML = rGal.read()
        galSoup = BeautifulSoup(galHTML,"lxml")
        galtable = galSoup.find('table')
        tableth = map(findString, galtable.find_all("th"))
        tabletd = map(findString, galtable.find_all("td"))
        m = [tableth,tabletd]
        t_m = zip(*m)
        tempra = 'nan nan nan'
        tempdec = 'nan nan nan'
        tempdis = 'nan'
        tempinc = 'nan'
        temppa = 'nan'
        for row in t_m:
            if "RA" in row[0]:
                tempra = row[1]
            if "DEC" in row[0]:
                tempdec = row[1]
            if "Distance" in row[0]:
                tempdis = row[1]
            if "Inclination" in row[0]:
                tempinc = row[1]
            if "PA" in row[0]:
                temppa = row[1]
        print tempra + " " + tempdec + " " + tempdis + " " + tempinc + " " + temppa
        thisgal = galaxy(link.string,tempra,tempdec,tempinc,temppa,tempdis)
        listofgals.append(thisgal)
        n = n + 1
        if n > 10:
            #break
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    listofcoords = map(galaxy.getangularmomentvect,listofgals)
    x = []
    y = []
    z = []
    u = []
    v = []
    w = []
    for coord in listofcoords:
        print coord
        if not float('nan') in coord.getstartend():
            x.append(coord.start[0])
            y.append(coord.start[1])
            z.append(coord.start[2])
            u.append(coord.end[0])
            v.append(coord.end[1])
            w.append(coord.end[2])
    
    with open('Coordinates.txt', 'a') as the_file:
        for gal in listofgals:
            the_file.write(str(gal))
    
    ax.scatter(x, y, z, c='r', marker='o')
    ax.quiver(x, y, z, u, v, w)
    
    ax.set_xlabel('X (parsecs)')
    ax.set_ylabel('Y (parsecs)')
    ax.set_zlabel('Z (parsecs)')
    ax.set_xlim([-500, 500])
    ax.set_ylim([-500, 500])
    #ax.set_zlim([-200, 200])
    for ii in xrange(0,360,1):
        ax.view_init(elev=10., azim=ii)
        plt.savefig("./movie/movie%d.png" % ii)
    
    br.close()
    
def main():
    global debugView
    global debugFileBool
    global debugLogFile
    global outputLogFile
    reload(sys)  
    sys.setdefaultencoding('utf8')
    debugLogFile = open("./debug_log.txt","w")
    outputLogFile = open("./output_log.txt","w")
    gatherData()
    debugLogFile.close()
    outputLogFile.close()
    print "All Done!"

if __name__ == "__main__":
    main()
